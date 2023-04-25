;;; mcmc-fitting.lisp
#|
create walker
advance it
visualize it
|#

;;; Example use
;; (defparameter woi (mcmc-fit :function (lambda (x &key m b &allow-other-keys) (+ b (* -3 m) (* (- m (/ b 60)) x))) :data '((-4 -1 2 5 10) (0 2 5 9 13)) :params '(:b -1 :m 2) :data-error 0.2))

;; (defparameter woi (mcmc-fit :function (list (lambda (x &key b m c d &allow-other-keys) (+ b (* m x) (* c x x) (* d x x x))) (lambda (x &key e m g &allow-other-keys) (+ e (* (+ m g) x)))) :data '(((0 1 2 3 4) (4 5 6 8 4)) ((10 20 30 40 50) (1 2 3 4 5))) :params '(:b -1 :m 2 :c 1 :d -1 :e 0.5 :g -2) :data-error 0.2))


;; (sb-ext:muffle-conditions sb-ext:compiler-note)

;; Full arrays > list of arrays > pure lists
;; just lists don't benefit from optimization
;; arrays do a lot
;; arrays aren't much better than lists without optimization

;; (ql:quickload :gsll)
;; (ql:quickload :antik)

#| For testing cholesky decomps vs gsl implementation
(defun eyes (n)
(do ((i 0 (1+ i))
(ret-array nil (push (do ((j 0 (1+ j))
(ret-list nil (push (if (= i j) 1 0) ret-list)))
((> j n) (reverse ret-list)))
ret-array)))
((> i n) (reverse ret-array))))

(defun matrix-add (mat1 mat2)
(let ((ret-mat nil))
(dolist (erc (mapcar #'list mat1 (transpose-matrix mat2)) (reverse ret-mat))
(push (mapcar #'+ (elt erc 0) (elt erc 1)) ret-mat))))

(defparameter randmat (repeat 3 (repeat 3 (random 10d0))))
(defparameter posdefmat (matrix-add (scale-matrix 10 (eyes 3)) (scale-matrix 0.5 (matrix-add randmat randmat))))
|#

(in-package :mcmc-fitting)


;;; Utility
(defun elts (tree &rest elts)
  "Apply multiple elt to a tree. I.e. (elt (elt (elt foo 3) 2) 1) is simply (elts foo 3 2 1) or (apply #'elts foo '(3 2 1))."
  (reduce (lambda (x y) (elt x y)) elts :initial-value tree))

(defun (setf elts) (new-value tree &rest elts)
  "Make elts setf-able. Usage is (setf (elts tree elt0 elt1 elt2 ...) new-value) or (setf (apply #'elts tree eltslist) new-value)."
  (setf (elt (apply #'elts tree (butlast elts)) (car (last elts))) (if (listp new-value) (copy-seq new-value) new-value)))

(defmacro return-this-but-also (expr &body body)
  "Returns the provided 'expr' while also performing manipulations on 'expr' provided in 'body' to allow for printing subsections, performing logic, etc. Uses 'it' as the key term in 'body' to refer to the result of the 'expr'.
i.e. (return-this-but-also (list 4 8 2 0 4 12 0) (print (list (count 0 it) (count-if #'evenp it) (position 8 it))))
==> (2 7 1)
==> (4 8 2 0 4 12 0)"
  (let ((it-g (gensym)))
    `(let ((,it-g (multiple-value-list ,expr)))
       (let ((it (lambda () (values-list ,it-g))))
	 ,@(subst '(funcall it) 'it body)
	 (funcall it)))))

(defun range (start-or-end &optional end)
  "Provides a list of numbers from 'start' to 'end' incrementing by 1. Always integers. Floors inputs.
(range 5.2) => '(0 1 2 3 4)
(range 2 5) => '(2 3 4)
(range 5 2) => '(5 4 3)"
  (let ((start (floor (if end start-or-end 0)))
	(end (floor (if end end start-or-end))))
    (let ((return-list nil) (sign (signum (- end start))))
      (do ((curr start (progn (push curr return-list) (+ sign curr))))
	  ((= curr end) (reverse return-list))))))

(defun slice (2d-list &optional rows cols)
  "Take a slice of a 2d list taking 'rows' (i.e. '(1 3 5) for rows 1 then 3 then 5) and 'cols' (i.e. '(20 40) for items 20 through 40)."
  (let ((rows (if rows rows (range (length 2d-list))))
	(cols (if cols cols (list 0 (length (elt 2d-list 0))))))
    (mapcar (lambda (x) (apply #'subseq (elt 2d-list x) cols)) rows)))

(defmacro mapcar-enum ((element index list) &body body)
  "mapcar with enumeration. Iterates over 'list' using 'element' to store each element, 'index' to store the index, and \"list\" to store the full list. Returns a list made of applying 'body' to each element of 'list'.
try: (mapcar-enum (e i '(30 20 10)) (print (list e i)))"
  (let ((return-list-g (gensym "RETURN-LIST")))
    `(let ((,return-list-g nil)
	   (list ,list))
       (do ((,index 0 (1+ ,index))
	    (,element nil))
	   ((>= ,index (length list)))
	 (setf ,element (elt list ,index))
	 (push ,@body ,return-list-g))
       (reverse ,return-list-g))))

(defun mapcar-tree (function tree)
  "Maps 'function' over the deepest elements of 'tree'. 'tree' can be composed of any combination of sequences.
E.g. (mapcar-tree #'1+ '((1 2) (2 3) (3 4) (100 (12 23 (324 54)))))
 =>  '((2 3) (3 4) (4 5) (101 (13 24 (325 55))))."
  (cond ((null tree)
	 nil)
	((and (atom tree) (or (stringp tree) (not (vectorp tree))))
	 (funcall function tree))
	(t
	 (map (type-of tree) (lambda (x) (mapcar-tree function x)) tree))))
(export 'mapcar-tree)

(defun plist-keys (plist)
  "Get all unique keys from 'plist'."
  (do ((keys (list (pop plist)) (progn (pop plist) (append (list (pop plist)) keys))))
      ((null (nthcdr 2 plist)) (reverse (remove-duplicates keys)))))

(defun plist-values (plist)
  "Get all values belonging to unique keys in 'plist'. Returns the first example of a value found."
  (let ((keys (plist-keys plist)))
    (mapcar (lambda (x) (getf plist x)) keys)))

(defun make-plist (keys values)
  "Riffles a set of 'keys' and 'values' into a plist."
  (apply #'append (mapcar #'list keys values)))

(defmacro repeat (count &body expression)
  `(mapcar (lambda (x) (declare (ignore x)) ,@expression) (make-list ,count)))

(defun mkstr (&rest args)
  (with-output-to-string (s)
    (dolist (a args) (princ a s))))

(defun symb (&rest args)
  (values (intern (apply #'mkstr args))))

(defun symb-keyword (&rest args)
  (values (intern (apply #'mkstr args) :keyword)))

(defun remove-consecutive-duplicates (list)
  (mapcon (lambda (x)
	    (if (eql (car x) (cadr x)) nil (list (car x))))
	  list))

(defun remove-duplicates-plist (plist)
  (let* ((keys (get-plist-keys plist))
	 (unique-keys (remove-duplicates keys)))
    (make-plist
     unique-keys
     (mapcar (lambda (x) (getf plist x)) unique-keys))))

(defun elements-between (sequence start-value end-value)
  (remove-if-not (lambda (x) (<= start-value x end-value)) sequence))

(defun linspace (start end &key (len 50) (step nil) (type 'double-float))
  "Provides a list from 'start' to 'end' (inclusive) that is 'len' in length. Or provides a list of number from 'start' until 'end' by 'step'. 'end' may not be reached if 'step' kwarg is provided.
Supported types for 'type' arg are 'float, 'double-float, 'rational, 'integer, 'bignum."
  (let ((step (if step step (/ (rational (- end start)) (1- len))))
	(len (if step (+ 1 (floor (- end start) step)) len))
	(return-list nil))
    (let ((rational-numbers
	    (do ((n 0 (1+ n))
		 (curr (rational start) (progn (push curr return-list) (+ curr step))))
		((>= n len) (values (reverse return-list) step)))))
      (case type
	((or float double-float rational) (mapcar (lambda (x) (coerce x type)) rational-numbers))
	(t (mapcar #'round rational-numbers))))))

(defun n-nums (n)
  (do ((m 0 (1+ m))
       (l nil (push m l)))
      ((= m n) (reverse l))
    (declare (fixnum m))))

(defun up-to (n &optional (start 0))
  (assert (>= n start))
  (do ((m start (1+ m))
       (l nil (push m l)))
      ((= m (1+ n)) (reverse l))
    (declare (fixnum m))))

(defun diff (list)
  (mapcon (lambda (x) (if (nthcdr 1 x) (list (- (cadr x) (car x))) nil)) list))

(defun diff-matrix (matrix)
  (mapcon (lambda (x) (if (nthcdr 1 x) (list (mapcar #'- (cadr x) (car x))) nil)) matrix))

(defun diff-lplist (lplist)
  (let* ((keys (get-plist-keys (elt lplist 0)))
	 (values (mapcar (lambda (x) (mapcar (lambda (y) (getf x y)) keys)) lplist)))
    (mapcar (lambda (x) (make-plist keys x)) (diff-matrix values))))

(defun partition (list n)
  "Pairs elements and removes any with full matches"
  (labels ((rec (l &optional (acc nil))
	     (if (nthcdr (- n 1) l)
		 (rec (subseq l n) (cons (subseq l 0 n) acc))
		 (cons l acc))))
    (reverse (cdr (rec list)))))

(defun flatten (list)
  "Remove all structure from list"
  (let ((return-list nil))
    (labels ((%flatten (list)
	       (cond ((null list)
		      nil)
		     ((consp list)
		      (mapcar #'%flatten list))
		     (t
		      (push list return-list)))))
      (%flatten list))
    (reverse return-list)))

(defun split-string (splitter string)
  "Split string by splitter, which is a char."
  (declare (optimize speed safety)
	   (character splitter))
  (labels ((inner-split (string &optional acc)
	     (declare (optimize speed safety)
		      (simple-string string)
		      (type (or null cons) acc))
	     (if (string/= string "")
		 (let ((pos (position splitter string :test #'char=)))
		   (if pos
		       (inner-split (subseq string (+ pos 1)) (cons (subseq string 0 pos) acc))
		       (inner-split "" (cons string acc))))
		 (reverse acc))))
    (inner-split string)))



;;; Log-liklihood and log-prior functions

(defun log-prior-flat (params data)
  (declare (ignore params data))
  0d0)

(defmacro prior-bounds-let ((&rest keys-low-high) &body body)
  (declare (sb-ext:muffle-conditions sb-ext:compiler-note))
  (flet ((expand-keys (key-low-high)
	   (let* ((key (car key-low-high))
		  (param-name (read-from-string (symbol-name key))))
	     `(,param-name (the double-float (getf ,(intern "PARAMS") ,key 0d0)))))
	 (expand-bounds (key-low-high)
	   (destructuring-bind (key low-expr high-expr) key-low-high
	     (let ((param-name (read-from-string (symbol-name key)))
		   (param-name-bound (read-from-string (concatenate 'string (symbol-name key) "-bound"))))
	       `(,param-name-bound (the double-float (if (< ,low-expr ,param-name ,high-expr)
							 0d0
							 (* -1d10 (- (exp (* (min (abs (- ,param-name ,high-expr)) (abs (- ,param-name ,low-expr))) 1d-5)) 1))))))))
	 (get-bound-name (key-low-high)
	   (destructuring-bind (key low-expr high-expr) key-low-high
	     (declare (ignore low-expr high-expr))
	     (let ((param-name-bound (read-from-string (concatenate 'string (symbol-name key) "-bound"))))
	       param-name-bound))))
    `(let* (,@(mapcar #'expand-keys keys-low-high)
	    ,@(mapcar #'expand-bounds keys-low-high)
	    (,(intern "BOUNDS-TOTAL") (+ ,@(mapcar #'get-bound-name keys-low-high))))
       ,@body)))

;; log-liklihood can depend on x, y, params, stddev (additional distribution parameters)
;; in this case, it is just log normal (need x and y and mu (func x) and sigma)
;; in other cases, it could be log poisson (need lambda (func y) and k (y)) (no x, no sigma)
(defun log-normal (x mu sigma)
  (declare (sb-ext:muffle-conditions sb-ext:compiler-note)
	   (optimize speed)
	   (double-float mu x)
	   (type (double-float 0d0 *) sigma))
  (+ (* -1/2 (log (* 2 pi))) (* -1 (log sigma)) (* -1/2 (expt (/ (- x mu) sigma) 2.0))))

(defun log-factorial (n)
  (reduce (lambda (x y) (+ x (log y))) (up-to n)))

(defun log-poisson (lambda k)
  (- (* k (log lambda)) lambda (log-factorial k)))

;; log liklihood is the function that determines how data is useful
;; It does the pattern matching and analysis on input variables
;; The liklihood defines the result
;; Each independent variable needs to have the appropriate distribution
;; Then they just all add together (for the log version)
;; The prior is then takked on at the end however the user sees fit
;; Since the prior does not change, it will just speed up the convergence to have a good one
;; prior could also depend on data, generally for bounds or such
(defun create-log-liklihood-normal-weighted (fn params data stddev data-column-x data-column-y)
  (declare (optimize speed)
	   (simple-vector data)
	   (function fn)
	   (simple-vector stddev))
  (let* ((x (elt data data-column-x))
	 (y (elt data data-column-y))
	 (stddev (if (= 1 (length stddev)) (make-array (length x) :initial-element (elt stddev 0)) stddev)))
    (declare (simple-vector x y stddev))
    (reduce #'+ (map 'vector (lambda (x y z) (log-normal y (apply fn x params) z)) x y stddev))))

(defun log-liklihood-normal (fn params data stddev)
  (create-log-liklihood-normal-weighted fn params data stddev 0 1))

(defun log-liklihood-normal-cutoff (fn params data stddev)
  (declare (optimize speed)
	   (simple-vector data)
	   (function fn)
	   (simple-vector stddev))
  (let* ((x (elt data 0))
	 (y (elt data 1))
	 (stddev (if (= 1 (length stddev)) (make-array (length x) :initial-element (elt stddev 0)) stddev)))
    (declare (simple-vector x y stddev))
    (reduce #'+ (map 'vector (lambda (x y z) (max -5000d0 (log-normal y (apply fn x params) z))) x y stddev))))
(export 'log-liklihood-normal-cutoff)



;;; plist functions
(defun get-plist-keys (plist &optional return-keys)
  (if (car plist)
      (get-plist-keys (cddr plist) (cons (car plist) return-keys))
      (reverse return-keys)))

(defun get-plist-values (plist &optional return-values)
  (if (car plist)
      (get-plist-values (cddr plist) (cons (cadr plist) return-values))
      (reverse return-values)))

(defun reduce-plists (function plist1 plist2 &optional result)
  (if (and (cadr plist1) (cadr plist2))
      (let ((key (car plist1)))
	(reduce-plists function (cddr plist1)
		       plist2
		       (nconc (list (funcall function (cadr plist1) (getf plist2 key)) key) result)))
      (reverse result)))

(defun map-plist (fn plist)
  (let* ((keys (get-plist-keys plist))
	 (values (get-plist-values plist))
	 (ret-values (mapcar fn values)))
    (make-plist keys ret-values)))

(defun scale-plist (scale plist)
  (map-plist (lambda (x) (* x scale)) plist))



;;; Walker and associated functions
(defstruct walker-step
  (prob most-negative-double-float :type float)
  (params nil :type list))
(export '(walker-step-prob walker-step-params))

(defstruct walker
  (function nil :type list)
  (param-keys nil :type list)
  (walk nil :type list)
  (length 1 :type integer)
  (most-likely-step (make-walker-step) :type walker-step)
  (last-step (make-walker-step) :type walker-step)
  (data nil :type (or list vector))
  (data-error nil :type list)
  (log-liklihood nil :type list)
  (log-prior nil :type list))
(export '(walker-function walker-param-keys walker-walk walker-length walker-most-likely-step walker-last-step walker-data walker-data-error walker-log-liklihood walker-log-prior))

;; Useful for inspecting when algorithm fails
(defun walker-check-for-complex-walks (walker &optional take)
  (let ((complex-walks-ps (mapcar (lambda (x) (some #'identity x)) (mapcar-tree #'complexp (lplist-to-l-matrix (diff-lplist (walker-get walker :get :unique-steps :take take)))))))
    (if (some #'identity complex-walks-ps) complex-walks-ps nil)))

(defun walker-get (walker &key (get (or :steps :unique-steps :most-likely-step :acceptance :param :most-likely-params :median-params :all-params :stddev-params :log-liklihoods :covariance-matrix :l-matrix)) take param)
  "Get any of the above items from a walker. Use the 'take' kwarg to limit the number of walks you are sampling over. Use the 'param' kwarg to specify which params for the :param :get option."
  (case get
    (:steps (subseq (walker-walk walker) 0 take))
    (:unique-steps (mapcon (lambda (x)
			     (if (equal (walker-step-prob (car x)) (when (cadr x) (walker-step-prob (cadr x))))
				 nil
				 (list (walker-step-params (car x)))))
			   (walker-get walker :get :steps :take take)))
    (:most-likely-step (reduce (lambda (x y)
				 (if (> (walker-step-prob x) (walker-step-prob y)) x y))
			       (walker-get walker :get :steps :take take)))
    (:acceptance
     (let* ((probs (mapcar #'walker-step-prob (walker-get walker :get :steps :take take))))
       (/ (length (remove-consecutive-duplicates probs)) (length probs))))
    (:param (mapcar (lambda (x) (getf (walker-step-params x) param)) (walker-get walker :get :steps :take take)))
    (:most-likely-params
     (let ((param-keys (walker-param-keys walker))
	   (most-likely-step (walker-most-likely-step walker)))
       (make-plist param-keys
		   (mapcar (lambda (key) (getf (walker-step-params most-likely-step) key)) param-keys))))
    (:median-params
     (let ((param-keys (walker-param-keys walker))
	   (steps (walker-get walker :get :steps :take take)))
       (make-plist param-keys
		   (mapcar (lambda (key) (median (mapcar (lambda (step) (getf (walker-step-params step) key)) steps))) param-keys))))
    ;; TODO add failsafe here if not enough valid steps to do cholesky
    (:stddev-params
     (let* ((param-keys (walker-param-keys walker)))
       (if (< (walker-length walker) 10)
	   (make-plist param-keys (make-list (length param-keys) :initial-element 0d0))
	   (let ((l-matrix (walker-get walker :get :l-matrix :take take)))
	     (values (make-plist param-keys
				 (mapcar (lambda (x)
					   (elt (elt l-matrix x) x))
					 (up-to (1- (length l-matrix)))))
		     l-matrix)))))
    (:log-liklihoods (mapcar #'walker-step-prob (walker-get walker :get :steps :take take)))
    (:covariance-matrix (lplist-covariance (walker-get walker :get :unique-steps :take take)))
    (:l-matrix (lplist-to-l-matrix (diff-lplist (walker-get walker :get :unique-steps :take take))))))
(export 'walker-get)


(defun walker-modify (walker &key (modify (or :add-step :add-walks :burn-walks :reset :delete)) new-step new-walks burn-number keep-number)
  (case modify
    (:add-step (progn (push new-step (walker-walk walker))
		      (setf (walker-last-step walker) new-step)
		      (setf (walker-length walker) (1+ (walker-length walker)))
		      (when (> (walker-step-prob new-step)
			       (walker-step-prob (walker-most-likely-step walker)))
			(setf (walker-most-likely-step walker) new-step))))
    (:add-walks (progn (nconc (reverse new-walks) (walker-walk walker))
		       (setf (walker-last-step walker) (last new-walks))
		       (setf (walker-length walker) (+ (length new-walks) (walker-length walker)))
		       (mapcar
			(lambda (new-step)
			  (when (> (walker-step-prob new-step)
				   (walker-step-prob (walker-most-likely-step walker)))
			    (setf (walker-most-likely-step walker) new-step)))
			new-walks)))
    (:burn-walks (progn (setf (walker-walk walker) (subseq (walker-walk walker) 0 (- (walker-length walker) burn-number)))
			(setf (walker-length walker) (- (walker-length walker) burn-number))))
    (:keep-walks (progn (setf (walker-walk walker) (subseq (walker-walk walker) 0 keep-number))
			(setf (walker-length walker) keep-number)))
    (:reset (progn (setf (walker-walk walker) (last (walker-walk walker)))
		   (setf (walker-last-step walker) (car (walker-walk walker)))
		   (setf (walker-length walker) 1)
		   walker))
    (:delete (progn (setf (walker-walk walker) nil)
		    (setf walker (make-walker))))))
(export 'walker-modify)

(declaim (inline cholesky-decomp))
(defun cholesky-decomp (covariance-matrix)
  "Generalization of the sqrt operation for postive definite matrices."
  (let* ((cov-array (make-array (list (length covariance-matrix) (length (car covariance-matrix))) :initial-contents covariance-matrix))
	 (l-array (make-array (array-dimensions cov-array) :initial-element 0d0))
	 (tmp-sum 0d0))
    (dotimes (i (array-dimension cov-array 0))
      (dotimes (k (+ i 1))
	(setq tmp-sum 0d0)
	(dotimes (j k)
	  (incf tmp-sum (* (aref l-array i j) (aref l-array k j))))
	(if (= i k)
	    (setf (aref l-array i k) (sqrt (- (aref cov-array i k) tmp-sum)))
	    (setf (aref l-array i k) (/ (- (aref cov-array i k) tmp-sum) (aref l-array k k))))))
    (array-matrix l-array)))

;;; matrix and covariance
;; An array is a lisp array. A matrix is a list of lists
(declaim (inline dot))
(defun dot (list1 list2)
  "Returns dot product of two lists."
  (reduce #'+ (mapcar #'* list1 list2)))

(defun transpose-matrix (matrix)
  (let ((num-columns (- (length (elt matrix 0)) 1)))
    (mapcar (lambda (el) (mapcar (lambda (row) (elt row el)) matrix)) (up-to num-columns))))

(defun scale-matrix (scale matrix)
  (mapcar (lambda (x) (mapcar (lambda (y) (* y scale)) x)) matrix))

(declaim (inline lplist-covariance))
(defun lplist-covariance (lplist)
  (let* ((n (length lplist))
	 (keys (get-plist-keys (elt lplist 0)))
	 (values (mapcar (lambda (y) (mapcar (lambda (x) (getf x y)) lplist)) keys))
	 (avg (mapcar (lambda (x) (/ (reduce #'+ x) n)) values))
	 (a (mapcar (lambda (x y)
		      (mapcar (lambda (z) (- z y))
			      x))
		    values avg))
	 (ns (make-list (length a) :initial-element (float n 0d0)))
	 (v (mapcar (lambda (x denom)
		      (mapcar (lambda (y)
				(/ (dot x y) denom))
			      a))
		    a ns)))
    v))


(defun lplist-to-l-matrix (lplist)
  (let* ((covariance (lplist-covariance lplist)))
    (cholesky-decomp covariance)))

(declaim (inline get-covariant-sample))
(defun get-covariant-sample (means l-matrix)
  (labels ((rotate-random-sample (sample l-matrix)
	     (mapcar (lambda (y) (dot y sample)) l-matrix)))
    (let ((samples (repeat (length l-matrix) (alexandria:gaussian-random))))
      (mapcar #'+ means (rotate-random-sample samples l-matrix)))))

(defun diagonal-covariance (list)
  (let ((len (length list))
	(i -1))
    (mapcar (lambda (x)
	      (declare (ignore x))
	      (incf i)
	      (maplist (lambda (y) (if (= (length y) (- len i)) (car y) 0)) list))
	    list)))

(defvar example-lplist '((:a 90 :b 60 :c 90)
			 (:a 90 :b 90 :c 30)
			 (:a 60 :b 60 :c 60)
			 (:a 60 :b 60 :c 90)
			 (:a 30 :b 30 :c 30)))
(defvar example-covariance (lplist-covariance example-lplist))
(defvar example-array (make-array '(3 3) :initial-contents example-covariance))
(defun array-matrix (array)
  (let ((dim0 (array-dimension array 0))
	(dim1 (array-dimension array 1))
	(i -1)
	(j -1))
    (repeat dim0
      (incf i)
      (setq j -1)
      (repeat dim1
	(incf j)
	(aref array i j)))))
(defvar example-matrix (array-matrix example-array))
(defvar example-l-matrix (lplist-to-l-matrix example-lplist))




;;; MCMC Metropolis-Hastings
(defun force-list (item)
  "Forces all non-list items into lists"
  (if (consp item)
      item
      (list item)))

(defun clean-data-error (stddev clean-data)
  "If you give one number, it uniformly applies it to all y-data. If you give a list/vector, it applies that to each fn."
  (labels ((get-y-data-length (data)
	     (mapcar (lambda (x) (length (elt x 1))) data)))
    (let ((y-lengths (get-y-data-length clean-data)))
      (cond ((numberp stddev)
	     (mapcar (lambda (x) (coerce (make-list x :initial-element stddev) 'vector)) y-lengths))
	    ((and (= (length stddev) (length clean-data)))
	     (mapcar-enum (y-length index y-lengths)
	       (if (= (length (elt stddev index)) y-length)
		   (coerce (elt stddev index) 'vector)
		   (coerce (make-list y-length :initial-element (elt stddev index)) 'vector))))
	    (t
	     (let ((stddev (list stddev)))
	       (mapcar-enum (y-length index y-lengths)
		 (if (= (length (elt stddev index)) y-length)
		     (coerce (elt stddev index) 'vector)
		     (coerce (make-list y-length :initial-element (elt stddev index)) 'vector)))))))))

(defun clean-data (data number-of-functions)
  "Forces data to be a list at top level of equal length to the number of functions, and everything beneath the top level to be vectors."
  (labels ((vectorize-everything (sequence)
	     (cond ((numberp (elt sequence 0))
		    (coerce sequence 'vector))
		   (t
		    (map 'vector #'vectorize-everything sequence)))))
    (cond ((= (length data) number-of-functions)
	   (map 'list #'vectorize-everything data))
	  (t
	   (map 'list #'vectorize-everything (list data))))))

(defun create-walker-data (data &rest columns)
  "Takes a larger dataset (contains x, y, stddev, freq, temp, etc.) and extracts just the columns you want into a walker-friendly format."
  (mapcar (lambda (x) (apply #'vector x))
	  (mapcar (lambda (y) (elt data y)) columns)))
(export 'create-walker-data)

(defun to-double-floats (l)
  "Converts all numbers and vectors in tree to double-float type"
  (cond ((numberp l)
	 (coerce l 'double-float))
	((symbolp l) ; for p-lists
	 l)
	((consp l)
	 (mapcar #'to-double-floats l))
	((vectorp l)
	 (map 'vector #'to-double-floats l))))

(defun log-prior-fixer (log-prior params data)
  "Checks if your prior returns a function (i.e. changes shape based on the data provided to it). If so, gets that function. If not, returns the prior you gave it."
  (let ((results (map 'list (lambda (lp d) (funcall lp params d)) log-prior data)))
    (mapcar (lambda (res fn) (if (numberp res) fn res)) results log-prior)))

(defun log-liklihood-fixer (log-liklihood fn params data error)
  "Checks if your liklihood returns a function (i.e. changes shape based on the data provided to it). If so, gets that function. If not, returns the liklihood you gave it."
  (let ((results (map 'list (lambda (ll f d e) (funcall ll f params d e)) log-liklihood fn data error)))
    (mapcar (lambda (res fn) (if (numberp res) fn res)) results log-liklihood)))



(defun walker-many-steps (the-walker n &optional l-matrix)
  (if (null l-matrix) (print (car (setq l-matrix (diagonal-covariance (get-plist-values (scale-plist 1e-2 (walker-get the-walker :get :median-params))))))))
  (dotimes (i n)
    (walker-take-step the-walker :l-matrix l-matrix)))

;; TODO remove and add functionality for troubleshooting
;; (print (car (diagonal-covariance (get-plist-values (scale-plist 1e-2 (walker-get-median-params woi 1000))))))
;; (print (standard-deviation (diff (remove-consecutive-duplicates (walker-get-param woi :w1-0-0)))))
;; (plt:plot (diff (remove-consecutive-duplicates (walker-get-param woi :w1-0-0))))


(defun walker-adaptive-steps-full (walker &optional (n 100000) (temperature 1d3) (auto t) l-matrix)
  (let* ((param-keys (walker-param-keys walker))
	 (shutting-down-p nil)
	 (temp-steps 10000)
	 (temperature (rational temperature))
	 (temps (nconc (linspace 1.0 temperature :len (floor temp-steps 2)) (linspace temperature 1.0 :len (floor temp-steps 2))))
	 (current-acceptance 0))
    (flet ((stable-probs-p (probs)
	     (< 5 (- (reduce #'max probs) (reduce #'min probs)) 15))
	   (get-optimal-mcmc-l-matrix (take)
	     (scale-matrix (/ (expt 2.38 2) (length param-keys))
			   (walker-get walker :get :l-matrix :take take))))
      (unless l-matrix
	(if (or (< (walker-length walker) 2000)
		(< (walker-get walker :get :acceptance :take 100) 0.1))
	    (setq l-matrix (diagonal-covariance (get-plist-values (scale-plist 1e-5 (walker-get walker :get :median-params)))))
	    (setq l-matrix (get-optimal-mcmc-l-matrix 2000))))
      (do ((i 0 (+ i 1))
	   (temp-index 0)) 
	  ((>= i n))
	(when (or
	       (and (not shutting-down-p) (< (- n i) 30000))
	       (and auto
		    (not shutting-down-p)
		    (= 0 (mod i 1000))
		    (> i 40000)
		    (< 0.2 (walker-get walker :get :acceptance :take 1000) 0.5)
		    (stable-probs-p (walker-get walker :get :log-liklihoods :take 2000))))
	  (setf temperature 1)
	  (setf shutting-down-p t)
	  (setf i (- n 10000)))
	;; catch complex stddev and retry with differnt opimal l matrix sampling
	(walker-take-step walker :l-matrix l-matrix :temperature temperature)
	;; Annealing
	(when (and (not shutting-down-p) (= 0 (mod (1+ (floor i temp-steps)) 2)) (< i (* temp-steps 95)))
	  (setf temperature (elt temps temp-index))
	  (when (>= temp-index (1- temp-steps)) (walker-modify walker :modify :add-step :new-step (walker-get walker :get :most-likely-step :take (* 2 temp-steps))))
	  (if (>= temp-index (1- temp-steps)) (setf temp-index 0) (incf temp-index)))
	;; Regular l-matrix updating
	(when (and (not shutting-down-p) (= 0 (mod i 1000)) (> i 2000))
	  (setf current-acceptance (walker-get walker :get :acceptance :take 1000))
	  (if (< current-acceptance 0.2)
	      (setf l-matrix (scale-matrix 0.5 l-matrix))
	      (if (> current-acceptance 0.5)
		  (setf l-matrix (scale-matrix 2 l-matrix))
		  (setf l-matrix (get-optimal-mcmc-l-matrix 2000)))))))))

;; (let ((temp-steps 100)) (mapcar (lambda (i) (and (= 0 (mod (1+ (floor i temp-steps)) 2)) (< i (* temp-steps 95)))) (range 500)))

(defun walker-adaptive-steps (walker &optional (n 100000))
  (walker-adaptive-steps-full walker n 1d3 t nil))

;; (defun walker-construct-print-list (walker &optional take)
;;   `(:fn ,(mapcar #'sb-kernel:%fun-name (walker-function walker))
;;     :data ,(walker-data walker)
;;     :param-keys ,(walker-param-keys walker)
;;     :stddev ,(walker-data-error walker)
;;     :log-liklihood ,(mapcar #'sb-kernel:%fun-name (walker-log-liklihood walker))
;;     :log-prior ,(mapcar #'sb-kernel:%fun-name (walker-log-prior walker))
;;     :walker ,(walker-walks walker take)))

;; (defun walker-save (walker filename &optional take)
;;   (with-open-file (out filename
;;                        :direction :output
;;                        :if-exists :supersede)
;;     (with-standard-io-syntax
;;       (null (write (walker-construct-print-list walker take) :stream out)))))

;; (defun walker-load (filename &key function log-liklihood log-prior quiet)
;;   (with-open-file (in filename
;; 		      :direction :input)
;;     (with-standard-io-syntax
;;       (let* ((full (read in))
;; 	     (data (getf full :data))
;; 	     (data-error (getf full :data-error))
;; 	     (walks (getf full :walker))
;; 	     (dummy-params (walker-step-params (elt walks 0))))
;; 	(unless quiet
;; 	  (format t "*Recommendations*~%function: ~s~%log-liklihood: ~s~%log-prior: ~s~%" (getf full :fn) (getf full :log-liklihood) (getf full :log-prior)))
;; 	(when (and function log-liklihood log-prior)
;; 	  (let ((walker (walker-create :function function :data data :params dummy-params :data-error data-error :log-liklihood log-liklihood :log-prior log-prior)))
;; 	    (walker-modify walker :modify :add-walks :new-walks walks)
;; 	    walker))))))


;;; Functionality for a set of similar walker instances
;; (defun walker-set-save (the-walker-set filename &optional take)
;;   (let ((walker-print-list (mapcar (lambda (x) (walker-construct-print-list x take)) the-walker-set)))
;;     (with-open-file (out filename
;; 			 :direction :output
;; 			 :if-exists :supersede)
;;       (with-standard-io-syntax
;; 	(null (print walker-print-list out))))))

;; (defun walker-set-load (filename fn log-liklihood log-prior)
;;   (with-open-file (in filename
;; 		      :direction :input)
;;     (with-standard-io-syntax
;;       (let* ((full (read in)))
;; 	(mapcar (lambda (x)
;; 		    (let* ((param-keys (getf x :param-keys))
;; 			   (dummy-params (make-plist param-keys (make-list (length param-keys) :initial-element 0.0d0)))
;; 			   (data (getf x :data))
;; 			   (stddev (getf x :stddev))
;; 			   (walks (getf x :walker))
;; 			   (the-walker (walker-init :fn fn :data data :params dummy-params :stddev stddev :log-liklihood log-liklihood :log-prior log-prior :init nil)))
;; 		      (walker-modify the-walker :modify :add-walks :new-walks walks)
;; 		      the-walker))
;; 		full)))))

(defun walker-set-delete (the-walker-set)
  (mapcar (lambda (x) (walker-modify x :modify :delete)) the-walker-set))

(defun walker-set-get-median-params (the-walker-set key &optional take)
  (mapcar (lambda (x) (getf (walker-get x :get :median-params :take take) key)) the-walker-set))

(defun walker-set-get-stddev-params (the-walker-set key &optional take)
  (mapcar (lambda (x) (getf (walker-get x :get :stddev-params :take take) key)) the-walker-set))

(defun walker-set-plot-param (the-walker-set key &optional take)
  (plt:plot (walker-set-get-median-params the-walker-set key take) "w l title \"Param\""))

;; median vs most likely
(defmacro walker-get-f (walker exp &key (take 1000))
  (let ((median-params (walker-get walker :get :median-params :take take)))
    (labels ((replace-param-names (item)
	       (cond ((null item) nil)
		     ((atom item) (if (char= #\: (char (prin1-to-string item) 0))
				      (getf median-params item)
				      item))
		     ((consp item)
		      (cons (replace-param-names (car item))
			    (replace-param-names (cdr item)))))))
      (let ((filled-exp (replace-param-names exp)))
	filled-exp))))

(defun walker-with-exp (walker exp &key (take 1000))
  (let ((median-params (walker-get walker :get :median-params :take take)))
    (labels ((replace-param-names (item)
	       (cond ((null item) nil)
		     ((atom item) (if (char= #\: (char (prin1-to-string item) 0))
				      (getf median-params item)
				      item))
		     ((consp item)
		      (cons (replace-param-names (car item))
			    (replace-param-names (cdr item)))))))
      (let ((filled-exp (replace-param-names exp)))
	(eval filled-exp)))))
(export 'walker-with-exp)


(defun walker-make-step (walker params)
  (make-walker-step :prob (+ (reduce #'+ (map 'list (lambda (ll f d e) (funcall ll f params d e)) (walker-log-liklihood walker) (walker-function walker) (walker-data walker) (walker-data-error walker)))
			     (reduce #'+ (map 'list (lambda (lp d) (funcall lp params d)) (walker-log-prior walker) (walker-data walker))))
		    :params params))

(defun walker-take-step (walker &key l-matrix (temperature 1))
  "The one that uses a full covariance matrix"
  (if (null l-matrix) (setq l-matrix (diagonal-covariance (plist-values (scale-plist 1e-2 (walker-get walker :get :median-params :take 1000))))))
  (let* ((previous-step (walker-last-step walker))
	 (prob0 (walker-step-prob previous-step))
	 (previous-params (walker-step-params previous-step))
	 (cov-values (get-covariant-sample (plist-values previous-params) l-matrix))
	 (next-params (make-plist (plist-keys previous-params) cov-values))
	 (next-step (walker-make-step walker next-params))
	 (prob1 (walker-step-prob next-step))
	 (accepted-step (if (or (> prob1 prob0)
				(> (/ (- prob1 prob0) temperature) (log (random 1.0d0))))
			    next-step
			    previous-step)))
    (walker-modify walker :modify :add-step :new-step accepted-step)))

;; TODO prior bounds implementation - it's a good idea
(defun walker-create (&key function data params data-error log-liklihood log-prior param-bounds)
  "Make a walker with the prescribed characteristics. Can run over multiple sets of functions and data by simply making lists of kwargs. i.e. :function #'my-function vs :function (list #'my-function #'other-function). data, params, etc. must also then be lists that correspond to the order found in the :function kwarg.
:function expects functions formatted with kwargs as parameters:
(lambda (x &key m b &allow-other-keys) (+ b (* -3 m) (* (- m (/ b 60)) x)))
or like this for multiple or linked independent variables
(lambda (x &key m b &allow-other-keys) (+ b (* -3 m (elt x 0)) (* (- m (/ b 60)) (elt x 1))))
:data expects formatting as (list x-data y-data)
:data-error can be a single number (for uniform error) or a list of the same size as the y-data
:params expects a plist like '(:b -1 :m 2)
"
  (declare (ignorable param-bounds))
  (let* ((function (force-list function))
	 (data (to-double-floats (clean-data data (length function))))
	 (data-error (to-double-floats (clean-data-error (if data-error data-error 1) data)))
	 (params (to-double-floats params))
	 (log-liklihood (log-liklihood-fixer (force-list (if log-liklihood log-liklihood (make-list (length (force-list function)) :initial-element #'log-liklihood-normal))) function params data data-error))
	 (log-prior (log-prior-fixer (force-list (if log-prior log-prior (make-list (length (force-list function)) :initial-element #'log-prior-flat))) params data))
	 (walker (make-walker :function function
			      :param-keys (plist-keys params)
			      :last-step (make-walker-step :prob (+ (reduce #'+ (map 'list (lambda (ll f d e) (funcall ll f params d e)) log-liklihood function data data-error))
								    (reduce #'+ (map 'list (lambda (lp d) (funcall lp params d)) log-prior data)))
							   :params params)
			      :walk (list (make-walker-step :prob (+ (reduce #'+ (map 'list (lambda (ll f d e) (funcall ll f params d e)) log-liklihood function data data-error))
								     (reduce #'+ (map 'list (lambda (lp d) (funcall lp params d)) log-prior data)))
							    :params params))
			      :data data
			      :data-error data-error
			      :log-liklihood log-liklihood
			      :log-prior log-prior)))
    walker))

(defun mcmc-fit (&key function data params data-error log-liklihood log-prior param-bounds)
  (declare (ignorable param-bounds))
  (let* ((walker (walker-create :function function
				:data data
				:params params
				:data-error data-error
				:log-liklihood log-liklihood
				:log-prior log-prior
				:param-bounds param-bounds)))
    (walker-adaptive-steps walker)
    walker))
(export 'mcmc-fit)

(defvar example-function (lambda (x &key m b &allow-other-keys) (+ b (* -3 m) (* (- m (/ b 60)) x))))
(defvar example-data '((-4 -1 2 5 10) (0 2 5 9 13)))
(defvar example-params '(:b -1 :m 2))
(defvar example-error 0.2)
(defvar example-log-liklihood #'log-liklihood-normal)
(defvar example-log-prior #'log-prior-flat)

;; (defparameter woi (mcmc-fit :function (lambda (x &key m b &allow-other-keys) (+ b (* -3 m) (* (- m (/ b 60)) x))) :data '((-4 -1 2 5 10) (0 2 5 9 13)) :params '(:b -1 :m 2) :data-error 0.2))

;; (defparameter woi (mcmc-fit :function (list (lambda (x &key b m c d &allow-other-keys) (+ b (* m x) (* c x x) (* d x x x))) (lambda (x &key e m g &allow-other-keys) (+ e (* (+ m g) x)))) :data '(((0 1 2 3 4) (4 5 6 8 4)) ((10 20 30 40 50) (1 2 3 4 5))) :params '(:b -1 :m 2 :c 1 :d -1 :e 0.5 :g -2) :data-error 0.2))

(defun walker-diagnose-params (walker &key (style (or :step-insert :plot-flow)) params)
  (case style
    (:step-insert (push (walker-make-step walker params) (walker-walk walker)))
    (:plot-flow (print "Implement some (10) steps of the walk in plot space to see evolution"))))
(export 'walker-diagnose-params)


;;; Walker and fit visualization
(defun walker-get-data-and-fit-no-stddev (walker &key (take 1000) (x-column 0) (y-column 1) (fn-number 0) (which-solution (or :most-likely :median)) x-shift y-shift)
  "Extract plottable data as well as fit and stddev data from a walker."
  (let* ((take (if (or (null take) (> take (walker-length walker)))
		   (walker-length walker)
		   take))
	 (fn (elt (walker-function walker) fn-number))
	 (data (elt (walker-data walker) fn-number))
	 (x-data (elt data x-column))
	 (x-data-shift (if (null x-shift) x-data (map 'vector (lambda (x) (+ x-shift x)) x-data)))
	 (y-data (elt data y-column))
	 (y-data-shift (if (null y-shift) y-data (map 'vector (lambda (y) (+ y-shift y)) y-data)))
	 (x-fit (linspace (reduce #'min x-data) (reduce #'max x-data) :len 1000))
	 (x-fit-shift (if (null x-shift) x-fit (map 'vector (lambda (x) (+ x-shift x)) x-fit)))
	 (params (case which-solution
		   (:most-likely (walker-step-params (walker-get walker :get :most-likely-step :take take)))
		   (:median (walker-get walker :get :median-params :take take))
		   (otherwise (walker-get walker :get :median-params :take take))))
	 (y-fit (mapcar (lambda (x) (apply fn x params)) x-fit))
	 (y-fit-shift (if (null y-shift) y-fit (map 'vector (lambda (y) (+ y-shift y)) y-fit))))
    (list x-fit-shift y-fit-shift x-data-shift y-data-shift params)))
(export 'walker-get-data-and-fit-no-stddev)

(defun walker-get-data-and-fit (walker &key (take 1000) (x-column 0) (y-column 1) (fn-number 0) (which-solution (or :most-likely :median)) x-shift y-shift)
  "Extract plottable data as well as fit and stddev data from a walker."
  (let* ((take (if (or (null take) (> take (walker-length walker)))
		   (walker-length walker)
		   take))
	 (fn (elt (walker-function walker) fn-number))
	 (data (elt (walker-data walker) fn-number))
	 (x-data (elt data x-column))
	 (x-data-shift (if (null x-shift) x-data (map 'vector (lambda (x) (+ x-shift x)) x-data)))
	 (y-data (elt data y-column))
	 (y-data-shift (if (null y-shift) y-data (map 'vector (lambda (y) (+ y-shift y)) y-data)))
	 (x-fit (linspace (reduce #'min x-data) (reduce #'max x-data) :len 1000))
	 (x-fit-shift (if (null x-shift) x-fit (map 'vector (lambda (x) (+ x-shift x)) x-fit)))
	 (params (case which-solution
		   (:most-likely (walker-step-params (walker-get walker :get :most-likely-step :take take)))
		   (:median (walker-get walker :get :median-params :take take))
		   (otherwise (walker-get walker :get :median-params :take take))))
	 (y-fit (mapcar (lambda (x) (apply fn x params)) x-fit))
	 (y-fit-shift (if (null y-shift) y-fit (map 'vector (lambda (y) (+ y-shift y)) y-fit)))
	 (previous-walks (walker-get walker :get :steps))
	 (sorted-walks (subseq (sort previous-walks #'> :key #'walker-step-prob) 0 (ceiling (* 0.66 take))))
	 (all-ys (mapcar (lambda (x) (mapcar (lambda (walk) (apply fn x (walker-step-params walk))) sorted-walks)) x-fit))
	 (max-ys (mapcar (lambda (y) (+ (if y-shift y-shift 0) (reduce #'max y))) all-ys))
	 (min-ys (mapcar (lambda (y) (+ (if y-shift y-shift 0) (reduce #'min y))) all-ys)))
    (list x-fit-shift max-ys min-ys y-fit-shift x-data-shift y-data-shift params)))
(export 'walker-get-data-and-fit)


(defun walker-plot-data-and-fit (walker &key (take 1000) (x-column 0) (y-column 1) (fn-number 0) (image-size '(1920 1080)) (which-solution (or :most-likely :median)) x-shift y-shift)
  (destructuring-bind (x-fit max-ys min-ys y-fit x-data y-data params) (walker-get-data-and-fit walker :take take :x-column x-column :y-column y-column :fn-number fn-number :which-solution which-solution :x-shift x-shift :y-shift y-shift)
    (plt:reset)
    (plt:send-command :terminal (apply #'format nil "qt size ~d,~d linewidth 3 font \"Arial,18\"" image-size)
		      :format "x \"%.4g\""
		      :format "y \"%.4g\""
		      :xlabel "\"x-data\""
		      :ylabel "\"y-data\"")
    (plt:plot x-fit max-ys "w l lc rgb \"green\" title \"fit stddev upper limit\""
	      x-fit min-ys "w l lc rgb \"green\" title \"fit stddev lower limit\""
	      x-fit y-fit "w l lc rgb \"red\" title \"fit\""
	      x-data y-data "w p pt 6 ps 2 lc rgb \"black\" title \"data\"")
    params))

(defun walker-plot-residuals (walker &key (take 1000) (x-column 0) (y-column 1) (fn-number 0))
  (let* ((take (if (or (null take) (> take (walker-length walker)))
		   (walker-length walker)
		   take))
	 (fn (elt (walker-function walker) fn-number))
	 (data (elt (walker-data walker) fn-number))
	 (stddev (elt (walker-data-error walker) fn-number))
	 (x-data (elt data x-column))
	 (y-data (elt data y-column))
	 (stddev (if (= 1 (length stddev)) (make-array (length x-data) :initial-element (elt stddev 0)) stddev))
	 (params (walker-get walker :get :median-params :take take))
	 (y-fit (map 'vector (lambda (x) (apply fn x params)) x-data))
	 (y-residuals (map 'vector (lambda (yf y) (- yf y)) y-fit y-data)))
    (plt:reset)
    (plt:send-command :terminal "qt size 1920,1080 linewidth 3 font \"Arial,18\""
		      :format "x \"%.4g\""
		      :format "y \"%.4g\""
		      :xlabel "\"x-data\""
		      :ylabel "\"y-data\"")
    (plt:plot x-data y-residuals "w p pt 6 ps 2 lc rgb \"black\" title \"residuals\""
	      x-data stddev "w p pt 2 ps 1 lc rgb \"red\" title \"point error\""
	      x-data (make-array (length x-data) :initial-element 0) "w l lc rgb \"red\" title \"baseline\"")))

(defun walker-catepillar-plots (walker &optional take (x-scale 1) (y-scale 1))
  (let* ((param-keys (walker-param-keys walker))
	 (n (length param-keys)))
    (plt:reset)
    (plt:send-command :terminal (format nil "pngcairo size ~d,~d linewidth 1 font \"Arial,18\"" (* x-scale 1920) (* y-scale 1080))
		      :output "\"temp.png\""
		      :format "y \"%.1e\""
		      :format "x \"%.1e\"")
    (apply #'plt:multiplot
	   (format nil "layout ~d,1" n)
	   (apply #'append
		  (mapcar-enum (key i param-keys)
		    (list (list (walker-get walker :get :param :param key :take take) (format nil "w l title \"~a\"" key))
			  (if (= i (1- n))
			      (list :xlabel "\"Step\"" :ylabel (format nil "\"~a\"" key))
			      (list :xlabel "" :ylabel (format nil "\"~a\"" key)))))))))

(defun walker-liklihood-plot (walker &optional take)
  (let ((probs (walker-get walker :get :log-liklihoods :take take)))
    (plt:reset)
    (plt:send-command :terminal "qt size 1920,1080 linewidth 1 font \"Arial,18\""
		      :xlabel "\"Step\""
		      :ylabel "\"Log Liklihood\"")
    (plt:plot probs "w l title \"Liklihood\"")))

(defun walker-2d-plot (walker keys &optional take)
  (let* ((x-key (first keys))
	 (y-key (second keys))
	 (x (walker-get walker :get :param :param x-key :take take))
	 (y (walker-get walker :get :param :param y-key :take take)))
    (plt:reset)
    (plt:send-command :terminal "qt size 1920,1080 linewidth 1 font \"Arial,18\""
		      :xlabel (format nil "\"~a\"" x-key)
		      :ylabel (format nil "\"~a\"" y-key))
    (plt:plot x y "w p pt 3 ps 4 title \"My plot\"")))

(defun walker-all-2d-plots (walker &key take (size '(3240 3240)))
  (labels ((permute-params (params)
	     (let ((return-list nil))
	       (do ((i 0 (1+ i)))
		   ((>= i (1- (length params))) (reverse return-list))
		 (do ((j (1+ i) (1+ j)))
		     ((>= j (length params)) nil)
		   (push (list (elt params i) (elt params j)) return-list)))))
	   (make-origin-list (params)
	     (let ((return-list nil)
		   (len (length params)))
	       (do ((i 0 (1+ i)))
		   ((>= i (1- (length params))) (reverse return-list))
		 (do ((j (1+ i) (1+ j)))
		     ((>= j (length params)) nil)
		   (push (list (* i (/ 1.0 (1- len))) (* (1- j) (/ 1.0 (1- len)))) return-list))))))
    (let* ((param-keys (walker-param-keys walker))
	   (all-params (apply #'append (partition (make-plist param-keys (mapcar (lambda (x) (walker-get walker :get :param :param x :take take)) param-keys)) 2)))
	   (key-pairs (permute-params param-keys))
	   (plot-data (mapcar (lambda (x) (destructuring-bind (key1 key2) x (list (getf all-params key1) (getf all-params key2)))) key-pairs))
	   (origin-list (make-origin-list param-keys)))
      (plt:send-plot-options :terminal (format nil "pngcairo size ~d,~d linewidth 1 font \"Arial,12\"" (elt size 0) (elt size 1)) :output "\"temp.png\"")
      (apply #'plt:multiplot ""
	     (apply #'nconc (mapcar-enum (element index plot-data)
			      (list (list (elt element 0) (elt element 1) "w p"
					  :origin (apply #'format nil "~,1f,~,1f" (elt origin-list index)) :size (format nil "~1,f,~1,f" (/ 1.0 (1- (length param-keys))) (/ 1.0 (1- (length param-keys)))) :xlabel (format nil "\"~a\"" (elt (elt key-pairs index) 0)) :ylabel (format nil "\"~a\"" (elt (elt key-pairs index) 1)))))))
      (plt:send-plot-options :output "unset"))))

(defun walker-param-histo (walker key &key (take 10000) (bins 20))
  (let* ((param-values (sort (walker-get walker :get :param :param key :take take) #'<))
	 (histo-x (make-histo-x param-values bins))
	 (histo (make-histo param-values bins)))
    (plt:reset)
    (plt:send-command :terminal "qt size 1920,1080 linewidth 3 font \"Arial,18\""
		      :xlabel (format nil "\"~a\"" key)
		      :ylabel "\"Counts\"")
    (plt:plot histo-x histo (format nil "w lp ps 3 title \"Histogram: ~a\"" key))))

(defun show ()
  "Shows the plots generated from 2d histogram and catepillar functions"
  (uiop:run-program "feh ./temp.png -3 5"))


;;; File Management

(defstruct parsed-file
  (filename "" :type string)
  (header-lines nil :type list)
  (data-lines nil :type list))

;; File searching
(defun walk-dirs (dir)
  (if (uiop:directory-exists-p dir) ;; it's a dir?
      (let ((items (nconc (uiop:subdirectories dir) (uiop:directory-files dir))))
	(mapcar (lambda (d) (walk-dirs d)) items))
      dir))

(defun get-filename (dir &key include exclude)
  "Returns all files under dir (and its subdirs) whose name (including directory name) matches ALL specified patterns"
  (let* ((include (when include (if (not (listp include)) (list include) include)))
	 (exclude (when exclude (if (not (listp exclude)) (list exclude) exclude)))
	 (all-files (flatten (walk-dirs dir)))
	 (files (remove-if-not (lambda (y)
				 (and (every (lambda (g) (search g (namestring y))) include)
				      (notany (lambda (g) (search g (namestring y))) exclude)))
			       all-files)))
    (if (= (length files) 1)
	(car files)
	files)))

;; File gathering, reading, and quick plotting
(defun read-file-lines (filename)
  "Read 'filename' line by line and return a list of lines"
  (with-open-file (s filename)
    (let ((data nil))
      (do ((line (read-line s nil nil) (read-line s nil nil)))
	  ((not line) (reverse data))
	(push line data)))))
(export 'read-file-lines)

(defun separate-header-and-data (file-data number-of-header-lines)
  (list (subseq file-data 0 number-of-header-lines)
	(subseq file-data number-of-header-lines)))
(export 'separate-header-and-data)

(defun auto-split-and-read-csv (unsplit-data-lines)
  (let* ((delim-counts (mapcar (lambda (x) (list x (count x (elt unsplit-data-lines 0)))) '(#\tab #\, #\; #\:)))
	 (most-likely-delim (elt (reduce (lambda (x y) (if (> (elt x 1) (elt y 1)) x y)) delim-counts) 0)))
    (remove-if
     (lambda (x) (every #'null x))
     (mapcar-tree
      (lambda (x) (read-from-string x nil nil))
      (plt::transpose
       (mapcar
	(lambda (x) (split-string most-likely-delim x))
	unsplit-data-lines))))))
(export 'auto-split-and-read-csv)

(defun file->file-specs (filename &key (delim #\tab))
  "Tells the user various metrics about the file. Lines, header lines, data lines, data rows, data sets."
  (with-open-file (in filename :direction :input)
    (labels ((get-lines (&optional (num-lines 0) (found-data nil) (data-length nil) (data-rows nil))
	       (let ((line (string-right-trim '(#\return) (read-line in nil nil)))) ; allow for Windows files
		 (cond ((string= line "NIL") ; end of file - return info ; ordered by freque
			(list :file-lines num-lines :header-lines found-data :data-length data-length :data-rows (if data-rows data-rows (- num-lines found-data)) :num-pages (if data-rows (floor (- num-lines found-data) data-rows) 1)))
		       ((and (string= line "") found-data (not data-rows)) ; set data rows
			(get-lines num-lines found-data data-length (- num-lines found-data)))
		       ((string= line "") ; ignore line
			(get-lines num-lines found-data data-length data-rows))
		       ((and (numberp (read-from-string (elt (split-string #\tab line) 0))) (not found-data)) ; first line w nums?
			(get-lines (+ num-lines 1) num-lines (length (split-string delim line)) data-rows))
		       (t ; regular line - just increase lines
			(get-lines (+ num-lines 1) found-data data-length data-rows))))))
      (get-lines))))

(defun make-data-into-pages (data file-specs)
  (let ((ret-tree (make-list (getf file-specs :num-pages))))
    (mapcar (lambda (x a b)
	      (setf x (mapcar (lambda (y)
				(subseq y a b))
			      data)))
	    ret-tree
	    (butlast (linspace 0 (* (getf file-specs :num-pages) (getf file-specs :data-rows)) :len (1+ (getf file-specs :num-pages)) :type 'integer))
	    (rest (linspace 0 (* (getf file-specs :num-pages) (getf file-specs :data-rows)) :len (1+ (getf file-specs :num-pages)) :type 'integer)))))

(defun read-file->data (filename &key (file-specs nil) (delim #\tab) (transpose t) pages)
  "Reads a file into a data list where 'elt 0' of the list is the first column by default.\n
If row based, set columns to nil"
  (unless file-specs
    (setq file-specs (file->file-specs filename :delim delim)))
  (let* ((header-lines (getf file-specs :header-lines))
	 (file-contents nil)
	 (line nil))
    (labels ((remove-header-lines (n stream)
	       (repeat n (read-line stream nil nil)))
	     (read-file (stream)
	       (when (setq line (read-line stream nil nil))
		 (setq line (string-right-trim '(#\return) line))
		 (let ((vals (mapcar #'read-from-string (split-string delim line))))
		   (when vals
		     (push vals file-contents)))
		 (read-file stream))))
      (with-open-file (in filename :direction :input)
	(remove-header-lines header-lines in)
	(read-file in)
	(let ((file-contents
		(if transpose
		    (transpose-matrix (reverse file-contents))
		    (reverse file-contents))))
	  (if pages
	      (make-data-into-pages file-contents file-specs)
	      file-contents))))))


(defun read-file->plot (filename &optional (x-column 0) (y-column 1))
  (let ((data (read-file->data filename)))
    (plt:plot (elt data x-column) (elt data y-column))))
(export 'read-file->plot)

(defun read-files->plot (filenames &optional (x-column 0) (y-column 1))
  (let ((data (mapcar #'read-file->data filenames)))
    (apply #'plt:plot (apply #'append (mapcar (lambda (e) (list (elt e x-column) (elt e y-column))) data)))))
(export 'read-files->plot)


;;; Statistics Functions
(defun multivariate-gaussian-random (covs)
  (mapcar (lambda (x) (* x (alexandria:gaussian-random))) covs))

(defun nth-percentile (n sequence &optional (sorted nil))
  (let* ((copy sequence)
	 (len (length sequence))
	 (n (* n (- len 1) 1/100)))
    (multiple-value-bind (pos rem) (floor n)
      (unless sorted
	(setq copy (sort (copy-seq sequence) #'<)))
      (if (= rem 0)
	  (elt copy pos)
	  (let ((e1 (elt copy pos))
		(e2 (elt copy (+ pos 1))))
	    (/ (+ e1 e2) 2))))))

(defun 95cr (sequence)
  (list (nth-percentile 2.5 sequence) (nth-percentile 97.5 sequence)))

(defun iqr (sequence &optional (sorted nil))
  (unless sorted (setq sorted (sort (copy-seq sequence) #'<)))
  (- (nth-percentile 75 sequence t) (nth-percentile 25 sequence t)))

(defun median (sequence &optional sorted)
  (nth-percentile 50 sequence sorted))

(defun mean (sequence)
  (/ (reduce #'+ sequence) (length sequence)))

(defun variance (sequence)
  (let* ((mean (/ (reduce #'+ sequence) (length sequence)))
	 (sum-sq (reduce #'+ (map 'list (lambda (x) (expt (- x mean) 2)) sequence))))
    (/ sum-sq (- (length sequence) 1))))

(defun standard-deviation (sequence)
  (sqrt (variance sequence)))

(defun standard-deviation-normal (sequence &optional sorted)
  (let* ((copy sequence))
    (unless sorted
      (setq copy (sort (copy-seq sequence) #'<)))
    (let ((middle-value (median copy t))
	  (upper-variance (nth-percentile 84.1 copy t)))
      (- upper-variance middle-value))))

(defun variance-normal (sequence &optional sorted)
  (expt (standard-deviation-normal sequence sorted) 2))


;;; binning / histogram
(defun make-histo (sequence &optional num-bins)
  (let* ((bottom (reduce #'min sequence))
	 (top (reduce #'max sequence))
	 (num-bins (if num-bins num-bins (floor (* (- top bottom) (expt (length sequence) 1/3)) (* 2 (iqr sequence)))))
	 (num-bins (+ 1 num-bins))
	 (boundaries (linspace bottom top :len num-bins)))
    (labels ((count-for-bins (n sequence &optional return-list)
	       (if (< n num-bins)
		   (let ((pos (position (elt boundaries n) sequence :test-not #'>=)))
		     (unless pos
		       (setq pos (length sequence)))
		     (count-for-bins (+ n 1)
				     (subseq sequence pos)
				     (cons pos return-list)))
		   (reverse return-list))))
      (count-for-bins 1 sequence))))

(defun make-histo-x (sequence &optional num-bins)
  (let* ((bottom (reduce #'min sequence))
	 (top (reduce #'max sequence))
	 (num-bins (if num-bins num-bins (floor (* (- top bottom) (expt (length sequence) 1/3)) (* 2 (iqr sequence)))))
	 (gap-size (/ (- top bottom) num-bins)))
    (linspace (+ bottom (/ gap-size 2)) top :len num-bins)))

