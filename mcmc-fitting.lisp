;;; mcmc-fitting.lisp

;;; Example basic use
;; Basic line fit
;; (defparameter woi (mfit:mcmc-fit :function (lambda (x &key m b &allow-other-keys) (+ b (* m x))) :data '((-4 -1 2 5 10) (0 2 5 9 13)) :params '(:b -1 :m 2) :data-error 0.2))

;; Single param list fit
;; (defparameter woi (mfit:mcmc-fit :function (lambda (x &key params &allow-other-keys) (+ (elt params 0) (* (elt params 1) x))) :data '((-4 -1 2 5 10) (0 2 5 9 13)) :params '(:params (-1 2)) :data-error 0.2))

;; Single param vector fit
;; (defparameter woi (mfit:mcmc-fit :function (lambda (x &key params &allow-other-keys) (+ (elt params 0) (* (elt params 1) x))) :data '((-4 -1 2 5 10) (0 2 5 9 13)) :params '(:params #(-1 2)) :data-error 0.2))

;; Single param array fit
;; (defparameter woi (mfit:mcmc-fit :function (lambda (x &key params &allow-other-keys) (+ (aref params 0 0) (* (aref params 1 0) x))) :data '((-4 -1 2 5 10) (0 2 5 9 13)) :params (list :params (make-array '(2 1) :element-type 'double-float :initial-contents '((-1d0) (2d0)))) :data-error 0.2))

;; Global parameter fit (polynomial and line)
;; (defparameter woi (mfit:mcmc-fit :function (list (lambda (x &key b m c d &allow-other-keys) (+ b (* m x) (* c x x) (* d x x x))) (lambda (x &key e m g &allow-other-keys) (+ e (* (+ m g) x)))) :data '(((0 1 2 3 4) (4 5 6 8 4)) ((10 20 30 40 50) (1 2 3 4 5))) :params '(:b -1 :m 2 :c 1 :d -1 :e 0.5 :g -2) :data-error 0.2))


;; test speed of sample space rotation via l-matrix (matrix multiplication)
;; Arrays are slower without optimization. 10x faster with optimization (basic array multiplication)
;; Performing calculation with arrays while transforming from lists to arrays is 1.37x faster
;; Vectors are slower without optimization. 2x faster with optimization
#|
(defun rotate-list (l-lists sample)
  (mapcar (lambda (x) (list (reduce #'+ (mapcar (lambda (s l) (* (elt s 0) l)) sample x)))) l-lists))

(defun rotate-vector (l-vectors sample)
  (declare (optimize speed debug (safety 0))
	   (simple-vector sample l-vectors))
  (map 'simple-vector (lambda (x) (declare (simple-vector x)) (list (reduce #'+ (map 'simple-vector (lambda (s l) (* (the double-float (elt s 0)) (the double-float l))) sample x)))) l-vectors))

(defun rotate-array (l-array sample)
  (declare (optimize speed debug (safety 0))
	   ((simple-array double-float *) l-array sample))
  (let* ((imax (array-dimension l-array 0))
	 (jmax (array-dimension l-array 1))
	 (kmax (array-dimension sample 1))
	 (return-array (make-array (list imax kmax))))
    (do ((i 0 (1+ i)))
	((>= i imax) return-array)
      (do ((k 0 (1+ k)))
	  ((>= k kmax))
	(setf (aref return-array i k)
	      (do ((j 0 (1+ j))
		   (mini-sum 0d0))
		  ((>= j jmax) mini-sum)
		(declare (double-float mini-sum))
		(incf mini-sum (* (aref l-array i j) (aref sample j k)))))))))

(progn
  (defparameter rank 3000)
  (defparameter l-lists (mfit::repeat rank (mfit::repeat rank (random 1d0))))
  (defparameter l-vectors (map 'simple-vector (lambda (x) (coerce x 'simple-vector)) l-lists))
  (defparameter l-array (make-array (list rank rank) :element-type 'double-float :initial-contents l-lists))
  (defparameter sample-list (mfit::repeat rank (list (random 1d0))))
  (defparameter sample-vector (coerce sample-list 'simple-vector))
  (defparameter sample-array (make-array (list rank 1) :element-type 'double-float :initial-contents sample-list))
  (defparameter rotated-sample (time (null(rotate-list l-lists sample-list))))
  (defparameter rotated-vectors (time (null (rotate-vector l-vectors sample-vector))))
  (defparameter rotated-array (time (null (rotate-array l-array sample-array))))
  (time (progn (rotate-array (make-array (list rank rank) :element-type 'double-float :initial-contents l-lists)
			     (make-array (list rank 1) :element-type 'double-float :initial-contents sample-list)))))
|#

;; testing function evaluation speed over a list of data vs arrays of data
;; vector is not faster than list
;; 2x speed up by optimizing for single/double floats
#|
(defun myfunc (x &key p1 p2 p3 &allow-other-keys)
  "Simple Gaussian"
  (declare (optimize (speed 3) (safety 1) (debug 1))
	   (double-float x p1 p2 p3))
  (* p3 (exp (- (expt (/ (- x p1) p2) 2)))))

(progn
  (defparameter params (list :p1 20d0 :p2 5d0 :p3 14d0))
  (defparameter x-list (mfit:linspace 0 50 :len 100 :type 'double-float))
  (defparameter x-vector (coerce x-list 'simple-vector))
  (time (null (mfit::repeat 100000 (mapcar (lambda (x) (apply #'myfunc x params)) x-list))))
  (time (null (mfit::repeat 100000 (let* ((func-x (lambda (x) (apply #'myfunc x params))) (vec-len (length x-vector)) (return-vec (make-array (list vec-len) :element-type 'double-float))) (do ((i 0 (1+ i))) ((>= i vec-len) return-vec) (setf (aref return-vec i) (funcall func-x (svref x-vector i)))))))))
|#

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
(defparameter posdefmat (matrix-add (scale-array 10 (eyes 3)) (scale-array 0.5 (matrix-add randmat randmat))))
|#

(in-package :mcmc-fitting)


;;; Utility
(defmacro br (&body body)
  `(break "break ~s" ,@body))

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

(defun thin (list n &optional (start 0))
  "Returns every 'n'th element from 'list' starting with 'start'."
  (if (= n 1)
      list
      (do ((i 0 (1+ i))
	   (return-list nil))
	  ((= i (length list)) (reverse return-list))
	(when (= 0 (mod (- i start) n))
	  (push (elt list i) return-list)))))

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

(defun map-tree (function tree)
  "Maps 'function' over the deepest elements of 'tree'. 'tree' can be composed of any combination of sequences.
E.g. (map-tree #'1+ '((1 2) (2 3) (3 4) (100 (12 23 (324 54)))))
 =>  '((2 3) (3 4) (4 5) (101 (13 24 (325 55))))."
  (cond ((null tree)
	 nil)
	((and (atom tree) (or (arrayp tree) (stringp tree) (not (vectorp tree))))
	 (funcall function tree))
	(t
	 (map (type-of tree) (lambda (x) (map-tree function x)) tree))))
(export 'map-tree)

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

(defun array-to-plist (keys values-array)
  (make-plist keys (mapcar (lambda (x) (aref values-array x 0)) (range (array-dimension values-array 0)))))

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
		((>= n len) (reverse return-list)))))
      (values (case type
		((or single-float float double-float rational) (mapcar (lambda (x) (coerce x type)) rational-numbers))
		(t (mapcar #'round rational-numbers)))
	      step))))

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

(defun diff-matrix (matrix)
  "Find the difference between consecutive items in a list of lists, vectors, or arrays. Returns a list of lists."
  (mapcon (lambda (x)
	    (if (nthcdr 1 x)
		(list (map 'list (lambda (x y)
				   (etypecase x
				     (number (- x y))
				     (list (mapcar #'- x y))
				     (vector (map 'list #'- x y))
				     (array (let ((return-list nil)) (dotimes (i (array-dimension x 0) return-list) (push (- (aref x i 0) (aref y i 0)) return-list))))))
			   (cadr x) (car x)))
		nil))
	  matrix))

(defun diff-lplist (lplist)
  (let* ((keys (get-plist-keys (elt lplist 0)))
	 (values (mapcar #'plist-values lplist)))
    (mapcar (lambda (x) (make-plist keys x)) (diff-matrix values))))

(defun partition (list n)
  "Pairs elements and removes any with full matches"
  (labels ((rec (l &optional (acc nil))
	     (if (nthcdr (- n 1) l)
		 (rec (subseq l n) (cons (subseq l 0 n) acc))
		 (cons l acc))))
    (reverse (cdr (rec list)))))

(defun transpose (xy-list)
  "Takes 'xy-list' of form ((x1 y1 z1 ...) (x2 y2 z2 ...) ...) and turns it into ((x1 x2 ...) (y1 y2 ...) (z1 z2 ...) ...). Also known as a transpose. It is its own inverse. Works for combinations of lists and vectors."
  (when xy-list
    (apply #'map 'list #'list xy-list)))

(defun list-of-arrays-transpose (l-o-array)
  "Transpose a list of column vectors (2d arrays)."
  (declare (optimize speed (safety 0))
	   (list l-o-array))
  (let ((return-list nil))
    (dotimes (j (array-dimension (the (array double-float (* *)) (car l-o-array)) 0) return-list)
      (push (mapcar
	     (lambda (x)
	       (declare ((simple-array double-float) x))
	       (aref x j 0))
	     l-o-array)
	    return-list))))

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
  "Always 0d0."
  (declare (ignore params data))
  0d0)

;; TODO let this work with positions as well
(defmacro prior-bounds-let ((&rest keys-low-high) &body body)
  "Helps construct useful bounds on parameters. -1d10 penalty for going outside of specified low and high values with an exponential gradient to help the algorithm find it's way back/prevent it from accidentally leaving. Anaphor 'bounds-total' should be used as part of the return value of the body. Assumes multiple kwarg usage of the walker (Needs variables names, not positions). Used as...
==> (prior-bounds-let ((x -10 10) (y 100 200) (z 1d-4 1d-2)) (+ bounds-total (* 4 x-bound)))"
  #+sbcl (declare (sb-ext:muffle-conditions sb-ext:compiler-note))
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
(defun log-normal (x mu sigma)
  #+sbcl (declare (sb-ext:muffle-conditions sb-ext:compiler-note)
	   (optimize speed)
	   (double-float mu x)
	   (type (double-float 0d0 *) sigma))
  (+ (* -1/2 (log (* 2 pi))) (* -1 (log sigma)) (* -1/2 (expt (/ (- x mu) sigma) 2d0))))

(defun log-factorial (n)
  (reduce (lambda (x y) (+ x (log y))) (up-to n)))

(defun log-poisson (lambda k)
  (- (* k (log lambda)) lambda (log-factorial k)))

;; log liklihood is the function that determines how data is useful
;; It does the pattern matching and analysis on input variables
;; The liklihood defines the result
;; Each independent variable needs to have the appropriate distribution
;; Then they just all add together (for the log version)
;; The prior is then tacked on at the end however the user sees fit
;; Since the prior does not change, it will just speed up the convergence to have a good one
;; prior could also depend on data for bounds or such
(defun log-liklihood-normal (fn params data stddev)
  (declare (optimize speed)
	   (list data stddev)
	   (function fn))
  (let* ((x (nth 0 data))
	 (y (nth 1 data)))
    (declare (list x y))
    (reduce #'+ (mapcar (lambda (x y z) (log-normal y (apply fn x params) z)) x y stddev))))

(defun create-log-liklihood-function (log-liklihood-function)
  "Input a function that determines the log-liklihood of some data point. 'log-liklihood-function' needs to accept 3 arguments. This function assumes your independent (x) variable is in the 0th column of the data and the dependent (y) variable is in the 1st column. If you have more elaborate needs, examine the implementation of #'create-log-liklihood-normal and modify to your needs.
The first is the actual y value of the point.
The second is the model predicted value of the point.
The third is the error of the point.
i.e. (lambda (y model error) (/ (- y model) error))
i.e. (lambda (y model error) (declare (ignore error)) (/ (- y model) 1)) ;; you can ignore any arguments you want, but the function must accept 3 arguments."
  (lambda (fn params data stddev)
    (declare (optimize speed)
	     (list data stddev)
	     (function log-liklihood-function fn))
    (let* ((x (nth 0 data))
	   (y (nth 1 data)))
      (declare (list x y))
      (reduce #'+ (mapcar (lambda (x y) (funcall log-liklihood-function y (apply fn x params) stddev)) x y)))))
(export 'create-log-liklihood-function)

(defun log-liklihood-normal-cutoff (fn params data stddev)
  (declare (optimize speed)
	   (list data stddev)
	   (function fn))
  (let* ((x (elt data 0))
	 (y (elt data 1)))
    (declare (list x y))
    (reduce #'+ (mapcar (lambda (x y z) (max -5000d0 (the double-float (log-normal y (apply fn x params) z)))) x y stddev))))
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
	 (ret-values (map-tree fn values)))
    (make-plist keys ret-values)))

(defun scale-plist (scale plist)
  (map-plist (lambda (x) (* x scale)) plist))



;;; Walker and associated functions
(defstruct walker-step
  (prob most-negative-double-float :type float)
  (params nil :type list))
(export '(walker-step make-walker-step walker-step-prob walker-step-params))

(defstruct walker
  (function nil :type list)
  (param-keys nil :type list)
  (param-style (or :multiple-kwargs :single-item) :type symbol)
  (walk nil :type list)
  (length 1 :type integer)
  (age 1 :type integer)
  (most-likely-step (make-walker-step) :type walker-step)
  (last-step (make-walker-step) :type walker-step)
  (data nil :type (or list vector))
  (data-error nil :type list)
  (log-liklihood nil :type list)
  (log-prior nil :type list))
(export '(walker make-walker walker-function walker-param-keys walker-walk walker-length walker-age walker-most-likely-step walker-last-step walker-data walker-data-error walker-log-liklihood walker-log-prior))

;; Useful for inspecting when algorithm fails
(defun walker-check-for-complex-walks (walker &optional take)
  (let ((complex-walks-ps (mapcar (lambda (x) (some #'identity x)) (map-tree #'complexp (lplist-to-l-matrix (diff-lplist (walker-get walker :get :unique-steps :take take)))))))
    (if (some #'identity complex-walks-ps) complex-walks-ps nil)))

(defun walker-get (walker &key (get (or :steps :unique-steps :forward-steps :most-likely-step :acceptance :param :params :most-likely-params :median-params :all-params :stddev-params :log-liklihoods :covariance-matrix :l-matrix)) take param)
  "Get any of the above items from a walker. Use the 'take' kwarg to limit the number of walks you are sampling over. Use the 'param' kwarg to specify which params for the :param :get option."
  (case get
    (:steps (subseq (walker-walk walker) 0 (when take (min (walker-length walker) take))))
    ;; TODO Bad names, should give step list, not just params
    (:unique-steps (mapcon (lambda (x)
			     (if (equal (walker-step-prob (car x)) (when (cadr x) (walker-step-prob (cadr x))))
				 nil
				 (list (walker-step-params (car x)))))
			   (walker-get walker :get :steps :take take)))
    (:forward-steps (mapcon (lambda (x)
			      (when (cadr x)
				(if (<= (walker-step-prob (car x)) (walker-step-prob (cadr x)))
				    nil
				    (list (walker-step-params (car x))))))
			    (walker-get walker :get :steps :take take)))
    (:most-likely-step (reduce (lambda (x y)
				 (if (> (walker-step-prob x) (walker-step-prob y)) x y))
			       (walker-get walker :get :steps :take take)))
    (:acceptance
     (let* ((probs (mapcar #'walker-step-prob (walker-get walker :get :steps :take take))))
       (/ (length (remove-consecutive-duplicates probs)) (length probs))))
    (:param (mapcar (lambda (x) (getf (walker-step-params x) param)) (walker-get walker :get :steps :take take)))
    (:params (mapcar #'walker-step-params (walker-get walker :get :steps :take take)))
    (:most-likely-params
     (let ((param-keys (walker-param-keys walker))
	   (most-likely-step (walker-most-likely-step walker)))
       (make-plist param-keys
		   (mapcar (lambda (key) (getf (walker-step-params most-likely-step) key)) param-keys))))
    (:median-params
     (let ((param-keys (walker-param-keys walker))
	   (steps (walker-get walker :get :steps :take take)))
       (make-plist param-keys
		   (typecase (cadr (walker-step-params (car steps)))
		     (number (mapcar (lambda (key) (median (mapcar (lambda (step) (getf (walker-step-params step) key)) steps))) param-keys))
		     ((or list vector) (list (mapcar #'median (transpose (mapcar (lambda (step) (getf (walker-step-params step) (elt param-keys 0))) steps)))))
		     (array (list (mapcar #'median (list-of-arrays-transpose (mapcar (lambda (step) (getf (walker-step-params step) (elt param-keys 0))) steps)))))))))
    ;; TODO add failsafe here if not enough valid steps to do cholesky
    (:stddev-params
     (let* ((param-keys (walker-param-keys walker)))
       (if (< (walker-length walker) 10)
	   (make-plist param-keys (make-list (length param-keys) :initial-element 0d0))
	   (let ((l-matrix (walker-get walker :get :l-matrix :take take)))
	     (values (if (eq :multiple-kwargs (walker-param-style walker))
			 (make-plist param-keys
				     (mapcar (lambda (x)
					       (elt (elt l-matrix x) x))
					     (up-to (1- (length l-matrix)))))
			 (make-plist param-keys
				     (list (mapcar (lambda (x)
						     (elt (elt l-matrix x) x))
						   (up-to (1- (length l-matrix)))))))
		     l-matrix)))))
    (:log-liklihoods (mapcar #'walker-step-prob (walker-get walker :get :steps :take take)))
    (:covariance-matrix (lplist-covariance (walker-get walker :get :unique-steps :take take)))
    ;; TODO Experimental - forward steps vs unique steps
    (:l-matrix (lplist-to-l-matrix (diff-lplist (walker-get walker :get :forward-steps :take take))))))
(export 'walker-get)


(defun walker-modify (walker &key (modify (or :add-step :add-walks :burn-walks :keep-walks :reset :reset-to-most-likely :delete)) new-step new-walks burn-number keep-number)
  (case modify
    (:add-step (progn (push new-step (walker-walk walker))
		      (setf (walker-last-step walker) new-step)
		      (setf (walker-length walker) (1+ (walker-length walker)))
		      (setf (walker-age walker) (1+ (walker-age walker)))
		      (when (> (walker-step-prob new-step)
			       (walker-step-prob (walker-most-likely-step walker)))
			(setf (walker-most-likely-step walker) new-step))))
    (:add-walks (progn (nconc (reverse new-walks) (walker-walk walker))
		       (setf (walker-last-step walker) (last new-walks))
		       (setf (walker-length walker) (+ (length new-walks) (walker-length walker)))
		       (setf (walker-age walker) (+ (length new-walks) (walker-age walker)))
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
    (:reset-to-most-likely
     (progn (setf (walker-walk walker) (list (walker-most-likely-step walker)))
	    (setf (walker-last-step walker) (car (walker-walk walker)))
	    (setf (walker-length walker) 1)
	    walker))
    (:delete (progn (setf (walker-walk walker) nil)
		    (setf walker (make-walker))))))
(export 'walker-modify)

(defun cholesky-decomp (covariance-matrix)
  "Generalization of the sqrt operation for postive definite matrices. Takes a 2d, double-float array 'covariance-matrix' and turns it into a 2d, double-float array that is often referred to a an l-matrix (lower triangular matrix)."
  (declare (optimize speed (safety 1))
	   ((simple-array double-float (* *)) covariance-matrix))
  (let ((l-array (make-array (array-dimensions covariance-matrix) :element-type 'double-float))
	(tmp-sum 0d0))
    (declare (double-float tmp-sum))
    (dotimes (i (array-dimension covariance-matrix 0))
      (dotimes (k (+ i 1))
	(setq tmp-sum 0d0)
	(dotimes (j k)
	  (incf tmp-sum (* (aref l-array i j) (aref l-array k j))))
	(if (= i k)
	    (setf (aref l-array i k) (sqrt (max 0d0 (- (aref covariance-matrix i k) tmp-sum))))
	    (setf (aref l-array i k) (/ (- (aref covariance-matrix i k) tmp-sum) (aref l-array k k))))))
    l-array))

;;; matrix and covariance
(defun transpose-matrix (matrix)
  (let ((num-columns (- (length (elt matrix 0)) 1)))
    (mapcar (lambda (el) (mapcar (lambda (row) (elt row el)) matrix)) (up-to num-columns))))

(defun scale-array (scale matrix)
  (declare (optimize speed (safety 1))
	   ((simple-array double-float (* *)) matrix)
	   (double-float scale))
  (dotimes (i (array-dimension matrix 0) matrix)
    (dotimes (j (array-dimension matrix 1))
      (setf (aref matrix i j) (* scale (aref matrix i j))))))

;; TODO Takes a long time....
(defun lplist-covariance (lplist)
  "Takes in a list of p lists (lplist) that either have the structure ((:a 1 :b 2 :c 3 ...) (:a 1 :b 2 :c 3 ...) ...) or ((:all-params #(2 3 4 4 5 6 2 ...)) (:all-params #(2 3 4 4 5 6 2 ...)) ...). Returns a covariance matrix as a 2d, double-float array."
  (declare (optimize speed (safety 0))
	   (cons lplist))
  (let* ((num-of-points (length lplist))
	 (single-item? (consp (cadar lplist)))
	 (num-of-params (if single-item? (length (the list (cadar lplist))) (floor (length (the list (car lplist))) 2)))
	 (values (mapcar #'get-plist-values lplist))
	 (values (transpose (if single-item? (map 'list #'first values) values)))
	 ;; Sets of each like parameter in an array
	 (values-array (make-array (list num-of-params num-of-points) :element-type 'double-float :initial-contents values))
	 ;; Average value of each like parameter
	 (avg (map 'vector (lambda (x) (/ (the double-float (reduce #'+ x)) num-of-points)) (the list values)))
	 ;; Value of each parameter subtraced from its average
	 (average-normalized (make-array (list num-of-params num-of-points) :element-type 'double-float))
	 ;; covariance output (dot a with each other a and divide by n)
	 (covariance (make-array (list num-of-params num-of-params) :element-type 'double-float)))
    (declare ((simple-array double-float (* *)) values-array)
	     ((simple-vector *) avg))
    (dotimes (i num-of-params average-normalized)
      (dotimes (j num-of-points)
	(setf (aref average-normalized i j) (- (aref values-array i j) (the double-float (aref avg i))))))
    (dotimes (i num-of-params covariance)
      (dotimes (j num-of-params)
	(setf (aref covariance i j)
	      (do ((k 0 (1+ k))
		   (mini-sum 0d0))
		  ((>= k num-of-points) mini-sum)
		(declare (double-float mini-sum))
		(incf mini-sum (/ (* (aref average-normalized i k) (aref average-normalized j k)) num-of-points))))))))

(defun lplist-covariance-old (lplist)
  "Takes in a list of p lists (lplist) that either have the structure ((:a 1 :b 2 :c 3 ...) (:a 1 :b 2 :c 3 ...) ...) or ((:all-params #(2 3 4 4 5 6 2 ...)) (:all-params #(2 3 4 4 5 6 2 ...)) ...). Returns a covariance matrix as a 2d, double-float array."
  (declare (optimize speed (safety 0))
	   (cons lplist))
  (let* ((num-of-points (length lplist))
	 (single-item? (consp (cadar lplist)))
	 (num-of-params (if single-item? (length (the list (cadar lplist))) (floor (length (the list (car lplist))) 2)))
	 (values (mapcar #'get-plist-values lplist))
	 ;; Sets of each like parameter
	 (values (transpose (if single-item? (map 'list #'first values) values)))
	 ;; Average value of each like parameter
	 (avg (map 'list (lambda (x) (/ (the double-float (reduce #'+ x)) num-of-points)) values))
	 ;; Value of each parameter subtraced from its average
	 (average-normalized (make-array (list num-of-params num-of-points) :element-type 'double-float))
	 ;; covariance output (dot a with each other a and divide by n)
	 (covariance (make-array (list num-of-params num-of-params) :element-type 'double-float)))
    (dotimes (i num-of-params average-normalized)
      (dotimes (j num-of-points)
	(setf (aref average-normalized i j) (- (the double-float (elts values i j))
					       (the double-float (elt avg i))))))
    (dotimes (i num-of-params covariance)
      (dotimes (j num-of-params)
	(setf (aref covariance i j)
	      (do ((k 0 (1+ k))
		   (mini-sum 0d0))
		  ((>= k num-of-points) mini-sum)
		(declare (double-float mini-sum))
		(incf mini-sum (/ (* (aref average-normalized i k) (aref average-normalized j k)) num-of-points))))))))


(defun lplist-to-l-matrix (lplist)
  "Thin wrapper around Cholesky routine."
  (cholesky-decomp (lplist-covariance lplist)))

(defun get-covariant-sample (means l-matrix)
  "Rotates a sample by an l-matrix"
  (declare (optimize speed debug (safety 0))
	   ((simple-array double-float (* *)) l-matrix means))
  (let* ((imax (array-dimension l-matrix 0))
	 (jmax (array-dimension l-matrix 1))
	 (kmax (array-dimension means 1))
	 (return-array (make-array (list imax kmax) :element-type 'double-float))
	 (samples-array (make-array (list (array-dimension means 0) 1) :element-type 'double-float :initial-contents (mapcar #'list (repeat (array-dimension means 0) (alexandria:gaussian-random))))))
    (declare ((simple-array double-float *) samples-array return-array))
    ;; matrix multiply ==> l-matrix·means - in that order
    (dotimes (i imax)
      (dotimes (k kmax)
	(setf (aref return-array i k)
	      (do ((j 0 (1+ j))
		   (mini-sum 0d0))
		  ((>= j jmax) mini-sum)
		(declare (double-float mini-sum))
		(incf mini-sum (* (aref l-matrix i j) (aref samples-array j k)))))))
    ;; Add previous matrix (shift matrix) to the provided means
    (dotimes (i imax return-array)
      (setf (aref return-array i 0) (+ (aref return-array i 0) (aref means i 0))))))

;; (mfit::get-covariant-sample #(5 9 2) (mfit::lplist-to-l-matrix mfit::example-array-lplist))
;; (defun get-covariant-sample (means l-matrix)
;;   "Takes a 2d, double-float array of previous means, a 2d, double-float array l-matrix. Takes a random sample with the same length as the means, rotates it by the l-matrix, and adds it to the previous means. Returns a 2d, double-float array (column vector)."
;;   (let ((samples (repeat (array-dimension means 0) (alexandria:gaussian-random))))
;;     (if (consp (car means))
;; 	(list (mapcar #'+ (car means) (rotate-random-sample samples l-matrix)))
;; 	(mapcar #'+ means (rotate-random-sample samples l-matrix)))))

(defun diagonal-covariance (set)
  "Take a set of means (list, vector, or column vector (2d array)) and return a 2d array with the means on the diagonal."
  (let* ((len (typecase (car set)
		(number (length set))
		((or list vector) (length (car set)))
		(array (array-dimension (car set) 0))
		(otherwise 0)))
	 (return-array (make-array (list len len) :element-type 'double-float)))
    (dotimes (i len return-array)
      (dotimes (j len)
	(setf (aref return-array i j)
	      (if (= i j)
		  (typecase (car set)
		    (number (elt set i))
		    ((or list vector) (elt (car set) i))
		    (array (aref (car set) i 0))
		    (otherwise 0d0))
		  0d0))))))

(defvar example-lplist '((:a 90d0 :b 60d0 :c 90d0)
			 (:a 90d0 :b 90d0 :c 30d0)
			 (:a 60d0 :b 60d0 :c 60d0)
			 (:a 60d0 :b 60d0 :c 90d0)
			 (:a 30d0 :b 30d0 :c 30d0)))
(defvar example-array-lplist
  (mapcar (lambda (x)
	    (list :params (make-array '(3) :element-type 'double-float :initial-contents x)))
	  '((90d0 60d0 90d0)
	    (90d0 90d0 30d0)
	    (60d0 60d0 60d0)
	    (60d0 60d0 90d0)
	    (30d0 30d0 30d0))))
(defparameter example-covariance (list (lplist-covariance example-lplist)
				       ;;(lplist-covariance example-array-lplist)
				       ))
;; #2A((504.0d0 360.0d0 180.0d0) (360.0d0 360.0d0 0.0d0) (180.0d0 0.0d0 720.0d0))
(defparameter example-l-matrix (list (lplist-to-l-matrix example-lplist)
				     ;;(lplist-to-l-matrix example-array-lplist)
				     ))
;; #2A((22.44994432064365d0 0.0d0 0.0d0)
;;     (16.035674514745462d0 10.141851056742201d0 0.0d0)
;;      (8.017837257372731d0 -12.677313820927745d0 22.248595461286993d0))


;;; MCMC Metropolis-Hastings
(defun force-list (item)
  "Forces all non-list items into lists"
  (if (consp item)
      item
      (list item)))

(defun get-depth (list-tree)
  "Returns the depth of the first element (as defined by a flattened sequence) and the length of the sequence containing that element. (list max-depth deepest-length).
==> (get-depth '(((4 5) (3 4))))
==> (3 2)"
  (cond ((null list-tree)
	 nil)
	((numberp list-tree)
	 0)
	((arrayp list-tree)
	 (array-rank list-tree))
	(t
	 (+ 1 (get-depth (elt list-tree 0))))))

(defun clean-data-error (stddev clean-data &optional (first t))
  "Copies the structure of the y-data and inserts provided sttdev value(s) into their place if they fit. If the stddev structure does not match the y-data structure, then just replace all numbers found in y with the first from stddev.
==> (clean-data-error '((0.1) (0.4)) '(((2 3) ((4) (5))) (() (3 4 5 6))))
==> (((0.1) (0.4)) (0.1 0.1 0.1 0.1))"
  (labels ((first-element (list-tree)
	     (cond ((null list-tree)
		    nil)
		   ((numberp list-tree)
		    list-tree)
		   (t
		    (first-element (car list-tree)))))
	   (eq-structure (list-tree1 list-tree2)
	     (cond ((and (numberp list-tree1) (numberp list-tree2))
		    t)
		   ((or (numberp list-tree1) (numberp list-tree2))
		    nil)
		   ((/= (length list-tree1) (length list-tree2))
		    nil)
		   (t
		    (every #'identity (mapcar #'eq-structure list-tree1 list-tree2))))))
    (let ((clean-data-ys (if first (mapcar #'second clean-data) clean-data))
	  (default-stddev (first-element stddev)))
      (cond ((null clean-data-ys)
	     nil)
	    ((numberp clean-data-ys)
	     default-stddev)
	    ((eq-structure stddev clean-data-ys)
	     stddev)
	    ((arrayp clean-data-ys)
	     (make-array (array-dimensions clean-data-ys) :element-type 'double-float :initial-element default-stddev))
	    (t
	     (mapcar (lambda (data) (clean-data-error stddev data nil)) clean-data-ys))))))

(defun clean-data (data number-of-functions)
  "Forces data to be lists only and of proper depth - unless they are arrays."
  (labels ((list-everything (list-tree)
	     (cond ((null list-tree)
		    nil)
		   ((numberp list-tree)
		    list-tree)
		   ((arrayp list-tree)
		    list-tree)
		   ((or (listp list-tree) (vectorp list-tree))
		    (map 'list #'list-everything list-tree)))))
    (cond ((= (get-depth data) 1)
	   (error "clean-data: data is of insufficient depth or improperly structured."))
	  ((= (get-depth data) 2)
	   (clean-data (list data) number-of-functions))
	  ((= (length data) number-of-functions)
	   (list-everything data))
	  (t
	   (error "clean-data: insufficient number of datasets, ~a, for the given number of functions, ~a." (length data) number-of-functions)))))

(defun create-walker-data (data &rest columns)
  "Takes a larger dataset (contains x, y, stddev, freq, temp, etc.) and extracts just the columns you want into a walker-friendly format."
  (mapcar (lambda (x) (apply #'vector x))
	  (mapcar (lambda (y) (elt data y)) columns)))
(export 'create-walker-data)

(defun to-double-floats (l)
  "Converts all numbers and vectors in tree to double-float type"
  (map-tree (lambda (x) (etypecase x (number (coerce x 'double-float)) (symbol x) ((array double-float) x))) l))

(defun log-prior-fixer (log-prior params data)
  "Checks if your prior returns a function (i.e. changes shape based on the data provided to it). If so, gets that function. If not, returns the prior you gave it."
  (let ((results (mapcar (lambda (lp d) (funcall lp params d)) log-prior data)))
    (mapcar (lambda (res fn) (if (numberp res) fn res)) results log-prior)))

(defun log-liklihood-fixer (log-liklihood fn params data error)
  "Checks if your liklihood returns a function (i.e. changes shape based on the data provided to it). If so, gets that function. If not, returns the liklihood you gave it."
  (let ((results (mapcar (lambda (ll f d e) (funcall ll f params d e)) log-liklihood fn data error)))
    (mapcar (lambda (res fn) (if (numberp res) fn res)) results log-liklihood)))



(defun walker-many-steps (the-walker n &optional l-matrix)
  "Takes 'n' steps with a constant l-matrix. No temperature or any other intelligent features."
  (if (null l-matrix) (setq l-matrix (diagonal-covariance (get-plist-values (scale-plist 1e-2 (walker-get the-walker :get :median-params))))))
  (dotimes (i n)
    (walker-take-step the-walker :l-matrix l-matrix)))

;; TODO remove and add functionality for troubleshooting
;; (print (car (diagonal-covariance (get-plist-values (scale-plist 1e-2 (walker-get-median-params woi 1000))))))
;; (print (standard-deviation (diff (remove-consecutive-duplicates (walker-get-param woi :w1-0-0)))))
;; (plt:plot (diff (remove-consecutive-duplicates (walker-get-param woi :w1-0-0))))

(defvar mfit-walker-estop nil)
(export 'mfit-walker-estop)
(defun walker-adaptive-steps-full (walker &key (n 100000) (temperature 1d3) (auto (or :prob-settle :slope-settle nil)) (sampling-optimization (or :covariance :best-value)) max-walker-length l-matrix)
  "Advances 'walker' through as many as 'n' steps. Intelligently adjusts the l-matrix (sampling rotation matrix for more efficient random sampling) to try to avoid too low and too high acceptance. Has an oscillating temperature between 1 and 'temperature' during the beginning of the journey. Can automatically terminate a journey if the probability has stabilized, indicating at least a point of local stability. 'max-walker-length' will prevent the walker from becomming too large if maxing out memory is an issue."
  (declare (optimize debug))
  (setf mfit-walker-estop nil)
  (let* ((n (floor n))
	 (reset-index 10000)
	 (max-walker-length (when max-walker-length (floor max-walker-length 2)))
	 (num-params (etypecase (cadr (walker-step-params (walker-most-likely-step walker)))
		       (number (length (walker-param-keys walker)))
		       ((or list vector) (length (cadr (walker-step-params (walker-most-likely-step walker)))))
		       (array (array-dimension (cadr (walker-step-params (walker-most-likely-step walker))) 0))))
	 (steps-to-settle (* 10 (max 50 num-params)))
	 (shutting-down-p nil)
	 (temp-steps (max n (* 10 steps-to-settle)))
	 (temperature (rational temperature))
	 ;; Temperature starts high and cycles every 5000 steps (the divisor to floor function)
	 (temps (mapcar (lambda (x) (max 1 (* (cos (* x pi (+ 1 (* 2 (floor temp-steps 5000))) (/ (* 2 temp-steps)))) temperature))) (range temp-steps)))
	 (current-acceptance 0))
    (flet ((stable-probs-p (probs)
	     (let* ((num-probs (length probs))
		    (early-max (reduce #'max (subseq probs 0 200)))
		    (late-max (reduce #'max (subseq probs (- num-probs 200)))))
	       (and (< (abs (- early-max late-max)) 0.5) ;; stable max values
		    (< 4 (- early-max (reduce #'min probs)) 9)))) ;; good variation
	   (stable-prob-slope-p (probs)
	     (< (getf (walker-step-params (walker-most-likely-step (mcmc-fit :function (lambda (x &key m b &allow-other-keys) (+ b (* b m x (/ (length probs))))) :data (list (list (thin (range (length probs)) 10) (thin probs 10))) :data-error 1d0 :params (list :m 10 :b -100)))) :m) 1)) ;; Is the log-liklihood slope less than 1%
	   (get-optimal-mcmc-l-matrix (take)
	     (if (eq sampling-optimization :covariance)
		 (scale-array (/ (expt 2.38d0 2) num-params)
			      (handler-case (walker-get walker :get :l-matrix :take take)
				(floating-point-overflow () (print "FPO in Cholesky") l-matrix)
				(division-by-zero () (print "DB0 is Cholesky") l-matrix)
				(type-error () (print "TE in Cholesky") l-matrix)))
		 (scale-array 1d-5 (diagonal-covariance (get-plist-values (walker-step-params (walker-most-likely-step walker))))))))
      (unless l-matrix
	(if (or (< (walker-length walker) steps-to-settle)
		(< (walker-get walker :get :acceptance :take 100) 0.1))
	    (setq l-matrix (diagonal-covariance (get-plist-values (walker-step-params (walker-most-likely-step walker)))))
	    (progn (setq l-matrix (diagonal-covariance (get-plist-values (walker-step-params (walker-most-likely-step walker)))))
		   (setq l-matrix (get-optimal-mcmc-l-matrix steps-to-settle)))))
      (do ((i 1 (+ i 1))
	   (temp-index 0))
	  ((or (>= i n) mfit-walker-estop))
	(when (or
	       (and (not shutting-down-p) (< (- n i) (max 2000 steps-to-settle)))
	       (and auto
		    (not shutting-down-p)
		    (= 0 (mod i 1000))
		    (> i (* 2 steps-to-settle))
		    (< 0.2 (walker-get walker :get :acceptance :take 1000) 0.5)
		    (or
		     (and (eq auto :prob-settle) (stable-probs-p (walker-get walker :get :log-liklihoods :take steps-to-settle)))
		     (and (eq auto :slope-settle) (stable-prob-slope-p (walker-get walker :get :log-liklihoods :take (min (walker-length walker) (max 2500 steps-to-settle))))))))
	  (setf temperature 1)
	  (setf shutting-down-p t)
	  (setf i (- n (max 2000 steps-to-settle))))
	(walker-take-step walker :l-matrix l-matrix :temperature temperature)
	;; Annealing
	(when (and (not shutting-down-p) (< i temp-steps))
	  (setf temperature (elt temps i)))
	;; Cleaning the walker
	(when (and max-walker-length (= i reset-index))
	  (if (> (walker-length walker) max-walker-length)
	      (progn (walker-modify walker :modify :keep-walks :keep-number max-walker-length)
		     (incf reset-index (1+ max-walker-length)))
	      (incf reset-index (1+ (- max-walker-length (walker-length walker))))))
	;; Regular l-matrix updating
	(when (and (> i 0) 
		   (or (and (= 0 (mod i 200)) (< (walker-get walker :get :acceptance :take 200) 0.2))
		       (and (= 0 (mod i 200)) (> (walker-get walker :get :acceptance :take 200) 0.4))
		       (and (not shutting-down-p) (= 0 (mod i (* 2 steps-to-settle))))))
	  (setf current-acceptance (walker-get walker :get :acceptance :take 200))
	  (if (< 0.2 current-acceptance 0.4)
	      (setf l-matrix (let ((tmp-l-matrix (get-optimal-mcmc-l-matrix steps-to-settle)))
			       (if (= num-params (array-dimension tmp-l-matrix 0))
				   tmp-l-matrix
				   l-matrix)))
	      (if (< current-acceptance 0.2)
		  (setf l-matrix (scale-array 0.1d0 l-matrix))
		  (if (> current-acceptance 0.4)
		      (setf l-matrix (scale-array 1.9d0 l-matrix))))))))))

;; (let ((temp-steps 100)) (mapcar (lambda (i) (and (= 0 (mod (1+ (floor i temp-steps)) 2)) (< i (* temp-steps 95)))) (range 500)))

(defun walker-adaptive-steps (walker &optional (n 30000))
  (walker-adaptive-steps-full walker :n n :temperature 10 :auto :prob-settle))

(defun walker-sample-region (walker &optional (initial-scale 1d-3))
  "Advances 'walker' through as many as 'n' steps. Intelligently adjusts the l-matrix (sampling rotation matrix for more efficient random sampling) to try to avoid too low and too high acceptance. Has an oscillating temperature between 1 and 'temperature' during the beginning of the journey. Can automatically terminate a journey if the probability has stabilized, indicating at least a point of local stability. 'max-walker-length' will prevent the walker from becomming too large if maxing out memory is an issue."
  (declare (optimize debug))
  (setf mfit-walker-estop nil)
  (let* ((num-params (etypecase (cadr (walker-step-params (walker-most-likely-step walker)))
		       (number (length (walker-param-keys walker)))
		       ((or list vector) (length (cadr (walker-step-params (walker-most-likely-step walker)))))
		       (array (array-dimension (cadr (walker-step-params (walker-most-likely-step walker))) 0))))
	 (steps-to-settle 3000)
	 ;; Temperature starts high and cycles every 5000 steps (the divisor to floor function)
	 (l-matrix (scale-array initial-scale (diagonal-covariance (get-plist-values (walker-step-params (walker-most-likely-step walker)))))))
    (do ((i 1 (+ i 1))
	 (temp-index 0))
	((or (>= i steps-to-settle) mfit-walker-estop))
      ;; catch complex stddev and retry with different opimal l matrix sampling
      (walker-pretend-take-step walker :l-matrix l-matrix)
      (when (and (> i 0) (= 0 (mod i 20)))
	(if (= (walker-get walker :get :acceptance :take 50) 1/50)
	    (setf l-matrix (scale-array 0.25d0 l-matrix))
	    (if (> (walker-get walker :get :acceptance :take 50) 4/50)
		(setf l-matrix (scale-array 1.7d0 l-matrix))))))))

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

(defun walker-set-get (the-walker-set &key (get (or :steps :unique-steps :most-likely-step :acceptance :param :most-likely-params :median-params :all-params :stddev-params :log-liklihoods :covariance-matrix :l-matrix)) take param)
  (mapcar (lambda (x) (walker-get x :get get :take take :param param)) the-walker-set))

(defun walker-set-delete (the-walker-set)
  (mapcar (lambda (x) (walker-modify x :modify :delete)) the-walker-set))

(defun walker-set-plot-param (the-walker-set key &optional take)
  (plt:plot (walker-set-get the-walker-set :get :median-params :param key :take take) "w l title \"Param\""))

;; median vs most likely
(defmacro walker-get-f (walker exp &key (take 1000))
  (let ((median-params (walker-get walker :get :most-likely-params :take take)))
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
  (let ((median-params (walker-get walker :get :most-likely-params :take take)))
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
  (if (null l-matrix) (setq l-matrix (diagonal-covariance (plist-values (scale-plist 1e-2 (walker-get walker :get :most-likely-params :take 1000))))))
  (let* ((previous-step (walker-last-step walker))
	 (prob0 (walker-step-prob previous-step))
	 (previous-params (walker-step-params previous-step))
	 (new-values (get-covariant-sample (etypecase (cadr previous-params)
					     (number (make-array (list (length (plist-keys previous-params)) 1) :element-type 'double-float :initial-contents (mapcar #'list (plist-values previous-params))))
					     (list (make-array (list (length (cadr previous-params)) 1) :element-type 'double-float :initial-contents (mapcar #'list (cadr previous-params))))
					     (vector (make-array (list (length (cadr previous-params)) 1) :element-type 'double-float :initial-contents (map 'list #'list (cadr previous-params))))
					     (array (cadr previous-params)))
					   l-matrix))
	 (next-params (etypecase (cadr previous-params)
			(number (array-to-plist (plist-keys previous-params) new-values))
			(list (list (car previous-params) (mapcar (lambda (x) (aref new-values x 0)) (range (array-dimension new-values 0)))))
			(vector (list (car previous-params) (map 'vector (lambda (x) (aref new-values x 0)) (range (array-dimension new-values 0)))))
			(array (list (car previous-params) new-values))))
	 (next-step (walker-make-step walker next-params))
	 (prob1 (walker-step-prob next-step))
	 (accepted-step (if (or (> prob1 prob0)
				(> (/ (- prob1 prob0) temperature) (log (random 1.0d0))))
			    next-step
			    previous-step)))
    (walker-modify walker :modify :add-step :new-step accepted-step)))

(defun walker-pretend-take-step (walker &key l-matrix)
  "Normal step taking routine, but always puts the previous step at the last step, so the walker doesn't actually move. Used for sampling the current location without moving the walker for proper l-matrix construction."
  (let ((previous-step (walker-last-step walker)))
    (if (null l-matrix) (setq l-matrix (diagonal-covariance (plist-values (scale-plist 1e-2 (walker-get walker :get :most-likely-params :take 1000))))))
  (let* ((previous-step (walker-last-step walker))
	 (prob0 (walker-step-prob previous-step))
	 (previous-params (walker-step-params previous-step))
	 (new-values (get-covariant-sample (etypecase (cadr previous-params)
					     (number (make-array (list (length (plist-keys previous-params)) 1) :element-type 'double-float :initial-contents (mapcar #'list (plist-values previous-params))))
					     (list (make-array (list (length (cadr previous-params)) 1) :element-type 'double-float :initial-contents (mapcar #'list (cadr previous-params))))
					     (vector (make-array (list (length (cadr previous-params)) 1) :element-type 'double-float :initial-contents (map 'list #'list (cadr previous-params))))
					     (array (cadr previous-params)))
					   l-matrix))
	 (next-params (etypecase (cadr previous-params)
			(number (array-to-plist (plist-keys previous-params) new-values))
			(list (list (car previous-params) (mapcar (lambda (x) (aref new-values x 0)) (range (array-dimension new-values 0)))))
			(vector (list (car previous-params) (map 'vector (lambda (x) (aref new-values x 0)) (range (array-dimension new-values 0)))))
			(array (list (car previous-params) new-values))))
	 (next-step (walker-make-step walker next-params))
	 (prob1 (walker-step-prob next-step))
	 (accepted-step (if (> prob1 prob0)
			    next-step
			    previous-step)))
    (walker-modify walker :modify :add-step :new-step accepted-step))
    ;;(walker-modify walker :modify :add-step :new-step previous-step)
    ))

(defun walker-force-take-step (walker)
  "Forces a new step onto the walker regardless of probability. Does not take a random sample, just uses the current param value. Useful for when swapping datasets within a walker."
  (let* ((previous-step (walker-last-step walker))
	 (previous-params (walker-step-params previous-step))
	 (next-step (walker-make-step walker previous-params)))
    (walker-modify walker :modify :add-step :new-step next-step)))

;; TODO prior bounds implementation - it's a good idea
(defun walker-create (&key function data params data-error log-liklihood log-prior param-bounds)
  "Make a walker with the prescribed characteristics. Can run over multiple sets of functions and data by simply making lists of kwargs. i.e. :function #'my-function vs :function (list #'my-function #'other-function). data, params, etc. must also then be lists that correspond to the order found in the :function kwarg.
:function expects functions formatted with kwargs as parameters:
(lambda (x &key m b &allow-other-keys) (+ b (* -3 m) (* (- m (/ b 60)) x)))
or like this for multiple or linked independent variables
(lambda (x &key m b &allow-other-keys) (+ b (* -3 m (elt x 0)) (* (- m (/ b 60)) (elt x 1))))
:data expects formatting as (list x-data y-data)
:data-error can be a single number (for uniform error) or a list of the same size as the y-data
:params expects a plist like '(:b -1 :m 2)"
  (declare (ignorable param-bounds))
  (let* ((function (force-list function))
	 (data (to-double-floats (clean-data data (length function))))
	 (data-error (to-double-floats (clean-data-error (if data-error data-error 1) data)))
	 (params (to-double-floats params))
	 (log-liklihood (log-liklihood-fixer (if (consp log-liklihood) log-liklihood (make-list (length (force-list function)) :initial-element (if log-liklihood log-liklihood #'log-liklihood-normal))) function params data data-error))
	 (log-prior (log-prior-fixer (if (consp log-prior) log-prior (make-list (length (force-list function)) :initial-element (if log-prior log-prior #'log-prior-flat))) params data))
	 (first-step (make-walker-step :prob (+ (reduce #'+ (map 'list (lambda (ll f d e) (funcall ll f params d e)) log-liklihood function data data-error))
						(reduce #'+ (map 'list (lambda (lp d) (funcall lp params d)) log-prior data)))
				       :params params))
	 (walker (make-walker :function function
			      :param-keys (plist-keys params)
			      :param-style (typecase (elt params 1)
					     ((or list vector array) :single-item)
					     (otherwise :mutliple-kwargs))
			      :last-step first-step
			      :most-likely-step first-step
			      :walk (list first-step)
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

;; Basic line fit
;; (defparameter woi (mfit:mcmc-fit :function (lambda (x &key m b &allow-other-keys) (+ b (* m x))) :data '((-4 -1 2 5 10) (0 2 5 9 13)) :params '(:b -1 :m 2) :data-error 0.2))

;; Single param list fit
;; (defparameter woi (mfit:mcmc-fit :function (lambda (x &key params &allow-other-keys) (+ (elt params 0) (* (elt params 1) x))) :data '((-4 -1 2 5 10) (0 2 5 9 13)) :params '(:params (-1 2)) :data-error 0.2))

;; Single param vector fit
;; (defparameter woi (mfit:mcmc-fit :function (lambda (x &key params &allow-other-keys) (+ (elt params 0) (* (elt params 1) x))) :data '((-4 -1 2 5 10) (0 2 5 9 13)) :params '(:params #(-1 2)) :data-error 0.2))

;; Single param array fit
;; (defparameter woi (mfit:mcmc-fit :function (lambda (x &key params &allow-other-keys) (+ (aref params 0 0) (* (aref params 1 0) x))) :data '((-4 -1 2 5 10) (0 2 5 9 13)) :params (list :params (make-array '(2 1) :element-type 'double-float :initial-contents '((-1d0) (2d0)))) :data-error 0.2))

;; Global parameter fit (polynomial and line)
;; (defparameter woi (mfit:mcmc-fit :function (list (lambda (x &key b m c d &allow-other-keys) (+ b (* m x) (* c x x) (* d x x x))) (lambda (x &key e m g &allow-other-keys) (+ e (* (+ m g) x)))) :data '(((0 1 2 3 4) (4 5 6 8 4)) ((10 20 30 40 50) (1 2 3 4 5))) :params '(:b -1 :m 2 :c 1 :d -1 :e 0.5 :g -2) :data-error 0.2))

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
    (plt:plot :terminal (apply #'format nil "qt size ~d,~d linewidth 3 font \"Arial,18\"" image-size)
	      :format "x \"%.4g\""
	      :format "y \"%.4g\""
	      :xlabel "\"x-data\""
	      :ylabel "\"y-data\""
	      x-fit max-ys "w l lc rgb \"green\" title \"fit stddev upper limit\""
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
    (plt:plot  :terminal "qt size 1920,1080 linewidth 3 font \"Arial,18\""
	       :format "x \"%.4g\""
	       :format "y \"%.4g\""
	       :xlabel "\"x-data\""
	       :ylabel "\"y-data\""
	       x-data y-residuals "w p pt 6 ps 2 lc rgb \"black\" title \"residuals\""
	       x-data stddev "w p pt 2 ps 1 lc rgb \"red\" title \"point error\""
	       x-data (make-array (length x-data) :initial-element 0) "w l lc rgb \"red\" title \"baseline\"")))

(defun walker-catepillar-plots (walker &optional take (x-scale 1) (y-scale 1))
  (let* ((param-keys (walker-param-keys walker))
	 (n (length param-keys)))
    (plt:reset)
    (plt:send-plot-options :terminal (format nil "pngcairo size ~d,~d linewidth 1 font \"Arial,18\"" (* x-scale 1920) (* y-scale 1080))
			   :output "\"temp.png\""
			   :format "y \"%.1e\""
			   :format "x \"%.1e\"")
    (apply #'plt:multiplot
	   (format nil "layout ~d,1" n)
	   (apply #'append
		  (mapcar-enum (key i param-keys)
		    (list (append (list (walker-get walker :get :param :param key :take take) (format nil "w l title \"~a\"" key))
				  (if (= i (1- n))
				      (list :xlabel "\"Step\"" :ylabel (format nil "\"~a\"" key))
				      (list :xlabel "" :ylabel (format nil "\"~a\"" key))))))))))
(export 'walker-catepillar-plots)
(export 'show)

(defun walker-liklihood-plot (walker &optional take)
  (let ((probs (walker-get walker :get :log-liklihoods :take take)))
    (plt:reset)
    (plt:plot :terminal "qt size 1920,1080 linewidth 4 font \"Arial,40\""
	      :xtics ""
	      :xlabel "\"Step\""
	      :ylabel "\"Log Liklihood\""
	      probs "w l title \"Liklihood\"")))

(defun walker-plot-one-corner (walker keys &optional take)
  (let* ((x-key (first keys))
	 (y-key (second keys))
	 (x (walker-get walker :get :param :param x-key :take take))
	 (y (walker-get walker :get :param :param y-key :take take)))
    (plt:reset)
    (plt:plot :terminal "qt size 1920,1080 linewidth 1 font \"Arial,18\""
	      :xlabel (format nil "\"~a\"" x-key)
	      :ylabel (format nil "\"~a\"" y-key)
	      x y "w p pt 3 ps 4 title \"My plot\"")))

(defun walker-plot-corner (walker &key take (size '(3240 3240)))
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
    (plt:plot :terminal "qt size 1920,1080 linewidth 3 font \"Arial,18\""
	      :xlabel (format nil "\"~a\"" key)
	      :ylabel "\"Counts\""
	      histo-x histo (format nil "w lp ps 3 title \"Histogram: ~a\"" key))))

(defun show ()
  "Shows the plots generated from 2d histogram and catepillar functions"
  (uiop:run-program "feh ./temp.png -3 5"))


;;; File Management
;; File searching
(defun walk-dirs (dir)
  (if (uiop:directory-exists-p dir) ;; is it a dir
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
  "Read 'filename' line by line and return a list of strings"
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
     (map-tree
      (lambda (x) (read-from-string x nil nil))
      (transpose
       (mapcar
	(lambda (x) (split-string most-likely-delim x))
	unsplit-data-lines))))))
(export 'auto-split-and-read-csv)

(defun file->file-specs (filename &key (delim #\tab))
  "Returns various metrics about the file. Lines, header lines, data lines, data rows, data sets."
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
  "Reads a file into a list assuming 'delim' separates data. Use 'pages' kwargs to make 3d list assuming each page of data is separated by an extra newline in the file."
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

(export '(log-normal prior-bounds-let standard-deviation walker-create log-liklihood-normal read-file->data walker-adaptive-steps walker-adaptive-steps-full walker-load log-prior-flat walker-liklihood-plot walker-plot-data-and-fit linspace walker-set-get get-filename))
