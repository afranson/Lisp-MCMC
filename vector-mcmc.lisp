;;; Trying to figure out how to do some basic data analysis in lisp
;;; Need to:
;;;  Read data into array - read x lines in info, the rest into data, handle empty lines
;;;  Manipulate arrays to get particular columns, parts of columns
;;;  Plot xy data (vgplot - wrapper over gnuplot session)
;;;  Fit xy data (bayesian)
;;;  Plot fits

;;; TODO add a frame to the analysis - this gets added to log liklihood *
;;; - this frame moves around randomly so that N datapoints get fit to
;;; - trying to emulate more human problem solving and reduce the scope of the problem
;;; - that way all the noise from the surrounding background doesn't get incorporated

;; Full arrays > list of arrays > pure lists
;; just lists don't benefit from optimization
;; arrays do a lot
;; arrays aren't much better than lists without optimization

;; (ql:quickload :gsll)
;; (ql:quickload :antik)

(in-package :mcmc-fitting)

;;; data manipulation
(defun diff-list (x-data y-data)
  (labels ((avg-x-values (x-data)
	     (/ (+ (car x-data) (cadr x-data)) 2))
	   (diff-y (x-data y-data)
	     (if (equal (cadr x-data) (car x-data))
		 0
		 (/ (- (cadr y-data) (car y-data)) (- (cadr x-data) (car x-data)))))
	   (build-diff (x-data y-data &optional return-x return-y)
	     (if (and (cadr x-data) (cadr y-data))
		 (build-diff (cdr x-data) (cdr y-data)
			     (nconc (list (avg-x-values x-data)) return-x)
			     (nconc (list (diff-y x-data y-data)) return-y))
		 (list (reverse return-x) (reverse return-y)))))
    (build-diff x-data y-data)))

(defun diff-vec (x-data y-data)
  (let* ((len (array-dimension x-data 0))
	 (return-x-vec (make-array (list (- len 1)) :fill-pointer 0))
	 (return-y-vec (make-array (list (- len 1)) :fill-pointer 0)))
      (labels ((avg-x-values (n x-data)
		 (/ (+ (aref x-data n) (aref x-data (+ n 1))) 2))
	       (diff-y (n x-data y-data)
		 (/ (- (aref y-data (+ n 1)) (aref y-data n))
		    (- (aref x-data (+ n 1)) (aref x-data n))))
	       (build-diff (&optional (n 0))
		 (when (< n (- len 1))
		   (vector-push (avg-x-values n x-data) return-x-vec)
		   (vector-push (diff-y n x-data y-data) return-y-vec)
		   (build-diff (+ n 1)))))
	(build-diff)
	(list return-x-vec return-y-vec))))


;;; Stats
(defun multivariate-gaussian-random (covs)
  (mapcar #'(lambda (x) (* x (alexandria:gaussian-random))) covs))

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
	 (sum-sq (reduce #'+ (map 'list #'(lambda (x) (expt (- x mean) 2)) sequence))))
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

(defun cholesky-decomp (covariance-matrix)
  (let* ((cov-array (make-array (list (length covariance-matrix) (length (car covariance-matrix))) :initial-contents covariance-matrix))
	 (l-array (make-array (array-dimensions cov-array) :initial-element 0d0))
	 (tmp-sum 0d0))
    (dotimes (i (array-dimension cov-array 0))
      (dotimes (k (+ i 1))
	(setq tmp-sum 0d0)
	(dotimes (j k)
	  (incf tmp-sum (* (aref l-array i j) (aref l-array k j))))
	(if (= i k)
	    (setf (aref l-array i k) (sqrt (abs (- (aref cov-array i k) tmp-sum)))) ; TODO remove abs
	    (setf (aref l-array i k) (/ (- (aref cov-array i k) tmp-sum) (aref l-array k k))))))
    (array-matrix l-array)))




;;; plist functions
(defun get-plist-keys (plist &optional return-keys)
  (if (car plist)
      (get-plist-keys (cddr plist) (cons (car plist) return-keys))
      (reverse return-keys)))

(defun get-plist-values (plist &optional return-values)
  (if (car plist)
      (get-plist-values (cddr plist) (cons (cadr plist) return-values))
      (reverse return-values)))

(defun make-plist (keys values &optional return-plist)
  (if (and (nthcdr 0 keys) (nthcdr 0 values))
      (make-plist (cdr keys) (cdr values) (nconc (list (car values) (car keys)) return-plist))
      (reverse return-plist)))

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
  (map-plist #'(lambda (x) (* x scale)) plist))





;;; matrix and covariance
;; An array is a lisp array. A matrix is a list of lists
(defun dot (list1 list2)
  (reduce #'+ (mapcar #'* list1 list2)))

(defun scale-matrix (scale matrix)
  (mapcar #'(lambda (x) (mapcar #'(lambda (y) (* y scale)) x)) matrix))

(defun lplist-covariance (lplist)
  (let* ((n (length lplist))
	 (keys (get-plist-keys (elt lplist 0)))
	 (values (mapcar #'(lambda (y) (mapcar #'(lambda (x) (getf x y)) lplist)) keys))
	 (avg (mapcar #'(lambda (x) (/ (reduce #'+ x) n)) values))
	 (a (mapcar #'(lambda (x y)
			(mapcar #'(lambda (z) (- z y))
				x))
		    values avg))
	 (ns (make-list (length a) :initial-element (float n 0d0)))
	 (v (mapcar #'(lambda (x denom)
			(mapcar #'(lambda (y)
				    (/ (dot x y) denom))
				a))
		    a ns)))
    v))

(defun map-array-row (fn row-n array)
  (let ((dim1 (array-dimension array 1))
	(i -1))
    (repeat dim1 (progn
		   (incf i)
		   (setf (aref array row-n i) (funcall fn (aref array row-n i)))))
    array))

(defun map-array (fn array)
  (let ((dim0 (array-dimension array 0))
	(i -1))
    (repeat dim0 (progn
		   (incf i)
		   (map-array-row fn i array)))
    array))

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

(defun map-matrix (fn matrix)
  (mapcar #'(lambda (x) (mapcar #'(lambda (y) (funcall fn y)) x)) matrix))

(defun transpose-matrix (matrix)
  (let ((num-columns (- (length (elt matrix 0)) 1)))
    (mapcar #'(lambda (el) (mapcar #'(lambda (row) (elt row el)) matrix)) (linspace 0 num-columns))))

(defun lplist-to-l-matrix (lplist)
  (let* ((covariance (lplist-covariance lplist)))
    (cholesky-decomp covariance)))

(defun rotate-random-sample (sample l-matrix)
  (mapcar #'(lambda (y) (dot y sample)) l-matrix))

(defun get-covariant-sample (means l-matrix)
  (let ((samples (repeat (length l-matrix) (alexandria:gaussian-random))))
    (mapcar #'+ means (rotate-random-sample samples l-matrix))))

(defun diagonal-covariance (list)
  (let ((len (length list))
	(i -1))
    (mapcar #'(lambda (x)
		(declare (ignore x))
		(incf i)
		(maplist #'(lambda (y) (if (= (length y) (- len i)) (car y) 0)) list))
	    list)))

(defvar example-lplist '((:a 90 :b 60 :c 90)
			 (:a 90 :b 90 :c 30)
			 (:a 60 :b 60 :c 60)
			 (:a 60 :b 60 :c 90)
			 (:a 30 :b 30 :c 30)))
(defvar example-covariance (lplist-covariance example-lplist))
(defvar example-array (make-array '(3 3) :initial-contents example-covariance))
(defvar example-matrix (array-matrix example-array))
(defvar example-l-matrix (lplist-to-l-matrix example-lplist))





;;; binning / histogram
(defun make-histo (sequence &optional num-bins)
  (let* ((bottom (reduce #'min sequence))
	 (top (reduce #'max sequence))
	 (num-bins (if num-bins num-bins (floor (* (- top bottom) (expt (length sequence) 1/3)) (* 2 (iqr sequence)))))
	 (num-bins (+ 1 num-bins))
	 (boundaries (linspace bottom top :steps num-bins)))
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
    (linspace (+ bottom (/ gap-size 2)) top :step-size gap-size)))




;;; MCMC Metropolis-Hastings
(defun log-prior-flat (params data)
  (declare (ignore params data))
  0d0)

(defmacro prior-bounds-let ((&rest keys-low-high) &body body)
  (flet ((expand-keys (key-low-high)
	   (let* ((key (car key-low-high))
		  (param-name (read-from-string (symbol-name key))))
	     `(,param-name (the double-float (getf params ,key 0d0)))))
	 (expand-bounds (key-low-high)
	   (destructuring-bind (key low-expr high-expr) key-low-high
	     (let ((param-name (read-from-string (symbol-name key)))
		   (param-name-bound (read-from-string (concatenate 'string (symbol-name key) "-bound"))))
	       `(,param-name-bound (the double-float (if (< ,low-expr ,param-name ,high-expr)
					0d0
					(* -1d100 (- (exp (* (min (abs (- ,param-name ,high-expr)) (abs (- ,param-name ,low-expr))) 1d-5)) 1))))))))
	 (get-bound-name (key-low-high)
	   (destructuring-bind (key low-expr high-expr) key-low-high
	     (declare (ignore low-expr high-expr))
	     (let ((param-name-bound (read-from-string (concatenate 'string (symbol-name key) "-bound"))))
	       param-name-bound))))
    `(let* (,@(mapcar #'expand-keys keys-low-high)
	    ,@(mapcar #'expand-bounds keys-low-high)
	    (bounds-total (+ ,@(mapcar #'get-bound-name keys-low-high))))
       ,@body)))

;; log-liklihood can depend on x, y, params, error (additional distribution parameters)
;; in this case, it is just log normal (need x and y and mu (func x) and sigma)
;; in other cases, it could be log poisson (need lambda (func y) and k (y)) (no x, no sigma)
(defun log-normal (x mu sigma)
  (declare (optimize speed)
	   (double-float mu x)
	   (type (double-float 0d0 *) sigma))
  (+ (* -1/2 (log (* 2 pi))) (* -1 (log sigma)) (* -1/2 (expt (/ (- x mu) sigma) 2.0))))

(defun log-factorial (n)
  (reduce #'(lambda (x y) (+ x (log y))) (linspace 0 n :step-size 1)))

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
(defun create-log-liklihood-normal-weighted (fn params data error data-column-x data-column-y)
  (declare (optimize speed)
	   (cons data)
	   (function fn)
	   (simple-vector error))
  (let* ((x (elt data data-column-x))
	 (y (elt data data-column-y))
	 (error (if (= 1 (length error)) (make-array (length x) :initial-element (elt error 0)) error)))
    (declare (simple-vector x y error))
    (reduce #'+ (map 'vector #'(lambda (x y z) (log-normal y (apply fn x params) z)) x y error))))

(defun log-liklihood-normal (fn params data error)
  (create-log-liklihood-normal-weighted fn params data error 0 1))


;; (defun matrix-list-of-vectors (matrix)
;;   (cond ((numberp matrix) (error "Too few elements in matrix, looking at ~s." matrix))
;; 	((atom matrix) matrix)
;; 	((numberp (car matrix)) (make-array (length matrix) :initial-contents (mapcar #'(lambda (x) (float x 1.0d0)) matrix)))
;; 	(t (mapcar #'matrix-list-of-vectors matrix))))

(defun force-list (item)
  (if (or (consp item) (arrayp item))
      item
      (list item)))

(defun log-prior-fixer (log-prior params data)
  (let ((results (mapcar #'(lambda (lp d) (funcall lp params d)) log-prior data)))
    (mapcar #'(lambda (res fn) (if (numberp res) fn res)) results log-prior)))

(defun log-liklihood-fixer (log-liklihood fn params data error)
  (let ((results (mapcar #'(lambda (ll f d e) (funcall ll f params d e)) log-liklihood fn data error)))
    (mapcar #'(lambda (res fn) (if (numberp res) fn res)) results log-liklihood)))


(defun error-clean (item)
  (cond ((numberp item)
	 (list (make-array 1 :initial-element (float item 0d0)))) ; just a single number
	((vectorp item)
	 (list (map 'vector #'(lambda (x) (float x 0d0)) item))) ; just an array
	((and (consp item) (vectorp (car item)))
	 (mapcar #'(lambda (x) (map 'vector #'(lambda (y) (float y 0d0)) x)) item)) ; already a list of vectors
	((and (consp item) (not (consp (car item))))
	 (list (map 'vector #'(lambda (x) (float x 0d0)) item)))
	; it's pobably a list of lists, so just make the inner lists into vectors
	(t
	 (mapcar #'(lambda (x) (map 'vector #'(lambda (y) (float y 0d0)) x)) item))))
; (mapcar #'error-clean '(4 #(1 2 3) (#(1 2) #(4 5)) (1 2 3) ((1 2 3) (4 5 6))))

;; (defun data-proper-depth (item)
;;   (cond ((numberp (car item)) (list item))
;; 	((numberp (caar item)) (list item))
;; 	((consp (caar item)) item)
;; 	(t item)))

(defmacro n-apply (x &rest functions)
  (let ((gx (gensym)))
    `(let ((,gx ,x))
       (mapcar #'(lambda (f) (funcall f ,gx)) ',functions))))

(defun s-combinator (f g x)
  (funcall #'(lambda (a) (funcall f x a)) (funcall g x)))

(defun p-combinator (f h g x)
  (funcall #'(lambda (a) (funcall f (funcall h x) a)) (funcall g x)))

(defun k-combinator (x y)
  (declare (ignore y))
  x)

; list of data sets = list of list of vectors
(defun data-clean (item)
  (cond ((numberp (caar item)) ; number in a list in a list (not deep enough) (single dataset given)
	 (list (mapcar #'(lambda (y) (map 'vector #'(lambda (z) (float z 0d0)) y)) item)))
	((vectorp (caar item)) ; vector in a list in a list
	 (mapcar #'(lambda (x) (mapcar #'(lambda (y) (map 'vector #'(lambda (z) (float z 0d0)) y)) x)) item))
	((numberp (caaar item)) ; number in a list in a list in a list
	 (mapcar #'(lambda (x) (mapcar #'(lambda (y) (map 'vector #'(lambda (z) (float z 0d0)) y)) x)) item))	
	((or (consp (caaar item)) (vectorp (caaar item))) ; list/vector in a list in a list in a list (x's and y's are lists, not just numbers)
	 (mapcar #'(lambda (x) (mapcar #'(lambda (y) (mapcar #'(lambda (z) (map 'vector #'(lambda (a) (float a 0d0)) z)) y)) x)) item))))

;;; walker stuff
;; Create a walker will all associated fields. Coverts data and errors to double float lists of vectors for efficiency of computation (speed, consistency, easy of use (keeps sequence structure)
;; assumes fully list inputs

;; need to add list around data
;; need to remove matrix-list-of-vectors
;; need to make more general so these things don't need correcting
(defun create-walker (fn params data error log-liklihood log-prior &key (init t))
  "Create a walker that can be used in walker-* functions to do fitting and other
probabilistic analysis."
  (let* ((walker nil)
	 (fn (force-list fn))
	 (log-liklihood (force-list log-liklihood))
	 (log-prior (force-list log-prior))
	 (error (error-clean error))
;	 (error (if (= 1 (length fn)) (error-proper-depth error) error))
;	 (error (mapcar #'matrix-list-of-vectors error))
;n	 (data (if (= 1 (length fn)) (list (data-proper-depth data)) data))
	 ; (data (mapcar #'matrix-list-of-vectors data))
	 (data (data-clean data))
	 (log-prior (log-prior-fixer log-prior params data))
	 (log-liklihood (log-liklihood-fixer log-liklihood fn params data error))
	 (param-keys (get-plist-keys params))
	 (walker-context (make-plist '(:fn :param-keys :data :error :log-liklihood :log-prior) (list fn param-keys data error log-liklihood log-prior))))
    (if (not (= (length fn) (length log-liklihood) (length log-prior)))
	(error "fn, log-liklihood, and log-prior need to have the same number of elements, either 1 (scalars) or lists of n items. Currently ~s, ~s, and ~s." fn log-liklihood log-prior))
    (labels ((make-step (params)
	       (let* ((prob (+ (reduce #'+ (mapcar #'(lambda (ll f d e) (funcall ll f params d e)) log-liklihood fn data error))
			       (reduce #'+ (mapcar #'(lambda (lp d) (funcall lp params d)) log-prior data)))))
		 (list (list :prob prob :params params))))
	     (init-walker (&optional init-params)
	       (unless init-params (setf init-params params))
	       (setq walker (make-step init-params)))
	     (delete-walker ()
	       (setq fn nil
		     log-liklihood nil
		     log-prior nil
		     params nil
		     data nil
		     error nil
		     walker nil
		     param-keys nil
		     walker-context nil)))
      ;; initialize the walker
      (when init
	(init-walker))
      (lambda (msg &rest args)
	(case msg
	  (:get-context walker-context)
	  (:set-error (setq error (car args)))
	  (:get-walker walker)
	  (:get-last-walk (car walker))
	  (:burn-walks (nbutlast walker (car args)))
	  (:remove-recents (setf walker (subseq walker (car args))))
	  (:reset-walker (init-walker (car args)))
	  (:empty-walker (delete-walker))
	  (:make-step (make-step (car args)))
	  (:add-step (setq walker (let ((r (random 50000))) (if (and (= r 1) (> (length walker) 200000)) (setq walker (subseq walker 0 100000))) (nconc (car args) walker)))))))))

;;; Functions that act on the walker closure ;; primitives
(defun walker-get (the-walker &optional take)
  (subseq (funcall the-walker :get-walker) 0 (when take (floor take))))

(defun walker-set-error (the-walker error)
  (funcall the-walker :set-error error))

(defun walker-get-last-walk (the-walker)
  (funcall the-walker :get-last-walk))

(defun walker-burn-walks (the-walker burn-in)
  (null (funcall the-walker :burn-walks burn-in)))

(defun walker-remove-recents (the-walker num-remove)
  (null (funcall the-walker :remove-recents num-remove)))

(defun walker-reset (the-walker &optional init-params)
  (funcall the-walker :reset-walker init-params))

(defun walker-delete (the-walker)
  (funcall the-walker :empty-walker))

(defun walker-get-context (the-walker)
  (funcall the-walker :get-context))

(defun walker-make-step (the-walker params)
  (funcall the-walker :make-step params))

(defun walker-add-walks (the-walker walk-list)
  (funcall the-walker :add-step walk-list))

(defun walker-dry-step (the-walker &optional l-matrix (temperature 1))
  "The one that uses a full covariance matrix"
  (if (null l-matrix) (setq l-matrix (diagonal-covariance (get-plist-values (scale-plist 1e-2 (walker-get-median-params the-walker 1000))))))
  (let* ((previous-walk (walker-get-last-walk the-walker))
	 (prob0 (print (getf previous-walk :prob)))
	 (previous-params (print (getf previous-walk :params)))
	 (cov-values (get-covariant-sample (get-plist-values previous-params) l-matrix))
	 (next-params (print (make-plist (get-plist-keys previous-params) cov-values)))
	 (next-walk (car (walker-make-step the-walker next-params)))
	 (prob1 (print (getf next-walk :prob)))
	 (accepted-step (if (or (> prob1 prob0)
				(> (/ (- prob1 prob0) temperature) (log (random 1.0d0))))
			    (list next-walk)
			    (list previous-walk))))
    (print accepted-step)))

(defun walker-step (the-walker &optional l-matrix (temperature 1))
  "The one that uses a full covariance matrix"
  (if (null l-matrix) (setq l-matrix (diagonal-covariance (get-plist-values (scale-plist 1e-2 (walker-get-median-params the-walker 1000))))))
  (let* ((previous-walk (walker-get-last-walk the-walker))
	 (prob0 (getf previous-walk :prob))
	 (previous-params (getf previous-walk :params))
	 (cov-values (get-covariant-sample (get-plist-values previous-params) l-matrix))
	 (next-params (make-plist (get-plist-keys previous-params) cov-values))
	 (next-walk (car (walker-make-step the-walker next-params)))
	 (prob1 (getf next-walk :prob))
	 (accepted-step (if (or (> prob1 prob0)
				(> (/ (- prob1 prob0) temperature) (log (random 1.0d0))))
			    (list next-walk)
			    (list previous-walk))))
    (walker-add-walks the-walker accepted-step)))

(defun walker-many-steps (the-walker n &optional l-matrix)
  (if (null l-matrix) (print (car (setq l-matrix (diagonal-covariance (get-plist-values (scale-plist 1e-2 (walker-get-median-params the-walker))))))))
  (dotimes (i n)
    (walker-step the-walker l-matrix)))

;; TODO remove and add functionality to language
;; (print (car (diagonal-covariance (get-plist-values (scale-plist 1e-2 (walker-get-median-params woi 1000))))))
;; (print (standard-deviation (diff (remove-consecutive-duplicates (walker-get-param woi :w1-0-0)))))
;; (vgplot:plot (diff (remove-consecutive-duplicates (walker-get-param woi :w1-0-0))))


(defun walker-get-probs (the-walker &optional take)
  (let ((walker (walker-get the-walker take)))
    (mapcar #'(lambda (x) (getf x :prob)) walker)))

(defun walker-get-param (the-walker key &optional take)
  (let ((walker (walker-get the-walker take)))
    (mapcar #'(lambda (x) (getf (getf x :params) key)) walker)))

(defun walker-get-params (the-walker &optional take)
  (mapcar #'(lambda (x) (getf x :params)) (walker-get the-walker take)))

(defun walker-get-median-params (the-walker &optional take)
  (let* ((context (walker-get-context the-walker))
	 (param-keys (getf context :param-keys)))
    (make-plist
     param-keys
     (mapcar #'(lambda (x) (median (walker-get-param the-walker x take))) param-keys))))

(defun walker-get-most-likely-step (the-walker &optional take)
  (reduce #'(lambda (x y)
	      (if (> (getf x :prob) (getf y :prob)) x y))
	  (walker-get the-walker take)))

(defun walker-get-length (the-walker)
  (length (walker-get the-walker)))

(defun walker-remove-failed-walks (the-walker &optional take)
  (mapcon #'(lambda (x)
	      (if (equal (getf (car x) :prob) (getf (cadr x) :prob))
		  nil
		  (list (getf (car x) :params))))
	  (walker-get the-walker take)))

(defun walker-get-covariance-matrix (the-walker &optional take)
  (lplist-covariance (walker-remove-failed-walks the-walker take)))

;; Super useful for checking
(defun walker-check-for-complex-walks (the-walker &optional take)
  (let ((complex-walks-ps (mapcar #'(lambda (x) (some #'identity x)) (map-matrix #'complexp (lplist-to-l-matrix (diff-lplist (walker-remove-failed-walks the-walker take)))))))
    (if (some #'identity complex-walks-ps) complex-walks-ps nil)))

;; TODO add failsafe here if not enough valid steps to do cholesky
(defun walker-get-l-matrix (the-walker &optional take)
  (lplist-to-l-matrix (diff-lplist (walker-remove-failed-walks the-walker take))))

(defun walker-get-stddev-params (the-walker &optional take)
  (let* ((context (walker-get-context the-walker))
	 (param-keys (getf context :param-keys)))
    (if (< (walker-get-length the-walker) 10)
	(make-plist param-keys (make-list (length param-keys) :initial-element 0d0))
	(let ((l-matrix (walker-get-l-matrix the-walker take)))
	  (make-plist
	   param-keys
	   (mapcar #'(lambda (x)
		       (elt (elt l-matrix x) x))
		   (linspace 0 (- (length l-matrix) 1) :step-size 1)))))))

(defun walker-get-acceptance (the-walker &optional take)
  (let* ((probs (walker-get-probs the-walker take))
	 (clean-probs (remove-consecutive-duplicates probs)))
    (/ (length clean-probs) (length probs))))


(defun walker-adaptive-steps-full (the-walker &optional (n 100000) (temperature 1d6) (auto t) l-matrix)
  (let* ((param-keys (getf (walker-get-context the-walker) :param-keys))
	 (temp-steps 10000)
	 (temps (linspace temperature 1 :steps temp-steps))
	 (temp-temps (copy-seq temps))
	 (current-acceptance 0))
    (flet ((stable-probs-p (probs)
	     (< 5 (- (reduce #'max probs) (reduce #'min probs)) 10))
	   (get-optimal-mcmc-l-matrix (take)
	     (scale-matrix (/ (expt 2.38 2) (length param-keys))
			   (walker-get-l-matrix the-walker take))))
      (unless l-matrix
	(if (or (< (walker-get-length the-walker) 2000)
		(< (walker-get-acceptance the-walker 100) 0.1))
	    (setq l-matrix (diagonal-covariance (get-plist-values (scale-plist 1e-5 (walker-get-median-params the-walker)))))
	    (setq l-matrix (get-optimal-mcmc-l-matrix 2000))))
      (do ((i 0 (+ i 1)))
	  ((or (>= i n)
	       (and auto
		    (= 0 (mod i 1000))
		    (> i 40000)
		    (< 0.2 (walker-get-acceptance the-walker 1000) 0.5)
		    (stable-probs-p (walker-get-probs the-walker 2000)))))
	;; catch complex error and retry with differnt opimal l matrix sampling
	(walker-step the-walker l-matrix temperature)
	;; Annealing
	(when (and (= 0 (mod (floor i temp-steps) 2)) (< i (* temp-steps 95)))
	  (when (null temp-temps) (setf temp-temps (copy-seq temps)))
	  (setq temperature (pop temp-temps))
	  (when (null temp-temps) (walker-add-walks the-walker (list (walker-get-most-likely-step the-walker temp-steps)))))
	;; Regular l-matrix updating
	(when (and (= 0 (mod i 1000)) (> i 2000))
	  (setf current-acceptance (walker-get-acceptance the-walker 1000))
	  (if (< current-acceptance 0.2)
	      (setf l-matrix (scale-matrix 0.5 l-matrix))
	      (if (> current-acceptance 0.5)
		  (setf l-matrix (scale-matrix 2 l-matrix))
		  (setf l-matrix (get-optimal-mcmc-l-matrix 2000)))))))))

(defun walker-adaptive-steps (the-walker &optional (n 100000))
  (walker-adaptive-steps-full the-walker n 1d6 t nil))

(defun walker-get-data-column (the-walker column &optional (data-set 0))
  (elt (elt (getf (walker-get-context the-walker) :data) data-set) column))


(defun walker-catepillar-plots (the-walker &optional take)
  (let* ((param-keys (getf (walker-get-context the-walker) :param-keys))
	 (n (length param-keys)))
    (vgplot:close-all-plots)
    (vgplot:format-plot t "reset")
    (vgplot:format-plot t "set terminal pngcairo size 1920,1920 linewidth 1 font \"Arial,16\"")
    (vgplot:format-plot t "set output \"temp.png\"")
    (vgplot:format-plot t "set format y \"%.1e\"")
    (vgplot:format-plot t "set format x")
    (mapcar #'(lambda (x)
		(vgplot:subplot n 1 (position x param-keys))
		(vgplot:xlabel "" :replot nil)
		(vgplot:ylabel "" :replot nil)
		(vgplot:plot (walker-get-param the-walker x take) (format nil ";~a; w l" x)))
	    param-keys)
    (vgplot:xlabel "Counts")
    (vgplot:close-all-plots)
    (vgplot:new-plot)
    (vgplot:format-plot t "set terminal qt size 1920,1920 linewidth 1 font \"Arial,16\"")))

(defun walker-liklihood-plot (the-walker &optional take)
  (let ((probs (walker-get-probs the-walker take)))
    (vgplot:close-all-plots)
    (vgplot:format-plot t "set terminal qt size 1920,1080 linewidth 1 font \"Arial,16\"")
    (vgplot:xlabel "Step" :replot nil)
    (vgplot:ylabel "Log Liklihood" :replot nil)
    (vgplot:plot probs ";Liklihood;w l")))

(defun walker-fit-error-max-min (the-walker &key (take 1000) (x-column 0) (fn-number 0))
  (let* ((context (walker-get-context the-walker))
	 (fn (elt (getf context :fn) fn-number))
	 (data (elt (getf context :data) fn-number))
	 (x-data (elt data x-column))
	 (x-fit (linspace (elt x-data 0)
			  (elt x-data (- (length x-data) 1))
			  :steps 1000))
	 (previous-walks (walker-get the-walker take))
	 (all-ys (mapcar #'(lambda (x) (mapcar #'(lambda (walk) (apply fn x (getf walk :params))) previous-walks)) x-fit))
	 (max-ys (mapcar #'(lambda (y) (reduce #'max y)) all-ys))
	 (min-ys (mapcar #'(lambda (y) (reduce #'min y)) all-ys)))
    (list max-ys min-ys)))

(defun walker-plot-data-and-fit (the-walker &key (take 1000) (x-column 0) (y-column 1) (fn-number 0))
  (let* ((take (if (or (null take) (> take (walker-get-length the-walker))) (walker-get-length the-walker) take))
	 (context (walker-get-context the-walker))
	 (fn (elt (getf context :fn) fn-number))
	 (data (elt (getf context :data) fn-number))
	 (x-data (elt data x-column))
	 (y-data (elt data y-column))
	 (x-fit (linspace (elt x-data 0)
			  (elt x-data (- (length x-data) 1))
			  :steps 1000))
	 (params (walker-get-median-params the-walker take))
	 (y-fit (map 'vector #'(lambda (x) (apply fn x params)) x-fit))
	 (previous-walks (walker-get the-walker take))
	 (sorted-walks (subseq (sort previous-walks #'> :key #'(lambda (x) (getf x :prob))) 0 (floor (* 0.66 take))))
	 (all-ys (mapcar #'(lambda (x) (mapcar #'(lambda (walk) (apply fn x (getf walk :params))) sorted-walks)) x-fit))
	 (max-ys (mapcar #'(lambda (y) (reduce #'max y)) all-ys))
	 (min-ys (mapcar #'(lambda (y) (reduce #'min y)) all-ys)))
    (vgplot:close-all-plots)
    (vgplot:format-plot t "set terminal qt size 1920,1080 linewidth 3 font \"Arial,\"")
    (vgplot:format-plot t "set format x \"%.4g\"")
    (vgplot:format-plot t "set format y \"%.4g\"")
    (vgplot:xlabel "x-data" :replot nil)
    (vgplot:ylabel "y-data" :replot nil)
    (vgplot:plot x-fit max-ys ";fit stddev upper limit; w l lc rgb \"green\""
	  x-fit min-ys ";fit stddev lower limit; w l lc rgb \"green\""
	  x-fit y-fit ";fit; w l lc rgb \"red\"" 
	  x-data y-data ";data; w p pt 6 ps 2 lc rgb \"black\"")
    params))

;; TODO Validate
;; (defun walker-plot-data-and-fit-only (the-walker &optional take (x-column 0) (y-column 1) (fn-number 0))
;;   (let* ((take (when take (if (> take (walker-get-length the-walker)) (walker-get-length the-walker) take)))
;; 	 (context (walker-get-context the-walker))
;; 	 (fn (getf context :fn))
;; 	 (fn (if (consp fn) (elt fn fn-number) fn))
;; 	 (data (getf context :data))
;; 	 (x-data (elt data x-column))
;; 	 (y-data (elt data y-column))
;; 	 (x-fit (linspace (- (elt x-data 0) 50)
;; 			  (+ (elt x-data (- (length x-data) 1)) 50)
;; 			  :steps 1000))
;; 	 (params (walker-get-median-params the-walker take))
;; 	 (y-fit (map 'vector #'(lambda (x) (nth 0 (apply fn (list x) params))) x-fit)))
;;     (vgplot:close-all-plots)
;;     (vgplot:format-plot t "set terminal qt size 1920,1080 linewidth 3 font \"Arial,16\"")
;;     (vgplot:format-plot t "set format x \"%.4g\"")
;;     (vgplot:format-plot t "set format y \"%.4g\"")
;;     (vgplot:xlabel "x-data" :replot nil)
;;     (vgplot:ylabel "y-data" :replot nil)
;;     (vgplot:plot x-fit y-fit ";fit; w l lc rgb \"red\""
;; 	  x-data y-data ";data; w p pt 6 ps 2 lc rgb \"black\"")
;;     params))

(defun walker-plot-residuals (the-walker &key (take 1000) (x-column 0) (y-column 1) (fn-number 0))
  (let* ((take (if (or (null take) (> take (walker-get-length the-walker))) (walker-get-length the-walker) take))
	 (context (walker-get-context the-walker))
	 (fn (elt (getf context :fn) fn-number))
	 (data (elt (getf context :data) fn-number))
	 (error (elt (getf context :error) fn-number))
	 (x-data (elt data x-column))
	 (y-data (elt data y-column))
	 (error (if (= 1 (length error)) (make-array (length x-data) :initial-element (elt error 0)) error))
	 (params (walker-get-median-params the-walker take))
	 (y-fit (map 'vector #'(lambda (x) (apply fn x params)) x-data))
	 (y-residuals (map 'vector #'(lambda (yf y) (- yf y)) y-fit y-data)))
    (vgplot:close-all-plots)
    (vgplot:format-plot t "set terminal qt size 1920,1080 linewidth 3 font \"Arial,16\"")
    (vgplot:format-plot t "set format x \"%.4g\"")
    (vgplot:format-plot t "set format y \"%.4g\"")
    (vgplot:xlabel "x-data" :replot nil)
    (vgplot:ylabel "y-data" :replot nil)
    (vgplot:plot x-data y-residuals ";residuals; w p pt 6 ps 2 lc rgb \"black\""
		 x-data error ";point error; w p pt 2 ps 1 lc rgb \"red\""
		 x-data (make-array (length x-data) :initial-element 0) ";baseline; w l lc rgb \"red\"")))

(defun walker-plot-data (the-walker &key (x-column 0) (y-column 1) (fn-number 0))
  (let* ((context (walker-get-context the-walker))
	 (data (elt (getf context :data) fn-number))
	 (x-data (elt data x-column))
	 (y-data (elt data y-column)))
    (vgplot:close-all-plots)
    (vgplot:format-plot t "set terminal qt size 1920,1080 linewidth 3 font \"Arial,16\"")
    (vgplot:xlabel "x-data" :replot nil)
    (vgplot:ylabel "y-data" :replot nil)
    (vgplot:plot x-data y-data ";data; w lp pt 6 ps 2 lc rgb \"black\"")))

(defun permute-params (params)
  (mapcon
   #'(lambda (y)
       (butlast (maplist
		 #'(lambda (x)
		     (list (car x) (car (last x))))
		 (reverse y))))
   params))

(defun walker-2d-plot (the-walker keys &optional take)
  (let* ((x-key (first keys))
	 (y-key (second keys))
	 (x (walker-get-param the-walker x-key take))
	 (y (walker-get-param the-walker y-key take)))
    (vgplot:close-all-plots)
    (vgplot:format-plot t "reset")
    (vgplot:format-plot t "set terminal qt size 1920,1080 linewidth 1 font \"Arial,16\"")
    (vgplot:xlabel (format nil "~a" x-key) :replot nil)
    (vgplot:ylabel (format nil "~a" y-key) :replot nil)
    (vgplot:plot x y ";My plot; w p pt 3 ps 4")))

;; TODO only extract each parameter once
(defun walker-all-2d-plots (the-walker &optional take)
  (let* ((param-keys (getf (walker-get-context the-walker) :param-keys))
	 (key-pairs (permute-params param-keys))
	 (len (- (length param-keys) 1))
	 (walker (walker-remove-failed-walks the-walker take)))
    (labels ((get-one-param (key)
	       (mapcar #'(lambda (x) (getf x key)) walker))
	     (get-key-pair-index (x-key y-key)
	       (let ((x-pos (position x-key param-keys))
		     (y-pos (position y-key param-keys)))
		 (+ (* len (- y-pos 1)) x-pos)))
	     (inner-plotter (key-pair)
	       (let* ((x-key (second key-pair))
		      (y-key (first key-pair))
		      (x (get-one-param x-key))
		      (y (get-one-param y-key)))
		 (vgplot:subplot len len (get-key-pair-index x-key y-key))
		 (vgplot:xlabel (format nil "~a" x-key) :replot nil)
		 (vgplot:ylabel (format nil "~a" y-key) :replot nil)
		 (vgplot:plot x y ";; w p pt 3 ps 4"))
	       ))
      (vgplot:close-all-plots)
      (vgplot:format-plot t "reset")
      (vgplot:format-plot t "set terminal pngcairo size 3840,3840 linewidth 1 font \"Arial,16\"")
      (vgplot:format-plot t "set output \"temp.png\"")
      (vgplot:format-plot t "set format y \"%.4g\"")
      (vgplot:format-plot t "set format x \"%.4g\"")
      (mapcar #'(lambda (x) (inner-plotter x)) key-pairs)
      (vgplot:close-all-plots))))

(defun walker-param-histo (the-walker key &optional take)
  (let* ((bins 20)
	 (param-values (sort (walker-get-param the-walker key take) #'<))
	 (histo-x (make-histo-x param-values bins))
	 (histo (make-histo param-values bins)))
    (vgplot:close-all-plots)
    (vgplot:format-plot t "reset")
    (vgplot:format-plot t "set terminal qt size 1920,1080 linewidth 1 font \"Arial,16\"")
    (vgplot:xlabel (format nil "~a" key) :replot nil)
    (vgplot:ylabel "Counts" :replot nil)
    (vgplot:plot histo-x histo ";Histogram; w l")))

(defun walker-construct-print-list (the-walker &optional take)
  (let ((context (walker-get-context the-walker)))
    `(:param-keys ,(getf context :param-keys)
      :data ,(getf context :data)
      :error ,(getf context :error)
      :walker ,(walker-get the-walker take))))

(defun walker-save (the-walker filename &optional take)
  (with-open-file (out filename
                       :direction :output
                       :if-exists :supersede)
    (with-standard-io-syntax
      (null (print (walker-construct-print-list the-walker take) out)))))

(defun walker-load (filename fn log-liklihood log-prior)
  (with-open-file (in filename
		      :direction :input)
    (with-standard-io-syntax
      (let* ((full (read in))
	     (data (getf full :data))
	     (error (getf full :error))
	     (walks (getf full :walker))
	     (dummy-params (getf (elt walks 0) :params))
	     (the-walker (create-walker fn dummy-params data error log-liklihood log-prior :init nil)))
	(walker-add-walks the-walker walks)
	the-walker))))

(defun walker-easy-commands (the-walker msg &key (take 1000) cat-probs-take (n 100000) (temp 1000) filename reset-params)
  (case msg
    (:plot (walker-plot-data-and-fit the-walker :take take))
    (:data (walker-plot-data the-walker))
    (:cat (walker-catepillar-plots the-walker cat-probs-take))
    (:accept (walker-get-acceptance the-walker take))
    (:reset (walker-reset the-walker reset-params))
    (:rm (walker-delete the-walker))
    (:astep (walker-adaptive-steps-full the-walker n temp t))
    (:save (walker-save the-walker filename take))
    (:2d (walker-all-2d-plots the-walker take))
    (:probs (walker-liklihood-plot the-walker cat-probs-take))
    (:medians (walker-get-median-params the-walker take))
    (:stddev (walker-get-stddev-params the-walker take))
    (:l-matrix (walker-get-l-matrix the-walker take))
    (:burn (walker-burn-walks the-walker take))
    (otherwise
     (format t "Accepted Commands:~%")
     (format t "~{~s~^, ~}" '(:plot :data :cat :accept :reset :rm :astep :save :2d :probs :medians :stddev :l-matrix :burn)))))



;;; Functionality for a set of similar walker instances
(defun walker-set-save (the-walker-set filename &optional take)
  (let ((walker-print-list (mapcar #'(lambda (x) (walker-construct-print-list x take)) the-walker-set)))
    (with-open-file (out filename
			 :direction :output
			 :if-exists :supersede)
      (with-standard-io-syntax
	(null (print walker-print-list out))))))

(defun walker-set-load (filename fn log-liklihood log-prior)
  (with-open-file (in filename
		      :direction :input)
    (with-standard-io-syntax
      (let* ((full (read in)))
	(mapcar #'(lambda (x)
		    (let* ((param-keys (getf x :param-keys))
			   (dummy-params (make-plist param-keys (make-list (length param-keys) :initial-element 0.0d0)))
			   (data (getf x :data))
			   (error (getf x :error))
			   (walks (getf x :walker))
			   (the-walker (create-walker fn dummy-params data error log-liklihood log-prior :init nil)))
		      (walker-add-walks the-walker walks)
		      the-walker))
		full)))))

(defun walker-set-delete (the-walker-set)
  (mapcar #'walker-delete the-walker-set))

(defun walker-set-get-median-params (the-walker-set key &optional take)
  (mapcar #'(lambda (x) (getf (walker-get-median-params x take) key)) the-walker-set))

(defun walker-set-get-stddev-params (the-walker-set key &optional take)
  (mapcar #'(lambda (x) (getf (walker-get-stddev-params x take) key)) the-walker-set))

(defun walker-set-plot-param (the-walker-set key &optional take)
  (vgplot:plot (walker-set-get-median-params the-walker-set key take)))

;; TODO is eval the only way???
(defmacro walker-get-f (the-walker exp &optional take)
  `(let ((median-params (walker-get-median-params ,the-walker ,take)))
     (labels ((replace-param-names (item)
		(cond ((null item) nil)
		      ((atom item) (if (char= #\: (char (format nil "~s" item) 0)) (getf median-params item) item))
		      ((consp item)
		       (cons (replace-param-names (car item))
			     (replace-param-names (cdr item)))))))
       (let ((filled-exp (replace-param-names ',exp)))
	  (eval filled-exp)))))

;; TODO eval again, boooo.
(defmacro walker-set-get-f (the-walker-set exp &optional take)
  `(mapcar #'(lambda (x) (eval (walker-get-f x ',exp ,take))) ,the-walker-set))

;; TODO make it actually work
(defun walker-set-plot-func (the-walker-set func)
  (vgplot:plot (mapcar #'(lambda (x) (apply func x)) the-walker-set)))

(defun walker-set-easy-commands (the-walker-set index msg
				 &key (take 1000) (n 100000) (temp 1d6) reset-params cat-probs-take plot-param filename)
  (if plot-param
      (walker-set-plot-param the-walker-set plot-param take)
      (if (eq msg :save)
	  (walker-set-save the-walker-set filename take)
       (let ((the-walker (elt the-walker-set index)))
	 (walker-easy-commands the-walker msg :take take :take take :cat-probs-take cat-probs-take :n n :temp temp :reset-params reset-params :filename filename)))))

;;; End of walker functionality



;;; Data extraction
(defun 3d-plot-file (filename &key (map t))
  (vgplot:format-plot t "reset")
  (vgplot:format-plot t "set terminal qt size 1920,1080 linewidth 1 font \"Arial,16\"")
  (vgplot:format-plot t "set ticslevel 0")
  (vgplot:format-plot t "set pm3d depthorder")
  (vgplot:format-plot t "set cblabel \"Field Offset (Oe)\"")
  (if map
      (vgplot:format-plot t "set view map")
      (vgplot:format-plot t "unset view"))
  (vgplot:format-plot t "unset key")
  (vgplot:format-plot t "unset grid")
  (vgplot:format-plot t "set pm3d at sb")
  (vgplot:xlabel "X Pos" :replot nil)
  (vgplot:ylabel "Y Pos" :replot nil)
  (vgplot:zlabel "Field Offset (Oe)" :replot nil)
  (vgplot:format-plot t (format nil "splot \"~a\" u 1:2:3 w pm3d" filename)))

(defun 3d-plot-admr-file (filename &key (map t))
  (vgplot:format-plot t "reset")
  (vgplot:format-plot t "set terminal qt size 1920,1080 linewidth 1 font \"Arial,40\"")
  (vgplot:format-plot t "set ticslevel 0")
  (vgplot:format-plot t "set pm3d depthorder")
  (vgplot:format-plot t "set cblabel \"Absorption (dB)\"")
  (if map
      (vgplot:format-plot t "set view map")
      (vgplot:format-plot t "unset view"))
  (vgplot:format-plot t "unset key")
  (vgplot:format-plot t "unset grid")
  (vgplot:format-plot t "set pm3d at sb")
  (vgplot:xlabel "X Field (Oe)" :replot nil)
  (vgplot:ylabel "Y Field (Oe)" :replot nil)
  (vgplot:zlabel "S21 (a.u.)" :replot nil)
  (vgplot:format-plot t (format nil "splot \"~a\" u 1:2:7 w pm3d, \"\" u (-1*$1):(-1*$2):8 w pm3d" filename)))

(defun file->file-specs (filename &key (delim #\tab))
  "Tells the user various metrics about the file. Lines, header lines, data lines, data rows, data sets."
  (with-open-file (in filename :direction :input)
    (labels ((get-lines (&optional (num-lines 0) (found-data nil) (data-length nil) (data-rows nil))
	       (let ((line (string-right-trim '(#\return) (read-line in nil nil)))) ; allow for Windows files
		 (cond ((string= line "NIL") ; end of file - return info ; ordered by freque
			(list :file-lines num-lines :header-lines found-data :data-length data-length :data-rows (unless data-rows (- num-lines found-data))data-rows :num-pages (if data-rows (floor (- num-lines found-data) data-rows) 1)))
		       ((and (string= line "") found-data (not data-rows)) ; set data rows
			(get-lines num-lines found-data data-length (- num-lines found-data)))
		       ((string= line "") ; ignore line
			(get-lines num-lines found-data data-length data-rows))
		       ((and (numberp (read-from-string (elt (split-string #\tab line) 0))) (not found-data)) ; first line w nums?
			(get-lines (+ num-lines 1) num-lines (length (split-string delim line)) data-rows))
		       (t ; regular line - just increase lines
			(get-lines (+ num-lines 1) found-data data-length data-rows))))))
      (get-lines))))

;;; Designing a data toolbox
;; Fundamentally, just stream all the numbers into one long array
;; then converts it to a 2d array, as that's what is usually dealt with
;; then the user could make it 1d or 3d or whatever they like after

;;; File stuff
(defun file->data-list (filename &key (file-specs nil) (delim #\tab) (columns t))
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
		 (push (mapcar #'read-from-string (split-string delim line)) file-contents)
		 (read-file stream))))
      (with-open-file (in filename :direction :input)
	(remove-header-lines header-lines in)
	(read-file in)
	(if columns
	    (transpose-matrix (reverse file-contents))
	    (reverse file-contents))))))

(defun file->walker (filename fn params error log-liklihood log-prior)
  "Generic function to read file and create walker based on specifying all inputs"
  (let* ((data (file->data-list filename)))
    (create-walker fn params data error log-liklihood log-prior)))


(defun file->plot (filename &key (x-column 0) (y-column 1))
  (let* ((data (file->data-list filename)))
    (vgplot:plot (elt data x-column) (elt data y-column))))


;; file searching
(defun walk-dirs (dir)
  (if (uiop:directory-exists-p dir) ;; it's a dir?
      (let ((items (nconc (uiop:subdirectories dir) (uiop:directory-files dir))))
	(mapcar #'(lambda (d) (walk-dirs d)) items))
      dir))

(defun get-filename (dir &key include exclude)
  "Returns all files under dir (and its subdirs) whose name (including directory name) matches ALL specified patterns"
  (let* ((all-files (alexandria:flatten (walk-dirs dir)))
	 (files (remove-if-not #'(lambda (y)
				   (and (every #'(lambda (g) (search g (namestring y))) include)
					(notany #'(lambda (g) (search g (namestring y))) exclude)))
			       all-files)))
    (if (= (length files) 1)
	(car files)
	files)))


(defun show ()
  "Shows the plots generated from 2d histogram and catepillar functions"
  (uiop:run-program "feh ./temp.png -3 5"))

;; (defparameter gnuplot-stream (ltk:do-execute "gnuplot" nil))
;; The difficult thing he takes care of is saving the data to a temporary file for plotting
;; Just need to make it easy to do this in a customizable way
;; Then you can just build plots from these files
;; a 'plot' commands builds all its input, but isn't sent to the buffer until to file is
;; 'execute'd.
;; All other things (format changes, labels, etc) are written at they are recieved
;; multiplot needs to be way better

;; Could write all my own commands to gnuplot this way
;; His wrapper is soo thin as to be nearly useless I think
;; Would be more useful to be able to add plot with additionl plot commands
;; then just make a catalogue of the most useful and helpful commands available

