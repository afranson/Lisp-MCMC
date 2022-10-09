;;; Log liklihood and log-normal

(in-package :mcmc-fitting)

(defun log-prior-flat (params data)
  (declare (ignore params data))
  0d0)

(defmacro prior-bounds-let ((&rest keys-low-high) &body body)
  (declare (sb-ext:muffle-conditions sb-ext:compiler-note))
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
					(* -1d10 (- (exp (* (min (abs (- ,param-name ,high-expr)) (abs (- ,param-name ,low-expr))) 1d-5)) 1))))))))
	 (get-bound-name (key-low-high)
	   (destructuring-bind (key low-expr high-expr) key-low-high
	     (declare (ignore low-expr high-expr))
	     (let ((param-name-bound (read-from-string (concatenate 'string (symbol-name key) "-bound"))))
	       param-name-bound))))
    `(let* (,@(mapcar #'expand-keys keys-low-high)
	    ,@(mapcar #'expand-bounds keys-low-high)
	    (bounds-total (+ ,@(mapcar #'get-bound-name keys-low-high))))
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
  (reduce #'(lambda (x y) (+ x (log y))) (up-to n)))

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
	   (cons data)
	   (function fn)
	   (simple-vector stddev))
  (let* ((x (elt data data-column-x))
	 (y (elt data data-column-y))
	 (stddev (if (= 1 (length stddev)) (make-array (length x) :initial-element (elt stddev 0)) stddev)))
    (declare (simple-vector x y stddev))
    (reduce #'+ (map 'vector #'(lambda (x y z) (log-normal y (apply fn x params) z)) x y stddev))))

(defun log-liklihood-normal (fn params data stddev)
  (create-log-liklihood-normal-weighted fn params data stddev 0 1))



