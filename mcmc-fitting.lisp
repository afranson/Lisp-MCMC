;;; mcmc-fitting.lisp
#|
create walker
advance it
visualize it
|#

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



;;; data manipulation
(defun diff-list (x-data y-data)
  "Returns two lists, the first is the center points of x-data, the second is the slopes."
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



;; Get the generalized sqrt of the variance (get the standard deviation of a multivariate covariance)
;; 5x faster than gsl implementation
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
	    (setf (aref l-array i k) (sqrt (- (aref cov-array i k) tmp-sum)))
	    (setf (aref l-array i k) (/ (- (aref cov-array i k) tmp-sum) (aref l-array k k))))))
    (array-matrix l-array)))


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
    (mapcar #'(lambda (el) (mapcar #'(lambda (row) (elt row el)) matrix)) (up-to num-columns))))

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


;;; MCMC Metropolis-Hastings
(defun force-list (item)
  (if (or (consp item))
      item
      (list item)))

(defun clean-stddev (stddev fn)
  "Assumes user competency. If you give one number, it uniformly applies it to all fns. If you give a list/vector, it applies that to each fn. If you supply a list of lists or list of vectors, it assumes you know what you're doing."
  (cond ((numberp stddev)
	 (repeat (length fn) (make-array 1 :initial-element stddev)))
	((and (consp stddev) (numberp (elt stddev 0)))
	 (list (repeat (length fn) (make-array (length stddev) :initial-contents stddev))))
	((vectorp stddev)
	 (repeat (length fn) stddev))
	((and (consp stddev) (consp (elt stddev 0)))
	 (mapcar #'(lambda (x) (make-array (length x) :initial-contents x)) stddev))))

(defun clean-data (data fn)
  "Fixes most common error of supply just a (list x y) as data and makes the corresponding x and y into vectors."
  (labels ((rec (data)
	     (cond ((and (consp data) (consp (car data)))
		    (mapcar #'(lambda (x) (rec x)) data))
		   ((and (consp data) (vectorp (car data)))
		    data)
		   (t
		    (apply #'vector data)))))
    (rec (if (= 1 (length fn)) (list data) data))))

(defun create-walker-data (data &rest columns)
  "Takes a larger dataset (contains x, y, stddev, freq, temp, etc.) and extracts just the columns you want into a walker-friendly format."
  (mapcar #'(lambda (x) (apply #'vector x))
	  (mapcar #'(lambda (y) (elt data y)) columns)))

(defun to-double-floats (l)
  "Converts all numbers and vectors in tree to double-float type"
  (cond ((numberp l)
	 (coerce l 'double-float))
	((symbolp l) ; for p-lists
	 l)
	((consp l)
	 (mapcar #'to-double-floats l))
	((vectorp l)
	 (map 'vector #'(lambda (x) (coerce x 'double-float)) l))))

(defun log-prior-fixer (log-prior params data)
  "Checks if your prior returns a function (i.e. changes shape based on the data provided to it). If so, gets that function. If not, returns the prior you gave it."
  (let ((results (mapcar #'(lambda (lp d) (funcall lp params d)) log-prior data)))
    (mapcar #'(lambda (res fn) (if (numberp res) fn res)) results log-prior)))

(defun log-liklihood-fixer (log-liklihood fn params data error)
  "Checks if your liklihood returns a function (i.e. changes shape based on the data provided to it). If so, gets that function. If not, returns the liklihood you gave it."
  (let ((results (mapcar #'(lambda (ll f d e) (funcall ll f params d e)) log-liklihood fn data error)))
    (mapcar #'(lambda (res fn) (if (numberp res) fn res)) results log-liklihood)))


(defmacro lambda-get-set-etc ((&rest vars) &rest message-args-commands)
  "Creates :set-var and :get-var functions for a let-over-lambda construct and allows for further messages as well.

Used as:
(lamda-get-set-etc
   (var1 var2 var3) ; these all receive get-var and set-var methods
   (:key-value0 (args0) (function0 ...))
   (:key-value1 (args1) (function1 ...))
   ...)"
  `(lambda (msg &rest args)
     (case msg
       ,@(mapcar #'(lambda (x) `(,(symb-keyword "GET-" x) ,x)) vars)
       ,@(mapcar #'(lambda (x) `(,(symb-keyword "SET-" x) (apply #'(lambda (input) (setf ,x input)) args))) vars)
       ,@(mapcar #'(lambda (x) `(,(car x) (apply #'(lambda ,@(cdr x)) args))) message-args-commands)
       (t (error "Bad msg '~s' with args: ~s" msg args)))))

;;; walker stuff
;; Create a walker will all associated fields. Coverts data and errors to double float lists of vectors for efficiency of computation (speed, consistency, easy of use (keeps sequence structure)
;; assumes fully list inputs

;; TODO name value cells for easier inspection of the closure
;; need to add list around data
;; need to remove matrix-list-of-vectors
;; need to make more general so these things don't need correcting
(defun walker-init (&key fn data params (stddev 1) (log-liklihood #'log-liklihood-normal) (log-prior #'log-prior-flat) (init t) trusted &allow-other-keys)
  "Create a walker that can be used in walker-* functions to do fitting and other
probabilistic analysis. If trusted, doesn't do data dn stddev checking."
  (let* ((walker nil)
	 (fn (force-list fn))
	 (stddev (if trusted stddev (to-double-floats (clean-stddev stddev fn))))
	 (data (if trusted data (to-double-floats (clean-data data fn))))
	 (params (to-double-floats params))
	 (log-liklihood (force-list log-liklihood))
	 (log-prior (force-list log-prior))
	 (log-prior (log-prior-fixer log-prior params data))
	 (log-liklihood (log-liklihood-fixer log-liklihood fn params data stddev))
	 (param-keys (get-plist-keys params)))
    (assert (= (length fn) (length log-liklihood) (length log-prior) (length data))
	    nil
	    "fn, data, log-liklihood, and log-prior need to have the same number of elements, either 1 (scalars) or lists of n items. Currently ~s, ~s, ~s, and ~s." fn data log-liklihood log-prior)
    (labels ((make-step (params)
	       (let* ((prob (+ (reduce #'+ (mapcar #'(lambda (ll f d e) (funcall ll f params d e)) log-liklihood fn data stddev))
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
		     stddev nil
		     walker nil
		     param-keys nil)))
      ;; initialize the walker
      (when init
	(init-walker))
      (lambda-get-set-etc
	  (walker fn data params stddev log-liklihood log-prior)
	  (:get-last-walk () (car walker))
	  (:burn-walks (n) (nbutlast walker n))
	  (:remove-recents (n) (setf walker (subseq walker n)))
	  (:reset-walker (init-params) (init-walker init-params))
	  (:empty-walker () (delete-walker))
	  (:make-step (params) (make-step params))
	  (:add-step (step) (setq walker (nconc step walker)))))))

(defmacro walker-create (name &rest rest &key fn data params (stddev 1) (log-liklihood #'log-liklihood-normal) (log-prior #'log-prior-flat) (init t) (docstring ""))
  "Generates a walker with 'name' that can be called as (name msg arg ...) to perform the 'msg' command or used as a variable name in functions."
  (declare (ignorable fn data params stddev log-liklihood log-prior init))
  `(progn
     (defparameter ,name (apply #'walker-init (list ,@rest)) ,docstring)
     (setf (symbol-function ',name) ,name)))

;; Calling the inner methods of the closure are as simple as (name :msg args)


;;; Functions that act on the walker closure ;; primitives for namespace convenience
(defun walker-get (the-walker &optional take)
  (subseq (funcall the-walker :get-walker) 0 (when take (floor take))))

(defun walker-get-param-keys (the-walker)
  (get-plist-keys (funcall the-walker :get-params)))

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

(defun walker-make-step (the-walker params)
  (funcall the-walker :make-step params))

(defun walker-add-walks (the-walker walk-list)
  (funcall the-walker :add-step walk-list))

(defun walker-step (the-walker &key l-matrix (temperature 1))
  "The one that uses a full covariance matrix"
  (if (null l-matrix) (setq l-matrix (diagonal-covariance (get-plist-values (scale-plist 1e-2 (walker-get-median-params the-walker 1000))))))
  (let* ((previous-walk (funcall the-walker :get-last-walk))
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
    (walker-step the-walker :l-matrix l-matrix)))

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
  (let* ((param-keys (walker-get-param-keys the-walker)))
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

;; Useful for inspecting when algorithm fails
(defun walker-check-for-complex-walks (the-walker &optional take)
  (let ((complex-walks-ps (mapcar #'(lambda (x) (some #'identity x)) (map-matrix #'complexp (lplist-to-l-matrix (diff-lplist (walker-remove-failed-walks the-walker take)))))))
    (if (some #'identity complex-walks-ps) complex-walks-ps nil)))

;; TODO add failsafe here if not enough valid steps to do cholesky
(defun walker-get-l-matrix (the-walker &optional take)
  (lplist-to-l-matrix (diff-lplist (walker-remove-failed-walks the-walker take))))

(defun walker-get-stddev-params (the-walker &optional take)
  (let* ((param-keys (walker-get-param-keys the-walker)))
    (if (< (walker-get-length the-walker) 10)
	(make-plist param-keys (make-list (length param-keys) :initial-element 0d0))
	(let ((l-matrix (walker-get-l-matrix the-walker take)))
	  (make-plist
	   param-keys
	   (mapcar #'(lambda (x)
		       (elt (elt l-matrix x) x))
		   (up-to (1- (length l-matrix)))))))))

(defun walker-get-acceptance (the-walker &optional take)
  (let* ((probs (walker-get-probs the-walker take))
	 (clean-probs (remove-consecutive-duplicates probs)))
    (/ (length clean-probs) (length probs))))


(defun walker-adaptive-steps-full (the-walker &optional (n 100000) (temperature 1d3) (auto t) l-matrix)
  (let* ((param-keys (walker-get-param-keys the-walker))
	 (temp-steps 10000)
	 (temps (linspace temperature 1 :num temp-steps))
	 (temp-temps (copy-seq temps))
	 (current-acceptance 0))
    (flet ((stable-probs-p (probs)
	     (< 5 (- (reduce #'max probs) (reduce #'min probs)) 15))
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
	;; catch complex stddev and retry with differnt opimal l matrix sampling
	(walker-step the-walker :l-matrix l-matrix :temperature temperature)
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
  (walker-adaptive-steps-full the-walker n 1d3 t nil))

(defun walker-catepillar-plots (the-walker &optional take)
  (let* ((param-keys (walker-get-param-keys the-walker))
	 (n (length param-keys)))
    (vgplot:close-all-plots)
    (vgplot:format-plot t "reset")
    (vgplot:format-plot t "set terminal pngcairo size 1920,1920 linewidth 1 font \"Arial,18\"")
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
    (vgplot:format-plot t "set terminal qt size 1920,1920 linewidth 1 font \"Arial,18\"")))

(defun walker-liklihood-plot (the-walker &optional take)
  (let ((probs (walker-get-probs the-walker take)))
    (vgplot:close-all-plots)
    (vgplot:format-plot t "set terminal qt size 1920,1080 linewidth 1 font \"Arial,18\"")
    (vgplot:xlabel "Step" :replot nil)
    (vgplot:ylabel "Log Liklihood" :replot nil)
    (vgplot:plot probs ";Liklihood;w l")))

(defun walker-fit-stddev-max-min (the-walker &key (take 1000) (x-column 0) (fn-number 0))
  (let* ((fn (elt (funcall the-walker :get-fn) fn-number))
	 (data (elt (funcall the-walker :get-data) fn-number))
	 (x-data (elt data x-column))
	 (x-fit (linspace (elt x-data 0)
			  (elt x-data (- (length x-data) 1))
			  :num 1000))
	 (previous-walks (walker-get the-walker take))
	 (all-ys (mapcar #'(lambda (x) (mapcar #'(lambda (walk) (apply fn x (getf walk :params))) previous-walks)) x-fit))
	 (max-ys (mapcar #'(lambda (y) (reduce #'max y)) all-ys))
	 (min-ys (mapcar #'(lambda (y) (reduce #'min y)) all-ys)))
    (list max-ys min-ys)))

(defun walker-plot-data-and-fit (the-walker &key (take 1000) (x-column 0) (y-column 1) (fn-number 0))
  (let* ((take (if (or (null take) (null (nthcdr take (funcall the-walker :get-walker))))
		   (walker-get-length the-walker)
		   take))
	 (fn (elt (funcall the-walker :get-fn) fn-number))
	 (data (elt (funcall the-walker :get-data) fn-number))
	 (x-data (elt data x-column))
	 (y-data (elt data y-column))
	 (x-fit (linspace (elt x-data 0)
			  (elt x-data (- (length x-data) 1))
			  :num 1000))
	 (params (walker-get-median-params the-walker take))
	 (y-fit (map 'vector #'(lambda (x) (apply fn x params)) x-fit))
	 (previous-walks (walker-get the-walker take))
	 (sorted-walks (subseq (sort previous-walks #'> :key #'(lambda (x) (getf x :prob))) 0 (floor (* 0.66 take))))
	 (all-ys (mapcar #'(lambda (x) (mapcar #'(lambda (walk) (apply fn x (getf walk :params))) sorted-walks)) x-fit))
	 (max-ys (mapcar #'(lambda (y) (reduce #'max y)) all-ys))
	 (min-ys (mapcar #'(lambda (y) (reduce #'min y)) all-ys)))
    (vgplot:close-all-plots)
    (vgplot:format-plot t "set terminal qt size 1920,1080 linewidth 3 font \"Arial,18\"")
    (vgplot:format-plot t "set format x \"%.4g\"")
    (vgplot:format-plot t "set format y \"%.4g\"")
    (vgplot:xlabel "x-data" :replot nil)
    (vgplot:ylabel "y-data" :replot nil)
    (vgplot:plot x-fit max-ys ";fit stddev upper limit; w l lc rgb \"green\""
	  x-fit min-ys ";fit stddev lower limit; w l lc rgb \"green\""
	  x-fit y-fit ";fit; w l lc rgb \"red\"" 
	  x-data y-data ";data; w p pt 6 ps 2 lc rgb \"black\"")
    params))

(defun walker-plot-data-and-fit-only (the-walker &key (take 1000) (x-column 0) (y-column 1) (fn-number 0))
  (let* ((take (if (or (null take) (null (nthcdr take (funcall the-walker :get-walker))))
		   (walker-get-length the-walker)
		   take))
	 (fn (elt (funcall the-walker :get-fn) fn-number))
	 (data (elt (funcall the-walker :get-data) fn-number))
	 (x-data (elt data x-column))
	 (y-data (elt data y-column))
	 (x-fit (linspace (- (elt x-data 0) 50)
			  (+ (elt x-data (- (length x-data) 1)) 50)
			  :num 1000))
	 (params (walker-get-median-params the-walker take))
	 (y-fit (mapcar #'(lambda (x) (apply fn x params)) x-fit)))
    (vgplot:close-all-plots)
    (vgplot:format-plot t "set terminal qt size 1920,1080 linewidth 3 font \"Arial,18\"")
    (vgplot:format-plot t "set format x \"%.4g\"")
    (vgplot:format-plot t "set format y \"%.4g\"")
    (vgplot:xlabel "x-data" :replot nil)
    (vgplot:ylabel "y-data" :replot nil)
    (vgplot:plot x-fit y-fit ";fit; w l lc rgb \"red\""
	  x-data y-data ";data; w p pt 6 ps 2 lc rgb \"black\"")
    params))

(defun walker-plot-residuals (the-walker &key (take 1000) (x-column 0) (y-column 1) (fn-number 0))
  (let* ((take (if (or (null take) (> take (walker-get-length the-walker)))
		   (walker-get-length the-walker)
		   take))
	 (fn (elt (funcall the-walker :get-fn) fn-number))
	 (data (elt (funcall the-walker :get-data) fn-number))
	 (stddev (elt (funcall the-walker :get-stddev) fn-number))
	 (x-data (elt data x-column))
	 (y-data (elt data y-column))
	 (stddev (if (= 1 (length stddev)) (make-array (length x-data) :initial-element (elt stddev 0)) stddev))
	 (params (walker-get-median-params the-walker take))
	 (y-fit (map 'vector #'(lambda (x) (apply fn x params)) x-data))
	 (y-residuals (map 'vector #'(lambda (yf y) (- yf y)) y-fit y-data)))
    (vgplot:close-all-plots)
    (vgplot:format-plot t "set terminal qt size 1920,1080 linewidth 3 font \"Arial,18\"")
    (vgplot:format-plot t "set format x \"%.4g\"")
    (vgplot:format-plot t "set format y \"%.4g\"")
    (vgplot:xlabel "x-data" :replot nil)
    (vgplot:ylabel "y-data" :replot nil)
    (vgplot:plot x-data y-residuals ";residuals; w p pt 6 ps 2 lc rgb \"black\""
		 x-data stddev ";point stddev; w p pt 2 ps 1 lc rgb \"red\""
		 x-data (make-array (length x-data) :initial-element 0) ";baseline; w l lc rgb \"red\"")))

(defun walker-plot-data (the-walker &key (x-column 0) (y-column 1) (fn-number 0))
  (let* ((data (elt (funcall the-walker :get-data) fn-number))
	 (x-data (elt data x-column))
	 (y-data (elt data y-column)))
    (vgplot:close-all-plots)
    (vgplot:format-plot t "set terminal qt size 1920,1080 linewidth 3 font \"Arial,18\"")
    (vgplot:xlabel "x-data" :replot nil)
    (vgplot:ylabel "y-data" :replot nil)
    (vgplot:plot x-data y-data ";data; w lp pt 6 ps 2 lc rgb \"black\"")))

(defun walker-2d-plot (the-walker keys &optional take)
  (let* ((x-key (first keys))
	 (y-key (second keys))
	 (x (walker-get-param the-walker x-key take))
	 (y (walker-get-param the-walker y-key take)))
    (vgplot:close-all-plots)
    (vgplot:format-plot t "reset")
    (vgplot:format-plot t "set terminal qt size 1920,1080 linewidth 1 font \"Arial,18\"")
    (vgplot:xlabel (format nil "~a" x-key) :replot nil)
    (vgplot:ylabel (format nil "~a" y-key) :replot nil)
    (vgplot:plot x y ";My plot; w p pt 3 ps 4")))

;; TODO only extract each parameter once
(defun walker-all-2d-plots (the-walker &optional take)
  (flet ((permute-params (params)
	   (mapcon
	    #'(lambda (y)
		(butlast (maplist
			  #'(lambda (x)
			      (list (car x) (car (last x))))
			  (reverse y))))
	    params)))
    (let* ((param-keys (walker-get-param-keys the-walker))
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
	(vgplot:format-plot t "set terminal pngcairo size 3840,3840 linewidth 1 font \"Arial,18\"")
	(vgplot:format-plot t "set output \"temp.png\"")
	(vgplot:format-plot t "set format y \"%.4g\"")
	(vgplot:format-plot t "set format x \"%.4g\"")
	(mapcar #'(lambda (x) (inner-plotter x)) key-pairs)
	(vgplot:close-all-plots)))))

(defun walker-param-histo (the-walker key &key (take 10000) (bins 20))
  (let* ((param-values (sort (walker-get-param the-walker key take) #'<))
	 (histo-x (print (make-histo-x param-values bins)))
	 (histo (print (make-histo param-values bins))))
    (vgplot:close-all-plots)
    (vgplot:format-plot t "reset")
    (vgplot:format-plot t "set terminal qt size 1920,1080 linewidth 1 font \"Arial,18\"")
    (vgplot:xlabel (format nil "~a" key) :replot nil)
    (vgplot:ylabel "Counts" :replot nil)
    (vgplot:plot histo-x histo ";Histogram; w l")))

(defun walker-construct-print-list (the-walker &optional take)
  `(:fn ,(mapcar #'sb-kernel:%fun-name (funcall the-walker :get-fn))
    :data ,(funcall the-walker :get-data)
    :param-keys ,(walker-get-param-keys the-walker)
    :stddev ,(funcall the-walker :get-stddev)
    :log-liklihood ,(mapcar #'sb-kernel:%fun-name (funcall the-walker :get-log-liklihood))
    :log-prior ,(mapcar #'sb-kernel:%fun-name (funcall the-walker :get-log-prior))
    :walker ,(walker-get the-walker take)))

(defun walker-save (the-walker filename &optional take)
  (with-open-file (out filename
                       :direction :output
                       :if-exists :supersede)
    (with-standard-io-syntax
      (null (write (walker-construct-print-list the-walker take) :stream out)))))

(defun walker-load (filename &key fn log-liklihood log-prior quiet)
  (with-open-file (in filename
		      :direction :input)
    (with-standard-io-syntax
      (let* ((full (read in))
	     (data (getf full :data))
	     (stddev (getf full :stddev))
	     (walks (getf full :walker))
	     (dummy-params (getf (elt walks 0) :params)))
	(unless quiet
	  (format t "*Recommendations*~%fn: ~s~%log-liklihood: ~s~%log-prior: ~s~%" (getf full :fn) (getf full :log-liklihood) (getf full :log-prior)))
	(when (and fn log-liklihood log-prior)
	  (let ((the-walker (walker-init :fn fn :data data :params dummy-params :stddev stddev :log-liklihood log-liklihood :log-prior log-prior :init nil :trusted t)))
	    (walker-add-walks the-walker walks)
	    the-walker))))))

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
			   (stddev (getf x :stddev))
			   (walks (getf x :walker))
			   (the-walker (walker-init :fn fn :data data :params dummy-params :stddev stddev :log-liklihood log-liklihood :log-prior log-prior :init nil)))
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

(defun file->walker (filename fn params stddev log-liklihood log-prior)
  "Generic function to read file and create walker based on specifying all inputs"
  (let* ((data (file->data-list filename)))
    (walker-init :fn fn :data data :params params :stddev stddev :log-liklihood log-liklihood :log-prior log-prior :init nil)))


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
