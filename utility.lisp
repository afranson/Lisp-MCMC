;;;; Utility functions to help throughout the mcmc-fitting package

(in-package :mcmc-fitting)

;;; Utility
(defmacro repeat (count &body expression)
  `(mapcar #'(lambda (x) (declare (ignore x)) ,@expression) (make-list ,count)))

(defmacro dlambda ((bindings) &body body)
  `(lambda (x)
     (destructuring-bind (,@bindings) x
       ,@body)))

(defun remove-consecutive-duplicates (sequence)
  (mapcon #'(lambda (x)
	      (if (eql (car x) (cadr x)) nil (list (car x))))
	  sequence))

(defun remove-duplicates-plist (plist)
  (let* ((keys (get-plist-keys plist))
	 (unique-keys (remove-duplicates keys)))
    (make-plist
     unique-keys
     (mapcar #'(lambda (x) (getf plist x)) unique-keys))))

(defun elements-between (sequence start-value end-value)
  (remove-if-not #'(lambda (x) (<= start-value x end-value)) sequence))

(defun linspace (start end &key steps (step-size 1))
  (cond (steps (setq step-size (/ (- end start) (- steps 1d0))))
	(step-size (setq steps (+ 1 (floor (- end start) step-size)))))
  (let ((return-list nil))
    (dotimes (i (- steps 1))
      (push (+ start (* step-size i)) return-list))
    (push end return-list)
    (reverse return-list)))

(defun diff (list)
  (mapcon #'(lambda (x) (if (nthcdr 1 x) (list (- (cadr x) (car x))) nil)) list))

(defun diff-matrix (matrix)
  (mapcon #'(lambda (x) (if (nthcdr 1 x) (list (mapcar #'- (cadr x) (car x))) nil)) matrix))

(defun diff-lplist (lplist)
  (let* ((keys (get-plist-keys (elt lplist 0)))
	 (values (mapcar #'(lambda (x) (mapcar #'(lambda (y) (getf x y)) keys)) lplist)))
    (mapcar #'(lambda (x) (make-plist keys x)) (diff-matrix values))))

(defun partition (list n)
  "Pairs elements and removes any with full matches"
  (labels ((rec (l &optional (acc nil))
	     (if (nthcdr (- n 1) l)
		 (rec (subseq l n) (cons (subseq l 0 n) acc))
		 (cons l acc))))
    (reverse (cdr (rec list)))))

(defun flatten (list)
  "Remove all structure from list"
  (labels ((rec (l &optional (acc nil))
	     (if l
		 (if (atom (car l))
		     (rec (cdr l) (cons (car l) acc))
		     (rec (cdr l) (rec (car l) acc)))
		 acc)))
    (reverse (rec list))))

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

