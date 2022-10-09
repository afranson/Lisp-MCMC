;;;; Utility functions to help throughout the mcmc-fitting package

(in-package :mcmc-fitting)

;;; Utility
(defmacro repeat (count &body expression)
  `(mapcar #'(lambda (x) (declare (ignore x)) ,@expression) (make-list ,count)))

(defmacro dlambda ((bindings) &body body)
  `(lambda (x)
     (destructuring-bind (,@bindings) x
       ,@body)))

(defun mkstr (&rest args)
  (with-output-to-string (s)
    (dolist (a args) (princ a s))))

(defun symb (&rest args)
  (values (intern (apply #'mkstr args))))

(defun symb-keyword (&rest args)
  (values (intern (apply #'mkstr args) :keyword)))

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

(defun linspace (start end &key (num 50))
  (let ((step (/ (- end start) num)))
    (do ((n 0 (1+ n))
	 (ret (list start) (push (+ (car ret) step) ret)))
	((= n num) (reverse ret)))))

(defun n-nums (n)
  (do ((m 0 (1+ m))
       (l nil (push m l)))
      ((= m n) (reverse l))
    (declare (fixnum m))))

(defun up-to (n &optional (start 0))
  (assert (> n start))
  (do ((m start (1+ m))
       (l nil (push m l)))
      ((= m (1+ n)) (reverse l))
    (declare (fixnum m))))

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
