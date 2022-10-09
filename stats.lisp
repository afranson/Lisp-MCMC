;;; Stats

(in-package :mcmc-fitting)

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


;;; binning / histogram
(defun make-histo (sequence &optional num-bins)
  (let* ((bottom (reduce #'min sequence))
	 (top (reduce #'max sequence))
	 (num-bins (if num-bins num-bins (floor (* (- top bottom) (expt (length sequence) 1/3)) (* 2 (iqr sequence)))))
	 (num-bins (+ 1 num-bins))
	 (boundaries (linspace bottom top :num num-bins)))
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
    (linspace (+ bottom (/ gap-size 2)) top :num num-bins)))

