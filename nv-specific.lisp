;;;; nv magnetometry specific walker fitting stuff

(in-package :mcmc-fitting)

(defun nv-data->separated (data)
  (mapcar #'(lambda (x) (list (elt data 0) x)) (subseq data 1)))

(defun nv-dir->data (directory)
  (let ((files (uiop:directory-files directory)))
    (mapcan #'(lambda (x) (nv-data->separated (read-file->data x :delim #\;))) files)))

(defun log-liklihood-nv (fn params data error)
  (declare (optimize speed)
	   (cons data)
	   (function fn))
  (let ((x (elt data 0))
	(y (elt data 1)))
    (declare (simple-vector x y))
    (reduce #'+ (map 'vector #'(lambda (x y) (log-normal (apply fn x params) error y)) x y))))

(defun log-prior-nv (params data)
  (declare (sb-ext:muffle-conditions sb-ext:compiler-note)
	   (optimize speed)
	   (ignorable data))
  (prior-bounds-let ((:scale1 1d-5 1d1)
		     (:scale2 1d-5 1d1)
	  	     (:mu1 2850 2870)
		     (:mu2 2870 2890)
		     (:sigma 9 20)
		     (:bg0 0 1d-5))
    (+ bounds-total
       (if (> mu1 mu2) -1e9 0e0)
       (if (< (- mu2 mu1) 6) -1e9 0e0)
       (if (not (< 0.9 (/ scale1 scale2) 1.1)) -1e9 0e0))))

(defun nv-data-std-dev (data)
  (let* ((y-data (elt data 1))
	 (y-length/10 (floor (length y-data) 10))
	 (y-start-std-dev (standard-deviation (subseq y-data 0 y-length/10)))
	 (y-end-std-dev (standard-deviation (subseq y-data (- (length y-data) y-length/10)))))
    (min y-start-std-dev y-end-std-dev)))

(defun guess-nv-params (data)
  (let* ((y (elt data 1))
	 (y-max (reduce #'max y))
	 (y-min (reduce #'min y))
	 (scale (/ (- y-max y-min) 4.4d-5)))
    (list :scale1 scale :scale2 scale :mu1 2863d0 :mu2 2873d0 :sigma 10d0 :bg0 (float y-min 0d0))))

(defun nv-walker (data)
  (walker-init :fn #'double-lorentzian-bg
	       :data data
	       :params (guess-nv-params data)
	       :stddev (nv-data-std-dev data)
	       :log-liklihood #'log-liklihood-nv
	       :log-prior #'log-prior-nv))

(defun dir->nv-walkers (dir)
  (let ((walkers (mapcar #'nv-walker (nv-dir->data dir))))
    (mapc #'(lambda (x) (walker-adaptive-steps x)) walkers)
    walkers))

(defun file->nv-walkers (filename)
  (let ((walkers (mapcar #'nv-walker (nv-data->separated (read-file->data filename :delim #\;)))))
    (mapc #'(lambda (x) (walker-adaptive-steps x)) walkers)
    walkers))

(defun walker-field-offset (the-walker)
  (walker-with-exp the-walker '(/ (- :mu2 :mu1) 2 2.8) :take 1000))

;;; TODO problem for another day
;; (defun walker-set-field-offset (the-walker-set background-splitting)
;;   (walker-set-get-f the-walker-set (- (/ (- :mu2 :mu1) 2 2.8) background-splitting) 1000))


(defmacro walker-set-make-file-3d-plot-exp (the-walker-set exp row-length &optional file-out take)
  (let ((walk-set (gensym))
	(row-len (gensym))
	(xs (gensym))
	(ys (gensym))
	(field-offsets (gensym))
	(filename-out (gensym)))
    `(let* ((,walk-set ,the-walker-set)
	    (,row-len ,row-length)
	    (,xs (mapcar #'(lambda (x) (mod x ,row-len)) (linspace 0 (- (length ,walk-set) 1))))
	    (,ys (mapcar #'(lambda (x) (floor x ,row-len)) (linspace 0 (- (length ,walk-set) 1))))
	    (,field-offsets (mapcar #' eval (walker-set-get-f ,walk-set ',exp ,take)))
	    (,filename-out ,file-out)
	    (,filename-out (if ,filename-out ,filename-out "./3d-temp-file.txt")))
       (with-open-file (out ,filename-out :direction :output :if-exists :supersede)
	 (mapcar #'(lambda (x)
		     (format out "~f ~f ~f~%" (elt x 0) (elt x 1) (elt x 2))
		     (when (= (elt x 0) (- ,row-len 1))
		       (terpri out)))
		 (mapcar #'list ,xs ,ys ,field-offsets))))))


(defun nv-pretty-heatmap (&key (map nil) (cbar-range '(0 "*")) (z-range '(-5 "*")))
  (plt:send-command :xlabel "\"X Pos\""
		    :ylabel "\"Y Pos\""
		    :zlabel "\"Field Offset (Oe)\" rotate parallel"
		    :cbrange (format nil "[~a:~a]" (elt cbar-range 0) (elt cbar-range 1))
		    :zrange (format nil "[~a:~a]" (elt z-range 0) (elt z-range 1))
		    :view (if map "map" "unset"))
  (plt:replot))
