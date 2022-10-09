;;; test.lisp
;;
;; Test the mcmc-fitting package

(asdf:load-system :mcmc-fitting)

(in-package :mcmc-fitting)

(defparameter data (file->data-list "example-data.xls"))

(walker-create woi
	       :fn #'lorder-mixed-bg
	       :data (list (list (make-array (length (elt data 1)) :initial-contents (elt data 1)) (make-array (length (elt data 4)) :initial-contents (elt data 4))))
	       :params '(:scale 4d-6 :linewidth 70d0 :x0 2800d0 :mix 0d0 :bg0 1d-9 :bg1 1d-10)
	       :stddev '(1d-7)
	       :log-liklihood #'log-liklihood-normal
	       :log-prior #'log-prior-flat)

(walker-adaptive-steps woi)

;;; errors concerning the depth of the data in the log-liklihood-normal function
