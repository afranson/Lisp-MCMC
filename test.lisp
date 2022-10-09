;;; test.lisp
;;
;; Test the mcmc-fitting package

(asdf:load-system :mcmc-fitting)

(in-package :mcmc-fitting)

;; try with data as the name, maybe breaks things?
(defparameter dat (file->data-list "example-data.xls"))

(walker-create woi
	       :fn #'lorder-mixed-bg
	       :data (list (list (make-array (length (elt dat 1)) :initial-contents (elt dat 1)) (make-array (length (elt dat 4)) :initial-contents (elt dat 4))))
	       :params '(:scale 1d-9 :linewidth 7 :x0 1400 :mix 0.9 :bg0 1d-7 :bg1 1d-9)
	       :stddev 1d-7
	       :log-liklihood #'log-liklihood-normal
	       :log-prior #'log-prior-lorder-mixed)

(walker-adaptive-steps woi)

;; 6.26 
