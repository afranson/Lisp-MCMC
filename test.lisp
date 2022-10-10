;;; test.lisp
;;
;; Test the mcmc-fitting package for a lorentzian lineshape

(asdf:load-system :mcmc-fitting)

(in-package :mcmc-fitting)

;; try with data as the name, maybe breaks things?
(defparameter dat (file->data-list "example-data.xls"))

(walker-create woi
	       :fn #'lorder-mixed-bg
	       :data (create-walker-data dat 1 4)
	       :params '(:scale 1d-9 :linewidth 7 :x0 1400
			 :mix 0.9 :bg0 1d-7 :bg1 1d-9)
	       :stddev 1d-7
	       :log-liklihood #'log-liklihood-normal
	       :log-prior #'log-prior-lorder-mixed
	       :docstring "Fully customized walker for fitting a single lorentzian")

;; 6.26 seconds for 1e5 steps
(walker-adaptive-steps woi 100000)
(walker-plot-data-and-fit woi :take 1000)

(defparameter woil (lorder-mixed-bg-walker :data dat :stddev 1d-7 :rows '(0 5)))

(walker-adaptive-steps woil 100000)
(walker-plot-data-and-fit woil :take 1000)

#|
(walker-save woil "walker001.wlk" 1000)
|#
