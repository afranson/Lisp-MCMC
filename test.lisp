;;; test.lisp
;;
;; Test the mcmc-fitting package for a lorentzian lineshape

(asdf:load-system :mcmc-fitting)

(in-package :mcmc-fitting)

;; Quickly gather a list of the files you want to process
(print (get-filename "." :include '("example" ".xls") :exclude '("test")))

;; try with data as the name, maybe breaks things?
(defparameter data (read-file->data "example-data.xls"))

(walker-create woi
	       :fn #'lorder-mixed-bg
	       :data (create-walker-data data 1 4)
	       :params '(:scale 1d-5 :linewidth 7 :x0 1400
				:mix 0.9 :bg0 1d-7 :bg1 1d-9)
               ; Also can give a list of stddev values for each individual data point
	       :stddev 1d-7
	       :log-liklihood #'log-liklihood-normal
	       :log-prior #'log-prior-lorder-mixed
	       :docstring "Fully customized walker for fitting a single lorentzian")

;; 6.26 seconds for 1e5 steps
(walker-adaptive-steps woi 100000)
(walker-plot-data-and-fit woi :take 1000)
(format t "Q Factor: ~,2e~%" (walker-get-f woi (/ :linewidth :x0) :take 1000))

(defparameter woil (lorder-mixed-bg-walker :data data :stddev 1d-7 :rows '(0 5)))

(walker-adaptive-steps woil 100000)
(walker-plot-data-and-fit woil :take 1000)

;;; Saving

;; (walker-save woil "walker001.wlk" 1000)

;; with no supplied functions, loading will suggest what functions to use (saving doesn't currently serialize functions, just saves their string names)
;; (walker-load "walker001.wlk")

;; with functions, it returns the walker	
;; (defparameter loaded-walker (walker-load "walker001.wlk"
;; 					 :fn #'lorder-mixed-bg
;; 					 :log-liklihood #'log-liklihood-normal
;; 					 :log-prior #'log-prior-lorder-mixed))


;;; Global fit
;; Shares the linewidth, x0, and mix parameters with the regular lorder-mixed-bg
(defun lorder-mixed-bg2 (x &key scale2 linewidth x0 mix (bg02 0d0) (bg12 0d0) &allow-other-keys)
  (lorder-mixed-bg x :scale scale2 :linewidth linewidth :x0 x0 :mix mix :bg0 bg02 :bg1 bg12))

;; Just create the normal walker and double, triple, etc. everything up for global fits of multiple data sets
(walker-create woig
	       :fn (list #'lorder-mixed-bg
			 #'lorder-mixed-bg2)
	       :data (list (create-walker-data data 1 4)
			   (create-walker-data data 1 5))
	       :params '(:scale 1d-6 :linewidth 100 :x0 2700
			 :mix 0.1 :bg0 1d-7 :bg1 1d-10
			 :scale2 1d-8 :bg02 1d-7 :bg12 1d-10)
	       :stddev (list (list 1d-7)
			     (list 1d-7))
	       :log-liklihood (list #'log-liklihood-normal
				    #'log-liklihood-normal)
	       :log-prior (list #'log-prior-lorder-mixed
				#'log-prior-lorder-mixed)
	       :docstring "Fully customized walker for fitting two lorentzians with shared parameters")

(walker-adaptive-steps woig 100000)
(walker-plot-data-and-fit woig :take 1000 :fn-number 0)
(walker-plot-data-and-fit woig :take 1000 :fn-number 1)
(walker-catepillar-plots woig)
(show)
(walker-all-2d-plots woig 1000) ; work in progress
(show)
