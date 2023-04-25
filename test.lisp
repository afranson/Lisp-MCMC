;;; test.lisp
;;
;; Test the mcmc-fitting package for a lorentzian lineshape

(require 'mcmc-fitting)

(in-package :mcmc-fitting)

;; Quickly gather a list of the files you want to process
(print (get-filename "." :include '("example" ".xls") :exclude '("test")))

(defparameter data (read-file->data "example-data.xls"))

(defparameter woi (walker-create :function #'lorder-mixed-bg
				 :data (create-walker-data data 1 4)
				 :params '(:scale 1d-5 :linewidth 7 :x0 2200
					   :mix 0.9 :bg0 1d-7 :bg1 1d-9)
					; Also can give a list of stddev values for each individual data point
				 :data-error 1d-7
				 :log-liklihood #'log-liklihood-normal
				 :log-prior #'log-prior-lorder-mixed))

;; 6.34 seconds for 1e5 steps
(walker-adaptive-steps woi 100000)
(print (walker-plot-data-and-fit woi :take 1000))
;; #S(THIS-WALKER-STEP
   ;; :PROB 4646.756030280576d0
   ;; :PARAMS (:SCALE -4.788638538682475d-6 :LINEWIDTH 121.09571484294366d0 :X0
   ;;          2784.6836516658504d0 :MIX 3.141546812249173d0 :BG0
   ;;          -1.0629009389997092d-6 :BG1 2.8207485034278606d-10))
(format t "Q Factor: ~,2e~%" (walker-with-exp woi '(/ :linewidth :x0) :take 1000))

(defparameter woil (lorder-mixed-bg-walker :data data :data-error 1d-7 :rows '(0 4)))

(walker-adaptive-steps woil 100000)
(print (walker-plot-data-and-fit woil :take 1000))

;;; Saving

;; (walker-save woil "walker001.wlk" 1000)

;; with no supplied functions, loading will suggest what functions to use (saving doesn't currently serialize functions, just saves their string names)
;; (walker-load "walker001.wlk")

;; with functions, it returns the walker	
;; (defparameter loaded-walker (walker-load "walker001.wlk"
;; 					 :function #'lorder-mixed-bg
;; 					 :log-liklihood #'log-liklihood-normal
;; 					 :log-prior #'log-prior-lorder-mixed))


;;; Global fit
;; Shares the linewidth, x0, and mix parameters with the regular lorder-mixed-bg
(defun lorder-mixed-bg2 (x &key scale2 linewidth x0 mix (bg02 0d0) (bg12 0d0) &allow-other-keys)
  (lorder-mixed-bg x :scale scale2 :linewidth linewidth :x0 x0 :mix mix :bg0 bg02 :bg1 bg12))

;; Just create the normal walker and double, triple, etc. everything up for global fits of multiple data sets
(defparameter woig (walker-create :function (list #'lorder-mixed-bg
						  #'lorder-mixed-bg2)
				  :data (list (create-walker-data data 1 4)
					      (create-walker-data data 1 5))
				  :params '(:scale 1d-6 :linewidth 100 :x0 2700
					    :mix 0.1 :bg0 1d-7 :bg1 1d-10
					    :scale2 1d-8 :bg02 1d-7 :bg12 1d-10)
				  :data-error (list (list 1d-7)
						    (list 1d-7))
				  :log-liklihood (list #'log-liklihood-normal
						       #'log-liklihood-normal)
				  :log-prior (list #'log-prior-lorder-mixed
						   #'log-prior-lorder-mixed)))

(walker-adaptive-steps woig 100000)
(walker-plot-data-and-fit woig :take 1000 :fn-number 0)
(walker-plot-data-and-fit woig :take 1000 :fn-number 1)
(walker-catepillar-plots woig)
(show)
(walker-all-2d-plots woig 1000) ;; work in progress
(show)
