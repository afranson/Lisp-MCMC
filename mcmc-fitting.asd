;;;; mcmc-fitting.asd

(asdf:defsystem :mcmc-fitting
  :serial t
  :description "MCMC-Based Library for Fitting Data"
  :author "Andrew Franson"
  :license "GPL"
  :depends-on (:vgplot)
  :components ((:file "package")
	       (:file "utility")
	       (:file "funcs")
               (:file "vector-mcmc")
	       (:file "fmr-specific")
	       (:file "nv-specific")))

;; (asdf:defsystem :mcmc-fitting-test
;;   :serial t
;;   :description "Test environment for mcmc-fitting"
;;   :author "Andrew Franson"
;;   :license "GPL"
;;   :depends-on (:mcmc-fitting :lisp-unit)
;;   :components ((:file "package-test")
;;                (:file "vgplot-test")))