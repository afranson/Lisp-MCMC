;;;; mcmc-fitting.asd

(asdf:defsystem :mcmc-fitting
  :serial t
  :description "MCMC-Based Library for Fitting Data"
  :author "Andrew Franson"
  :license "GPLv3"
  :depends-on (:cl-gnuplot
	       :alexandria)
  :components ((:file "package")
               (:file "mcmc-fitting")))
