;;;; package.lisp

(defpackage :mcmc-fitting
  (:nicknames :mfit)
  (:use :cl)
  (:documentation "# mcmc-fitting

Provides an interface for using Markov Chain Monte Carlo to do fitting of various functions.

Makes it easy to generate walker and probability distributions as well as advancing and visualizing the walkers.

GNU General Public License
"))
