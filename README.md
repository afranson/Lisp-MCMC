# Lisp MCMC

A framework for doing Markov chain Monte Carlo (MCMC) based fitting with Common Lisp (SBCL).

## Basic Usage

A quick example can be found in test.lisp.

Most functions are softly namespaced by starting with 'walker-', e.g. walker-adaptive-steps, walker-plot-data-and-fit, walker-save, walker-liklihood-plot, etc.

Generally, you will want to abstract the 'walker-init' function to specialize to your specific fitting case (lorentzians, gaussians, whatever function).

You will need to specify several functions before your MCMC journey can begin. First, a fitting function, or more generally, a function that assists in computing the liklihood that you want to maximize:
```
(defun my-function (x &key param0 much-better-param-name1 param2 &allow-other-keys) ... body)
```
Next, a liklihood function (or generally, a 'loss' function) that you want to maximize:
```
(defun log-liklihood-normal-weighted (fn params data error)
  (declare (optimize speed)
	   (cons data)
	   (function fn))
  (let* ((error (if (= 1 (length error)) (make-array (length x) :initial-element (elt error 0)) error)))
    (declare (simple-vector x y error))
    (reduce #'+ (map 'vector #'(lambda (x y z) (log-normal y (apply fn x params) z)) x y error))))
```
The current built-in functions are log-normal and log-poisson.

Finally, a prior which establishes the initial environment of the MCMC simulation, or if your just interested in some basic fitting, just use a flat prior (included as log-prior-flat):
```
(defun my-log-prior (params data)
    (declare (ignore data))
    (prior-bounds-let ((:param0 0 1d4)
                       (:much-better-param-name1 100 150)
                       (:param2 -25 -24))
      bounds-total))
```
Above we make use of the anaphoic macro 'prior-bounds-let' to bound the values of the parameters from our function. The macro generates the 'bounds-total' variable that holds to penalty if the MCMC algorithm attempts to step outside those bounds. IT also generates values such as 'param0-bounds' which holds the penalty for just the :param0 term.

Once all that's together, you can create your MCMC walker:
```
(defvar my-walker)
(setq my-walker (walker-init :fn my-function :data data :params initial-params :stddev stddev :log-liklihood #'log-liklihood-normal-weighted :log-prior #'my-log-prior))
```
and then run it
```
(walker-adaptive-steps my-walker)
```
and then see if it produced anything useful (requires vgplot package and gnuplot with qt)
```
(walker-plot-data-and-fit my-walker)
(walker-plot-liklihood-plot my-walker)
(walker-catepillar-plots my-walker)
```

## Global Parameter Fitting

To do global parameter fitting, you just need to use the walker-create function with a list of datasets, list of error values, list of fitting functions, list of log-liklihoods, and a list of log-priors. That's it. Everything else runs as normal (except walker-plot-data-and-fit, there you have to add the :fn-number kwarg to pick which function/data/fit you want to plot).
