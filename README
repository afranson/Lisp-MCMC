# Lisp MCMC

A framework for doing Markov chain Monte Carlo (MCMC) based fitting with Common Lisp.

## Basic Usage

Most functions are softly namespaced by starting with 'walker-', e.g. walker-adaptive-steps, walker-plot-data-and-fit, walker-save, walker-liklihood-plot, etc. The significant anomoly is 'create-walker', to initialize a walker.

Generally, you will want to abstract the 'create-walker' function to specialize to your specific fitting case (lorentzians, gaussians, whatever function).

First you'll need to create your fitting function. The usual format will be like:
```
	(defun gaussian (x &key scale mu sigma &allow-other-keys) ... body)
```

Then you likely need to generate a log-liklihood function and a log prior. For simple data fitting (x, y data, gaussian error), the default functions 'log-liklihood-normal', and 'log-prior-flat' will usually be fine as long as your error is reasonable and your initial parameters are 'close'.

If you need to set limits on your fitting parameters, that is done with the 'prior-bounds-let' macro. This works similar to a regular 'let' expression, but instead of (let ((varname expr)...) body), it works like (prior-bounds-let ((symbol min-value max-value)...) body).

Once that's all squared away, use your walker object via (walker-adaptive-steps walker-object) to advance the MCMC algorithm and then (walker-liklihood-plot walker-object) and (walker-plot-data-and-fit walker-object) to visualize the progress of the walker. Then run stepping function repeatedly until you find the algorithm has stabilized.

## Global Parameter Fitting

To do global parameter fitting, you just need to use the create-walker function with a list of datasets, list of error values, list of fitting functions, list of log-liklihoods, and a list of log-priors. That's it. Everything else runs as normal (except walker-plot-data-and-fit, there you have to add the :fn-number kwarg to pick which function/data/fit you want to plot).