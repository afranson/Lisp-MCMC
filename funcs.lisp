;;; Functions for Fitting and Analyzing Data Trends


(in-package :mcmc-fitting)

;;; Background
(defun bg (x &rest args &key (bg0 0d0) (bg1 0d0) (bg2 0d0) &allow-other-keys)
  "Generic background function of bg0 + bg1*x + bg2*x^2 + ..."
  (declare (sb-ext:muffle-conditions sb-ext:compiler-note)
	   (optimize speed)
	   (double-float x bg0 bg1 bg2))
  (reduce #'+ (mapcar #'(lambda (a e) (declare (double-float a) (fixnum e)) (* a (expt x e))) (list bg0 bg1 bg2) (up-to (length args)))))


;;; Gaussians
(defun gaussian (x &key scale mu sigma &allow-other-keys)
  (declare (sb-ext:muffle-conditions sb-ext:compiler-note)
	   (optimize speed)
	   (double-float x mu sigma scale))
  (* scale (exp (- (expt (/ (- x mu) sigma) 2)))))

(defun double-gaussian-bg (x &key scale mu1 mu2 sigma bg0 &allow-other-keys)
  (declare (sb-ext:muffle-conditions sb-ext:compiler-note)
	   (optimize speed)
	   (double-float x scale mu1 mu2 sigma bg0))
  (+ (the double-float (bg x :bg0 bg0))
     (the double-float (gaussian x :scale scale :mu mu1 :sigma sigma))
     (the double-float (gaussian x :scale scale :mu mu2 :sigma sigma))))



;;; Lorentzians
(defun lorentzian-i (x &key scale linewidth x0 &allow-other-keys)
  "Imaginary part of lorentzian response function, the absorption part"
  (let ((numer (* scale linewidth x))
	(denom (+ (expt (* x linewidth) 2) (expt (- (expt x 2) (expt x0 2)) 2))))
    (/ numer denom)))

(defun lorentzian-r (x &key scale linewidth x0 &allow-other-keys)
  "Real part of lorentzian response function, the phase part"
  (let ((numer (* scale (- (expt x0 2) (expt x 2))))
	(denom (+ (expt (* x linewidth) 2) (expt (- (expt x 2) (expt x0 2)) 2))))
    (+ 1d0 (/ numer denom))))

(defun lorder-i (x &key scale linewidth x0 &allow-other-keys)
  "Derivative of lorentzian absorption function"
  (declare (sb-ext:muffle-conditions sb-ext:compiler-note)
	   (optimize speed)
	   (double-float x scale linewidth x0))
  (let* ((sq-diff (- (expt x 2) (expt x0 2)))
	 (scale-norm (* (* (expt linewidth 2) x0) -0.7698d0))
	 (numer (- (* scale-norm scale linewidth (+ (expt (* x linewidth) 2) (* 4d0 sq-diff (expt x 2)) (- (expt sq-diff 2))))))
	 (denom (expt (+ (expt (* x linewidth) 2) (expt sq-diff 2)) 2)))
    (/ numer denom)))

(defun lorder-r (x &key scale linewidth x0 &allow-other-keys)
  "Derivative of lorentzian phase function"
  (declare (sb-ext:muffle-conditions sb-ext:compiler-note)
	   (optimize speed)
	   (double-float x scale linewidth x0))
  (let* ((sq-diff (- (expt x 2) (expt x0 2)))
	 (scale-norm (* (* (expt linewidth 2) x0) -0.7698d0))
	 (numer (- (* 2d0 scale-norm scale x (- (* linewidth (sqrt (the (double-float 0d0 *)  (- (expt x 2) sq-diff)))) sq-diff) (+ (* linewidth (sqrt (the (double-float 0d0 *) (- (expt x 2) sq-diff)))) sq-diff))))
	 (denom (expt (+ (expt (* x linewidth) 2) (expt sq-diff 2)) 2)))
    (/ numer denom)))

(defun double-lorentzian-bg (x &key scale1 scale2 mu1 mu2 sigma bg0 &allow-other-keys)
  (declare (sb-ext:muffle-conditions sb-ext:compiler-note)
	   (optimize speed)
	   (double-float x scale1 scale2 mu1 mu2 sigma bg0))
  (+ (the double-float (bg x :bg0 bg0))
     (the double-float (lorentzian-i x :scale scale1 :x0 mu1 :linewidth sigma))
     (the double-float (lorentzian-i x :scale scale2 :x0 mu2 :linewidth sigma))))

(defun lorder-mixed (x &key scale linewidth x0 mix &allow-other-keys)
  "Combination of lorder-i and lorder-r with a cosine-like mixing of the two"
  (declare (sb-ext:muffle-conditions sb-ext:compiler-note)
	   (optimize speed)
	   (double-float x scale linewidth x0 mix))
  (+ (* (cos mix) (the double-float (lorder-i x :scale scale :linewidth linewidth :x0 x0)))
     (* (sin mix) (the double-float (lorder-r x :scale scale :linewidth linewidth :x0 x0)))))

(defun lorder-mixed-bg (x &key scale linewidth x0 mix (bg0 0d0) (bg1 0d0) &allow-other-keys)
  "Combination of lorder-i and lorder-r with a cosine-like mixing of the two also with constant and linear background terms"
  (declare (sb-ext:muffle-conditions sb-ext:compiler-note)
	   (optimize speed)
	   (double-float x scale linewidth x0 mix bg0 bg1))
  (+ (the double-float (lorder-mixed x :scale scale :linewidth linewidth :x0 x0 :mix mix))
     (the double-float (bg x :bg0 bg0 :bg1 bg1))))

;;; proves the full and simple match with the same input parameters
;; (defparameter pt (list :scale 20d0 :linewidth 1d0 :x0 500d0 :mix 0d0 :bg0 0d0 :bg1 0d0))
;; (defparameter xt (linspace 0d0 1000 :steps 500000))
;; (defparameter yt (mapcar #'(lambda (x) (apply #'lorder-mixed-bg x pt)) xt))
;; (defparameter yt2 (mapcar #'(lambda (x) (apply #'lorder+bg x pt)) xt))
;; (vgplot:plot xt yt "full" xt yt2 "simple")

(defun example-lorder-params ()
  (list :scale 1d-5 :bg0 1d-6 :linewidth 20d0 :x0 800d0 :bg1 1d-11))

(defun in-plane-fmr-happ-p (phi &key c-w Meff gamma Hu phi-offset %-complete &allow-other-keys)
  "Input phi (in plane rotation angle), output resonant field"
  (let* ((phi (* %-complete (+ phi phi-offset)))
	 (cos2phi (cos (* 2 phi)))
	 (cosphi (cos phi))
	 (gamma2 (expt gamma 2))
	 (t1 (+ (* 4 (expt c-w 2)) (expt (* Meff gamma) 2)))
	 (t2 (* Hu Hu gamma2 (expt cos2phi 2)))
	 (t3 (* cos2phi (- 0 (* 2 Hu Hu gamma2 (expt cosphi 2)) (* 2 Hu Meff gamma2))))
	 (t4 (* (expt cosphi 4) Hu Hu gamma2))
	 (t5 (* (expt cosphi 2) 2 Hu Meff gamma2))
	 )
    (/ (- (sqrt (+ t1 t2 t3 t4 t5))
	  (* Hu gamma cos2phi)
	  (* Hu gamma (expt cosphi 2))
	  (* Meff gamma))
       2 gamma)))

(defun in-plane-fmr-happ-w (w &key c-phi-offset Meff gamma Hu &allow-other-keys)
  "Input frequency, output resonant field"
  (let* ((phi c-phi-offset)
	 (cos2phi (cos (* 2 phi)))
	 (cosphi (cos phi))
	 (gamma2 (expt gamma 2))
	 (t1 (+ (* 4 (expt w 2)) (expt (* Meff gamma) 2)))
	 (t2 (* Hu Hu gamma2 (expt cos2phi 2)))
	 (t3 (* cos2phi (- 0 (* 2 Hu Hu gamma2 (expt cosphi 2)) (* 2 Hu Meff gamma2))))
	 (t4 (* (expt cosphi 4) Hu Hu gamma2))
	 (t5 (* (expt cosphi 2) 2 Hu Meff gamma2))
	 )
    (/ (- (sqrt (+ t1 t2 t3 t4 t5))
	  (* Hu gamma cos2phi)
	  (* Hu gamma (expt cosphi 2))
	  (* Meff gamma))
       2 gamma)))

(defun in-plane-fmr-w (xs &key Meff gamma Hu phi-offset)
  "Input resonant field and phi, output frequency"
  (destructuring-bind (Happ phi) xs
    (let* ((phi (+ phi phi-offset))
	   (t1 (+ Meff Happ (* Hu (expt (cos phi) 2))))
	   (t2 (+ Happ (* Hu (cos (* 2 phi))))))
      (sqrt (* (expt gamma 2) t1 t2)))))

(defun oop-fmr-damping (w &key alpha gamma dH0 &allow-other-keys)
  "Takes in frequency, outputs FWHM linewidth"
  (+ dH0 (/ (* 4 pi alpha w) gamma)))

;; (defun in-plane-fmr-w-error (Happ phi &key Meff gamma Hu Happ-error phi-error)
;;   (let* ((Happc (complex Happ 1e-9))
;; 	 (phic (complex phi 1e-9))
;; 	 (w-Happ-error (* Happ-error (/ (imagpart (in-plane-fmr-w Happc phi :Meff Meff :gamma gamma :Hu Hu)) 1e-9)))
;; 	 (w-phi-error (* phi-error (/ (imagpart (in-plane-fmr-w Happ phic :Meff Meff :gamma gamma :Hu Hu)) 1e-9))))
;;     (sqrt (+ (expt w-happ-error 2) (expt w-phi-error 2)))))
