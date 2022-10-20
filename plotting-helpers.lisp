;;; Data extraction

(in-package :mcmc-fitting)

(defun 3d-plot-file (filename &key (map t))
  (vgplot:format-plot t "reset")
  (vgplot:format-plot t "set terminal qt size 1920,1080 linewidth 1 font \"Arial,16\"")
  (vgplot:format-plot t "set ticslevel 0")
  (vgplot:format-plot t "set pm3d depthorder")
  (vgplot:format-plot t "set cblabel \"Field Offset (Oe)\"")
  (if map
      (vgplot:format-plot t "set view map")
      (vgplot:format-plot t "unset view"))
  (vgplot:format-plot t "unset key")
  (vgplot:format-plot t "unset grid")
  (vgplot:format-plot t "set pm3d at sb")
  (vgplot:xlabel "X Pos" :replot nil)
  (vgplot:ylabel "Y Pos" :replot nil)
  (vgplot:zlabel "Field Offset (Oe)" :replot nil)
  (vgplot:format-plot t (format nil "splot \"~a\" u 1:2:3 w pm3d" filename)))

(defun 3d-plot-admr-file (filename &key (map t))
  (vgplot:format-plot t "reset")
  (vgplot:format-plot t "set terminal qt size 1080,1080 linewidth 1 font \"Arial,25\"")
  (vgplot:format-plot t "set ticslevel 0")
  (vgplot:format-plot t "set pm3d depthorder")
  (vgplot:format-plot t "set cblabel \"Absorption (dB)\"")
  (if map
      (vgplot:format-plot t "set view map")
      (vgplot:format-plot t "unset view"))
  (vgplot:format-plot t "unset key")
  (vgplot:format-plot t "unset grid")
  (vgplot:format-plot t "set pm3d at sb")
  (vgplot:axis (list -80 80 -80 80))
  (vgplot:xlabel "X Field (Oe)" :replot nil)
  (vgplot:ylabel "Y Field (Oe)" :replot nil)
  (vgplot:zlabel "S21 (a.u.)" :replot nil)
  (vgplot:format-plot t (format nil "splot \"~a\" u 1:2:7 w pm3d, \"\" u (-1*$1):(-1*$2):8 w pm3d" filename)))
