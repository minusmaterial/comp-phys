(in-package #:comp-phys)

;;A simple simulation of a rocket flying to the moon.
;;all units in m, s, kg.
;;Vectors are in 3d for posterity, but all the action is happening in the x-y plane.

(defparameter *big-g* 6.6726e-11)
(defparameter *mass-earth* 5.9742e24)
(defparameter *mass-moon* 7.35e22)
(defparameter *d-earth-moon* 3.844e8)
(defparameter *orb-rad-earth* 
  (* *d-earth-moon* 
     (/ *mass-moon*
        (+ *mass-moon*
           *mass-earth*))))

(defparameter *orb-angular-v-earth* 
    (sqrt (/ (* *big-g* 
                (+ *mass-moon* 
                   *mass-earth*))
             (expt *orb-rad-earth* 3))))

(defparameter *orb-rad-moon* 
    (* *d-earth-moon* 
      (/ *mass-earth*
         (+ *mass-moon*
             *mass-earth*))))

(defparameter *orb-angular-v-moon*
    (sqrt (/ (* *big-g* 
                (+ *mass-moon* 
                   *mass-earth*))
             (expt *orb-rad-moon* 3))))

(defparameter *earth-init-ang* 0)
(defparameter *moon-init-ang* pi)

(defparameter *rocket* (make-instance 'body
                                      'pos (v! 2e8 0 0)
                                      'v (v! 0 1000 0)
                                      'mass 1000))

(defclass body nil
    ((pos :initarg pos
          :type rtg-math.base-vectors:v3!
          :accessor body-pos)
     (v :initarg v
        :type rtg-math.base-vectors:v3!
        :accessor body-v)
     (mass :initarg mass
           :type long-float
           :accessor body-mass)))

(defun iterate-body (body dt &rest forces)
    (when (equal nil forces)
      (setf forces `(,(v! 0 0 0))))
    (let ((net-force (v:/ (apply #'v:+ forces) (body-mass body))))
      net-force))

(defun runge-kutta (x v a dt)
    (let* ((z1 (+ x (* (/ dt 2) v)))
           (z1-dot (+ x (* (/ dt 2) a)))
           )))

(defun d (vec1 vec2)
    (rtg-math.vectors:distance vec1 vec2))


;;;In this proof-of-concept version, this is a toy model.  Only the rocket is 'properly' modelled as a full body; the earth and moon are analytically-derived smoke and mirrors.
(defun pos-earth (cur-time)
       (v:* (v! 
              (cos (+ (* cur-time 
                         *orb-angular-v-earth*)
                         *earth-init-ang*))
              (sin (+ (* cur-time 
                         *orb-angular-v-earth*)
                         *earth-init-ang*))
              0e0)
            *orb-rad-earth*))

(defun pos-moon (cur-time)
       (v:* (v! 
              (cos (+ (* cur-time 
                         *orb-angular-v-moon*)
                         *moon-init-ang*))
              (sin (+ (* cur-time 
                         *orb-angular-v-moon*)
                         *moon-init-ang*))
              0e0)
            *orb-rad-moon*))



(defun plot-planets (&optional (start 0) 
                               (end (* 3600 24 30)) 
                               (dt 1) 
                               (plot-steps 1000))
    (let* (
           (timevec (range start end (/ (- end start) plot-steps)))
           (sim-steps (round (/ (- end start) dt)))
           (plot-interval (round (/ sim-steps plot-steps )))
           (e-x (make-array plot-steps))
           (e-y (make-array plot-steps))
           (m-x (make-array plot-steps))
           (m-y (make-array plot-steps))
           )
      (iter (for i from 0 to (1- sim-steps))
            (when (= (mod i plot-interval) 0)
              (let* ((j (/ i plot-interval))
                   (e-v (pos-earth (aref timevec j)))
                   (m-v (pos-moon (aref timevec j))))
                (setf (aref e-x j) (x e-v))
                (setf (aref e-y j) (y e-v))
                (setf (aref m-x j) (x m-v))
                (setf (aref m-y j) (y m-v)))))
      (progn 
        (3d-plot e-x e-y timevec "earth" m-x m-y timevec "moon" )
        (format-plot nil 
                     "set xrange[-~a:~a]" 
                     *orb-rad-moon* 
                     *orb-rad-moon*)
        (format-plot nil 
                     "set yrange[-~a:~a]"
                     *orb-rad-moon*
                     *orb-rad-moon*)
        ;(format-plot nil "set size ratio -1")
        ;(format-plot nil "set size square")
        )
      ))


