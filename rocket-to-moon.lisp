(in-package #:comp-phys)

;;A simple simulation of a rocket flying to the moon.
;;all units in m, s, kg.
;;Vectors are in 3d for posterity, but all the action is happening in the x-y plane.

(defparameter *big-g* 6.6726d-11)
(defparameter *mass-earth* 5.9742d24)
(defparameter *rad-earth* 6.3781d6)
(defparameter *mass-moon* 7.35d22)
(defparameter *rad-moon* 1738.1d3)
(defparameter *d-earth-moon* 3.844d8)

(defun days (x)
    (* 3600 24 x))

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

(defun pos-moon (cur-time)
       (v* (vec 
              (cos (+ (* cur-time 
                         *orb-angular-v-moon*)
                         *moon-init-ang*))
              (sin (+ (* cur-time 
                         *orb-angular-v-moon*)
                         *moon-init-ang*))
              0d0)
            *orb-rad-moon*))

(defun pos-earth (cur-time)
       (v* (vec 
              (cos (+ (* cur-time 
                         *orb-angular-v-earth*)
                         *earth-init-ang*))
              (sin (+ (* cur-time 
                         *orb-angular-v-earth*)
                         *earth-init-ang*))
              0d0)
            *orb-rad-earth*))

(defun v-earth (cur-time)
  (v* (vec 
              (* -1 *orb-angular-v-earth*
                 (sin (+ (* cur-time 
                         *orb-angular-v-earth*)
                         *earth-init-ang*)))
              (* *orb-angular-v-earth*
                 (cos (+ (* cur-time 
                         *orb-angular-v-earth*)
                         *earth-init-ang*)))
              0d0)
            *orb-rad-earth*))

(defparameter *d-moon-l2*
    (* *orb-rad-moon*
       (expt (/ *mass-moon*
                (* 3 *mass-earth*))
             1/3)))

(defparameter *new-d-moon-l2*
    (+ *d-moon-l2* 3.739046d6))

(defparameter *initial-l2-pos* 
    (v+ (pos-moon 0) 
       (v* *d-moon-l2* 
          (vunit (pos-moon 0)))))

(defparameter *adjusted-l2-pos*
    (v+ (pos-moon 0) 
       (v* *new-d-moon-l2*
          (vunit (pos-moon 0)))))

;(defparameter *made-up-l2-pos* (vec3 0 -4.405093102104164d8 0))
(defparameter *made-up-l2-pos* (vec3  -4.412282404881194d8 0 0))

(defstruct (body 
            (:constructor %make-body)) 
    (pos (vec 0 0 0) :type vec3 )
    (v (vec 0 0 0) :type vec3 )
    (mass 0d0 :type long-float))

(defun make-body (pos v mass)
    (setf mass (coerce mass 'double-float))
    (setf pos (vec (vx pos) (vy pos) (vz pos)))
    (setf v (vec (vx v) (vy v) (vz v)))
    (%make-body :pos pos :v v :mass mass)
    )

(defparameter *rocket* (make-body 
                         (vec -4.54038613d8
                              0
                              0)
                         (vec 0 
                              (* *orb-rad-moon* 
                                 *orb-angular-v-moon*
                                 -1)
                              0)
                         1000))

(defparameter *rocket-2* (make-body 
                           (v+ (pos-earth 0) (v* (vunit (pos-earth 0)) (+ *rad-earth* 200d3)))
                           ;(v/ (v- (pos-earth 0.1) (pos-earth 0)) 0.1)
                           (v* (v-earth 0) 1)
                           

                           1000))

(defun lagrange-approx (r &optional 
                          (m1 *mass-earth*) 
                          (m2 *mass-moon*) 
                          (big-r *d-earth-moon*))
    (- (+ (/ m1 
             (expt (+ r big-r) 2))
          (/ m2
             (expt r 2)))
       (+ (/ m1 
             (expt big-r 2))
          (/ (* r
                (+ m1 m2))
             (expt big-r 3)))))

(defun simple-acceleration (cur-time)
    (setf cur-time (coerce cur-time 'double-float))
    (lambda (pos dt) 
      (setf dt (coerce dt 'double-float))
      (funcall 
        (lambda (pos dt) 
          (declare (type double-float dt cur-time)
                   (type vec3  pos))
          (let ((vec-to-earth (v- (pos-earth (+ cur-time dt)) 
                                              pos))
                (vec-to-moon (v- (pos-moon (+ cur-time dt)) 
                                            pos)))
            (v+ (v/ (v* (vunit vec-to-earth) 
                           (* *big-g* *mass-earth*))
                      (expt (vlength vec-to-earth) 2))
                 (v/ (v* (vunit vec-to-moon) 
                           (* *big-g* *mass-moon* ))
                      (expt (vlength vec-to-moon) 2))))) 
        pos dt)))

(defun iterate-body-r4k (body dt accel-fun)
   ;Force-fun should take one argument, a number between 0 and dt.
   ;It's needed to do r4k by recalculating the forces applied at various advanced timesteps.
   (declare (optimize (safety 1) (speed 3)) )
   (declare (type function accel-fun )
            (type double-float dt)
            (type body body))
    (let* ((half-dt (/ dt 2))
           (orig-v (body-v body))
           (orig-a (funcall accel-fun (body-pos body) 0d0))
           (z1 (v+ (body-pos body) 
                    (v* (body-v body) 
                         half-dt)))
           (z1-dot (v+ (body-v body) 
                        (v* orig-a
                             half-dt)))
           (z1-a (funcall accel-fun z1 half-dt))
           (z2 (v+ (body-pos body) 
                    (v* z1-dot 
                         half-dt)))
           (z2-dot (v+ (body-v body) 
                        (v* z1-a
                             half-dt)))
           (z2-a (funcall accel-fun z2 half-dt))
           (z3 (v+ (body-pos body) 
                    (v* z2-dot 
                         dt)))
           (z3-dot (v+ (body-v body) 
                        (v* z2-a
                             dt)))
           (z3-a (funcall accel-fun z3 dt))
           )
      ;(print "BEGIN DIAGNOSTICS")
      ;(print (list z1 z1-dot z1-a z2 z2-dot z2-a z3 z3-dot z3-a))
      ;(print "END DIAGNOSTICS")

      (setf (body-v body)
            (v+ (body-v body)
               (v* (/ dt 6)
                  (v+ orig-a
                     (v* 2 z1-a)
                     (v* 2 z2-a)
                     z3-a))))
      (setf (body-pos body)
            (v+ (body-pos body)
               (v* (/ dt 6)
                  (v+ orig-v
                     (v* 2 z1-dot)
                     (v* 2 z2-dot)
                     z3-dot))))
      
      ))

;;;In this proof-of-concept version, this is a toy model.  Only the rocket is 'properly' modelled as a full body; the earth and moon are analytically-derived smoke and mirrors.




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
                (setf (aref e-x j) (vx e-v))
                (setf (aref e-y j) (vy e-v))
                (setf (aref m-x j) (vx m-v))
                (setf (aref m-y j) (vy m-v)))))
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
        (format-plot nil "set view equal xy")
        ;(format-plot nil "set size square")
        )
      ))

(defun sim-planets-and-rocket (&optional (rock-orig *rocket*)
                                         (start 0) 
                                         (end (* 3600 24 7)) 
                                         (dt 20) 
                                         (n-plot-samples 1000))
    ;(declare (optimize (safety 0) (speed 3)) )
    (close-all-plots)
    (let* (
           (rocket (copy-structure rock-orig))
           (timevec (range start end (/ (- end start) n-plot-samples)))
           (sim-steps (round (/ (- end start) dt)))
           (plot-ratio (/ sim-steps n-plot-samples))
           (plot-coords
             (nreverse
               (cons `(,(1+ sim-steps) ,(1+ n-plot-samples))
                 (nreverse (iter 
                             (for i from 0 to (1- n-plot-samples))
                             (collecting (list 
                                           (min 
                                             (round (* i plot-ratio))
                                             sim-steps)
                                           i)))))))
           (e-x (make-array n-plot-samples))
           (e-y (make-array n-plot-samples))
           (m-x (make-array n-plot-samples))
           (m-y (make-array n-plot-samples))
           (r-x (make-array n-plot-samples))
           (r-y (make-array n-plot-samples))
           (d-r-m (make-array n-plot-samples))
           
           )
      
      
      (iter (for i from 0 to (1- sim-steps))
            (let* ((curtime (* dt i)))

             (when (= i (caar plot-coords))
               ;(print i)
              (let* ((j (cadr (pop plot-coords)))
                   (e-v (pos-earth (aref timevec j)))
                   (m-v (pos-moon (aref timevec j))))
                (setf (aref e-x j) (vx e-v))
                (setf (aref e-y j) (vy e-v))
                (setf (aref m-x j) (vx m-v))
                (setf (aref m-y j) (vy m-v))
                (setf (aref r-x j) (vx (body-pos rocket)))
                (setf (aref r-y j) (vy (body-pos rocket)))
                (setf (aref d-r-m j) (vdistance (body-pos rocket)
                                                e-v))
                ))
              (iterate-body-r4k rocket 
                                (coerce dt 'double-float) 
                                (simple-acceleration curtime))
              (when (or (< (vdistance (body-pos rocket)
                                      (pos-moon curtime))
                           *rad-moon*)
                        (< (vdistance (body-pos rocket)
                                      (pos-earth curtime))
                           *rad-earth*))
                    (progn (format t "COLLISION!!!~%")
                           (iter:finish)))))
      
      (progn 
        (vgplot:close-all-plots)
        ;(3d-plot e-x e-y timevec "earth" m-x m-y timevec "moon" r-x r-y timevec "rocket")
        (3d-plot e-x e-y timevec "earth" r-x r-y timevec "rocket")
        ;(let ((size (coerce *orb-rad-moon* 'short-float)))
         ; (format-plot nil 
         ;             "set xrange[-~a:~a]" 
         ;             size 
         ;             size)
         ;(format-plot nil 
         ;             "set yrange[-~a:~a]"
         ;             size 
         ;             size)
         ;) 
         (format-plot nil "set view equal xy")
         (xlabel "test")
         ;(format-plot nil "set size ratio -1")
         ;(format-plot nil "set size square")
         )
       (read )
      (progn
        (vgplot:close-all-plots)
        ;(setf timevec (map 'vector 
                           ;(lambda (x) (/ x (/ (* 2 pi) *orb-angular-v-moon*)))
                           ;timevec))
        (plot timevec d-r-m timevec (make-array (list (length timevec)) :initial-element 100000d3)))
      (list (subseq r-x (- (length r-x) 200) (length r-x))
            (subseq r-y (- (length r-y) 200) (length r-y))
            )
      ))
