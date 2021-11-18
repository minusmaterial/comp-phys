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
(defparameter *orb-rad-moon* 
    (* *d-earth-moon* 
      (/ *mass-earth*
         (+ *mass-moon*
             *mass-earth*))))

(defparameter *orb-angular-v-earth* 
    (sqrt (/ (* *big-g* 
                (+ *mass-moon* 
                   *mass-earth*))
             (expt (+ *orb-rad-earth* *orb-rad-moon* ) 3))))


(defparameter *orb-angular-v-moon*
    (sqrt (/ (* *big-g* 
                (+ *mass-moon* 
                   *mass-earth*))
             (expt (+ *orb-rad-earth* *orb-rad-moon*) 3))))

(assert (equal *orb-angular-v-moon* *orb-angular-v-earth*))

(defparameter *orb-period-moon*
    (/ (* 2 pi)
       *orb-angular-v-moon*))

(defparameter *earth-init-ang* 0)
(defparameter *moon-init-ang* pi)

(defparameter *ldist* 6.451858358131048d7 )

(defparameter *e-m-l2-k*
    (/ (+ *orb-rad-moon* *ldist*)
       (+ *orb-rad-moon* *orb-rad-earth*)))

(defparameter *e-m-a* 
    (a (mu *mass-earth* *mass-moon*)
       *e-m-l2-k*))
(defparameter *e-m-b* 
    (b (mu *mass-earth* *mass-moon*)
       *e-m-l2-k*))

(defparameter *e-m-mat*
  (3d-matrices:mat 
    0 0 1 0
    0 0 0 1
    *e-m-a* 0 0 2
    0 (* -1 *e-m-b*) -2 0))

;m is rows, n is columns
;Also, YOU ENTER VALUES COLUMN BY COLUMN!!!!!!
;access with (ref mat nrow ncol)
;(defparameter *e-m-mat-working*
;    (magicl:make-complex-matrix 
;      4 
;      4 
;      (list 0 0 *e-m-a* 0
;            0 0 0 (* -1 *e-m-b*)
;            1 0 0 -2
;            0 1 2 0
;            )))

;(defun weird-to-rot-mat (m)
;    (multiple-value-bind (vals vects)
;      (magicl:eig m)
;    (let* ()
;      
;      (magicl:multiply-complex-matrices vects 
;                                        (magicl:make-complex-matrix
;                                          4
;                                          1
;                                          (list 1 2 -1 2)))
;                                          ))
;      ;;(format t "M~%")
;      ;(print-matrix m)
;      ;(format t "Eigenvalues LAMBDA~%")
;      ;(print-matrix val-diag)
;      ;(format t "Eigenvectors V~%")
;      ;(print-matrix vects)
;      ;(format t "M*V~%")
;      ;(print-matrix (multiply-complex-matrices m vects))
;      ;(format t "V*LAMBDA~%")
;      ;(print-matrix (multiply-complex-matrices vects val-diag))
;      )


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

;(defun rotate-back (vec)
    ;(3d-matrices:m*
      ;(3d-matrices:mrotation (vec 0 0 1) )))

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

(defun v-moon (cur-time)
    (v* (vec 
              (* -1 *orb-angular-v-moon*
                 (sin (+ (* cur-time 
                         *orb-angular-v-moon*)
                         *moon-init-ang*)))
              (* *orb-angular-v-moon*
                 (cos (+ (* cur-time 
                         *orb-angular-v-moon*)
                         *moon-init-ang*)))
              0d0)
            *orb-rad-moon*))


(defparameter *d-moon-l2-approx*
    (* *orb-rad-moon*
       (expt (/ *mass-moon*
                (* 3 *mass-earth*))
             1/3)))


(defstruct (body 
            (:constructor %make-body)) 
    (pos (vec 0 0 0) :type vec3 )
    (v (vec 0 0 0) :type vec3 )
    (mass 0d0 :type long-float)
    (name "" :type string)
    (grav-significant T :type boolean))

(defun make-body (pos v mass 
                  &optional (name (write-to-string (random 999999)))
                            (grav-significant T))
    (setf mass (coerce mass 'double-float))
    (setf pos (vec (vx pos) (vy pos) (vz pos)))
    (setf v (vec (vx v) (vy v) (vz v)))
    (%make-body :pos pos :v v :mass mass :name name :grav-significant grav-significant)
    )

(defparameter *rocket* (make-body 
                         (vec -4.54038643d8
                              0
                              0)
                         (v-moon 0)
                         1000))


(defparameter *final-rocket* 
  (make-body 
    (v+ (pos-moon 0)
        (vec (* -1 *ldist*) 0 0)
        )
    (v* (v-moon 0) (/ (+ *orb-rad-moon* *ldist*) *orb-rad-moon*))
    1000))

(defparameter *L2-rocket* 
  (make-body 
    (v+ (pos-moon 0)
        (vec (* -1 *ldist*) 0 0)
        )
    (v* (v-moon 0) (/ (+ *orb-rad-moon* *ldist*) *orb-rad-moon*))
    1000))

(defparameter *phase* (- (*   *orb-angular-v-moon* 
                         (hm-transfer-period *mass-earth* 
                                             (+ *rad-earth* 200d3) 
                                             (+ *orb-rad-moon* *orb-rad-earth* ))) 
                         0.0785))

(defparameter *simple-to-moon*
  (make-body 
    (v+ (pos-earth 0)
        (v* (+ *rad-earth* 200d3)
            (vunit (vec (cos *phase*) (sin *phase*) 0)))
        )
    
    (v* (+ (orbital-v *mass-earth* (+ *rad-earth* 200d3)) 
           (hm-dv-1 *mass-earth* (+ *rad-earth* 200d3) (+ *orb-rad-moon* *orb-rad-earth* ))
           5)
        (vunit (vec (* -1 (sin *phase*)) (* 1 (cos *phase*)) 0))
        )
    1000))


(defun body-jacobi-em (body)
    ())



(defparameter *free-rocket* 
  (make-body 
    (vec 0 10d8 0)
    (vec 100 0 0)
    1000) )


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

(defun grav (r1 r2 m)
    ())

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
   ;accel-fun should take one argument, a number between 0 and dt.
   ;It's needed to do r4k by recalculating the forces applied at various advanced timesteps.
   (declare (optimize (safety 1) (speed 3)))
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

(defun energy-earth-moon (curtime)
    (+ (expt (vlength (v-earth curtime)) 2) 0))


(defun parameter-distance (bod1 bod2)
    (+
     (vdistance (body-pos bod1)
                (body-pos bod2))
     (vdistance (body-v bod1)
                (body-v bod2))))

(defun just-simulate (&optional (rock-orig *rocket*)
                                (start 0)
                                (end (* 3600 24 7))
                                (dt 20))
    (let* (
           (rocket (copy-structure rock-orig))
           (sim-steps (round (/ (- end start) dt))))
      
      (iter (for i from 0 to (1- sim-steps))
            (let* ((curtime (* dt i)))

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
                           (setf rocket (make-body (vec 0 0 0) (vec 0 0 0) 10))
                           (iter:finish)))))
      
      rocket))

(defun gravity-accel (this-pos)
    (let* ((to-earth (v- (pos-earth 0) this-pos))
           (to-moon (v- (pos-moon 0) this-pos)))
      (v+ (v* (vunit to-earth) 
              *big-g*
              *mass-earth* 
              (expt (vlength to-earth) -2))
          (v* (vunit to-moon) 
              *big-g*
              *mass-moon* 
              (expt (vlength to-moon) -2)))))

(defparameter *grav-from-e-at-m*
    (v* (vunit (v- (pos-earth 0) (pos-moon 0))) 
                *big-g*
                *mass-earth* 
                (expt (vlength (v- (pos-earth 0) (pos-moon 0)))
                      -2)))

(defun plot-grav (&optional (x-range   100000) 
                            (x-fidelity 5000)
                            (x-guess  -4.5679d8)
                            (y-range   100000) 
                            (y-fidelity 5000)
                            (y-guess 0))
    (let* (
           (x-vector (range (- x-guess (/ x-range 2)) 
                            (+ x-guess (/ x-range 2))
                              x-fidelity))
           (y-vector (range (- y-guess (/ y-range 2)) 
                            (+ y-guess (/ y-range 2))
                              y-fidelity))
           (nresults (* (length x-vector)
                        (length y-vector)))
           (x-results (make-array nresults))
           (y-results (make-array nresults))
           (grav-results (make-array nresults))
           (i 0)
           )
      (assert (= nresults (* (length x-vector) (length y-vector))))
      (iter (for x in-vector x-vector)
            (iter (for y in-vector y-vector)
                  (setf (aref x-results i) (abs (- x (vx (pos-moon 0)))))
                  (setf (aref y-results i) y)
                  (setf (aref grav-results i)
                        (abs (- (vlength (gravity-accel (vec x y 0)))
                           (vlength *grav-from-e-at-m*))))
                  (incf i)
                  ))

        (format t "~A~%" (let ((the-i 
                               (iterate (for i from 0 to (1- nresults))
                                 (finding i minimizing 
                                          (aref grav-results i)))))
                         (format nil "pos: ~A vel: ~A dist: ~A" 
                                 (aref x-results the-i)
                                 (aref y-results the-i)
                                 (aref grav-results the-i))))

      
      ;;time to draw
      (close-all-plots)
      (3d-plot x-results y-results grav-results)
      ;(format-plot nil "set hidden3d")
      (format-plot nil "set dgrid3d 50,50")
      (format-plot nil "set pm3d " )
      ;(format-plot nil "set lines")
      ;(format-plot nil "set contour base")
      ;(format-plot nil "set view square xy")
      (format-plot nil "set tics font \", 10\"")
      (format-plot nil "set xlabel font \", 15\"")
      (format-plot nil "set ylabel font \", 15\"")
      (format-plot nil "set title font \", 18\"")
      (xlabel "Distance away from moon (m)")
      (ylabel "Distance perpendicular to earth-moon axis (m)")
      (title "Difference in gravitational pull relative to lunar orbit about estimated location of L2 point (m/s^2)")
      ))

(defun sim-planets-and-rocket (&optional (rock-orig *rocket*)
                                         (start 0) 
                                         (end (* 3600 24 7)) 
                                         (dt 20) 
                                         (n-plot-samples 1000))
    ;(declare (optimize (safety 0) (speed 3)) )
    (close-all-plots)
    (let* (
           (not-yet-encountered-moon t)
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
           (d-r-e (make-array n-plot-samples))
           )
      
      (iter (for i from 0 to (1- sim-steps))
            (let* ((curtime (* dt i)))

             ;this is all just plotting code
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
                                                m-v))
                (setf (aref d-r-e j) (vdistance (body-pos rocket)
                                                e-v))
                ))


              ;the actual core simulation code
              (iterate-body-r4k rocket 
                                (coerce dt 'double-float) 
                                (simple-acceleration curtime))
              ;basically the scripting area
              (when (and (< (vdistance (body-pos rocket)
                                      (pos-moon curtime))
                           (+ *rad-moon* 100d3))
                         not-yet-encountered-moon)
                  (setf not-yet-encountered-moon nil)
                  (format t "close!~%")
                  (let* ((lunar-orbital-v (orbital-v *mass-moon*
                                                     (vlength (v- (body-pos rocket)
                                                                  (pos-moon curtime)))))
                         (target-v (v+ (v-moon curtime)
                                       (v* 
                                         (vc (vec 0 0 1)
                                             (vunit 
                                               (v- (body-pos rocket)
                                                   (pos-moon curtime))))
                                         lunar-orbital-v))
                                   )
                         
                         (delta-v (vlength (v- (body-v rocket) target-v)))
                         (dv-r-m (v- (body-v rocket) (v-moon curtime)))
                         
                         )
                   (format t "Needed orbital velocity relative to moon: ~a~%"
                             lunar-orbital-v)
                   (format t "Current velocity:  ~a, ~a, ~a~%"
                           (vx (body-v rocket))(vy (body-v rocket))(vz (body-v rocket)))
                   (format t "Target velocity for lunar orbit:  ~a, ~a, ~a~%"
                           (vx target-v) (vy target-v) (vz target-v))
                   (format t "Delta-v needed to reach target vector: ~a~%"
                             delta-v)
                   (setf (body-v rocket ) 
                         target-v) )
                  )
             
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
        ;(format-plot nil "set zlabel \"Time elapsed (s)\"")
       
        
        
        (3d-plot e-x e-y timevec "earth" 
                 m-x m-y timevec "moon" 
                 r-x r-y timevec "k;rocket" )
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
         
         (format-plot nil "set key left bottom")
         
         (format-plot nil "set title \"Hohmann-approximating transfer orbit with lunar orbital insertion\"")
         (format-plot nil "set title bottom")
         (format-plot nil "set view equal xy") 
         
         (format-plot nil "set size ratio -1")
         (format-plot nil "set size square")
         )
       ;(read )
      ))

;(when nil (progn
        ;(vgplot:close-all-plots)
        ;;(setf timevec (map 'vector 
                           ;(lambda (x) (/ x (/ (* 2 pi) *orb-angular-v-moon*)))
                           ;timevec))
        ;(format-plot nil "set xlabel \"Time elapsed (lunar periods)\"")
        ;(format-plot nil "set ylabel \"Rocket-Moon distance (m)\"")
      

        ;(text 0 105000000 "Distance Limit")
        ;(text 0  79000000 "Lagrange point guess")
        ;(text 0   5000000 "Radius of Moon")
        ;(let ((dtv (map 'vector
                        ;(lambda (x)
                                ;(/ x *orb-period-moon*))
                        ;timevec)))
          ;(plot dtv d-r-m 
                ;dtv (make-array (list (length dtv)) 
                                ;:initial-element 100000d3)
                ;dtv (make-array (list (length dtv)) 
                                ;:initial-element 7.43103725d7)
                ;dtv (make-array (list (length dtv)) 
                                ;:initial-element *rad-moon*)))))


(defun pcr3brun ()
    (sim-planets-and-rocket *simple-to-moon* 0 (days 4.4) 5 3000 ))

(defun lagrange-over-range (&optional (pos-range 100) 
                                      (pos-fidelity 1)
                                      (pos-guess (vec -4.54041001d8
                                                      0
                                                      0))
                                      (v-range 0.01) 
                                      (v-fidelity 0.0001)
                                      (v-guess (vec 0 
                                                    -1030.828735
                                                    0)))
    (let* (
           (posmag-guess (vlength pos-guess))
           (pos-vector (range (- posmag-guess (/ pos-range 2)) 
                              (+ posmag-guess (/ pos-range 2))
                              pos-fidelity))
           (pos-guess (vunit pos-guess))
           (vmag-guess (vlength v-guess))
           (v-guess (vunit v-guess))
           (v-vector   (range (- vmag-guess (/ v-range 2)) 
                              (+ vmag-guess (/ v-range 2))
                              v-fidelity))
           (nresults (* (length v-vector)
                        (length pos-vector)))
           (pos-results (make-array nresults))
           (v-results (make-array nresults))
           (dist-results (make-array nresults))
           (i 0)
           (sim-iterations (round *orb-period-moon*))
           )
      (assert (= nresults (* (length pos-vector) (length v-vector))))
      

      (iter (for p in-vector pos-vector)
            (iter (for v in-vector v-vector)
              (setf (aref pos-results i) p)     
              (setf (aref v-results i) v)
              (setf (aref dist-results i)
                    (parameter-distance 
                      (make-body (v* pos-guess
                                     p)
                                 (v* v-guess v)
                                 1000)
                      (just-simulate 
                        (make-body (v* pos-guess
                                       p)
                                   (v* v-guess v)
                                   100)
                        0
                        sim-iterations
                        100)))
              (incf i)))
      (format t "~A~%" (let ((the-i 
                               (iterate (for i from 0 to (1- nresults))
                                 (finding i minimizing 
                                          (aref dist-results i)))))
                         (format nil "pos: ~A vel: ~A dist: ~A" 
                                 (aref pos-results the-i)
                                 (aref v-results the-i)
                                 (aref dist-results the-i))))
      ;;time to draw
      (close-all-plots)
      (3d-plot pos-results v-results dist-results)
      (format-plot nil "set hidden3d")
      (format-plot nil "set dgrid3d")
      (format-plot nil "set pm3d")
      (format-plot nil "set view square xy")
      (format-plot nil "set tics font \", 10\"")
      (format-plot nil "set xlabel font \", 15\"")
      (format-plot nil "set ylabel font \", 15\"")
      (xlabel "Initial position")
      (ylabel "Initial velocity")
      )) 


(defun lhs (r)
    (let ((big-r *d-earth-moon*)
          (M1 *mass-earth*)
          (M2 *mass-moon*)
          )
      (+ (/ M1 
           (expt (+ big-r r) 2))
         (/ M2
            (expt r 2)))
         ))

(defun d-lhs (r)
    (let ((big-r *d-earth-moon*)
          (M1 *mass-earth*)
          (M2 *mass-moon*)
          )
      (+ (/ (* -2 M1) 
           (expt (+ big-r r) 3))
         (/ (* M2 -2)
            (expt r 3)))
         ))


(defun rhs (r)
    (let ((big-r *d-earth-moon*)
          (M1 *mass-earth*)
          (M2 *mass-moon*)
          )
      (+ (/ M1
            (expt big-r 2))
         (/ (* r
               (+ M1 M2))
            (expt big-r 3)))))
(defun d-rhs (r)
    (let ((big-r *d-earth-moon*)
          (M1 *mass-earth*)
          (M2 *mass-moon*)
          )
      
      (/ (* 1
            (+ M1 M2))
         (expt big-r 3))))

(defun diff (r)
  (- (lhs r) (rhs r)))

(defun d-diff (r)
  (- (d-lhs r) (d-rhs r)))

(defun nr-next (func deriv guess)
    (- guess
       (/ (funcall func guess)
          (funcall deriv guess))))

(defun nr-approx (func deriv guess)
    (let ((next-guess (nr-next func deriv guess)))
        (do ()
            ((< (abs (/ (- guess next-guess) 
                        next-guess))
                0.01))
            (setf guess next-guess)
            (setf next-guess (nr-next func deriv guess)))
        next-guess))


;(defparameter *a-e-m*
    ;(a (mu *mass-earth* *mass-moon*)
       ;0.5 ))
;

;(defparameter *lagrangedist* 
  ;(nr-approx #'diff #'d-diff 1))
;

(defun tot (x)
    (v* (+ *orb-rad-moon* *orb-rad-earth*)
        x))

(defun pos-l2 (curtime)
    (v+ (pos-moon curtime)
        (v* (vunit (pos-moon curtime))
            *ldist*)))

(defun v-l2 (curtime)
  (v* (v-moon curtime) (/ (+ *orb-rad-moon* *ldist*) *orb-rad-moon*)))


(defparameter *last-rocket* 
  (make-body  (v+ (pos-l2 0)  (tot (vec 1.326  0 0)))
              (rotating-to-vel 
                (v+ (pos-l2 0)  (tot (vec 1.326  0 0))) 
                (vec 0 4.301 0)
                *orb-angular-v-moon* )
                
              
              1000
  
  ))

(defparameter *last-rocket-2* 
  (make-body  (v+ (pos-l2 0) (vec 10 0 0))
              (v+ (v-l2 0) (vec 10 0 0))
              1000
  ))

    ;(declare (optimize (safety 0) (speed 3)) )
;300 280 240 200 160 80 40 20 10
(defun last-sim (&optional (rock-orig *rocket*)
                                         (start 0) 
                                         (end (* 3600 24 7)) 
                                         (n-plot-samples 1000))
    (close-all-plots)
    (let ((rockets nil)
          (timevecs nil)
          (dtlist (list 50 -50 )))
     (dolist (this-dt dtlist)
      (let* (
           (dt this-dt) 
           (not-yet-encountered-moon t)
           (not-yet-boost t)
           (not-yet-encountered-earth t)
           (rocket (copy-structure rock-orig))
           (timevec (range start end (/ (- end start) n-plot-samples)))
           (sim-steps (round (/ (- end start) (abs dt))))
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
           (d-r-e (make-array n-plot-samples))
           )
      
      (iter (for i from 0 to (1- sim-steps))
            (let* ((curtime (* dt i)))

             ;this is all just plotting code
             (when (= i (caar plot-coords))
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
                                                m-v))
                (setf (aref d-r-e j) (vdistance (body-pos rocket)
                                                e-v))
                ))

              ;the actual core simulation code
              (iterate-body-r4k rocket 
                                (coerce dt 'double-float) 
                                (simple-acceleration curtime))
              ;basically the scripting area
              
              (when (and (< (vdistance (body-pos rocket)
                                      (pos-moon curtime))
                           (+ *rad-moon* 7500d3))
                         not-yet-encountered-moon)
                  (setf not-yet-encountered-moon nil)
                  (format t "close!~%")
                  (let* ((lunar-orbital-v (orbital-v *mass-moon*
                                                     (vlength (v- (body-pos rocket)
                                                                  (pos-moon curtime)))))
                         (target-v (v+ (v-moon curtime)
                                       (v* 
                                         (vc (vec 0 0 1)
                                             (vunit 
                                               (v- (body-pos rocket)
                                                   (pos-moon curtime))))
                                         lunar-orbital-v))
                                   )
                         
                         (delta-v (vlength (v- (body-v rocket) target-v)))
                         (dv-r-m (v- (body-v rocket) (v-moon curtime)))
                         (to-moon (vlength (v- (pos-moon curtime) (body-pos rocket))))
                         )
                   (format t "Needed orbital velocity relative to moon: ~a~%"
                             lunar-orbital-v)
                   (format t "Current velocity:  ~a, ~a, ~a~%"
                           (vx (body-v rocket))(vy (body-v rocket))(vz (body-v rocket)))
                   (format t "Target velocity for lunar orbit:  ~a, ~a, ~a~%"
                           (vx target-v) (vy target-v) (vz target-v))
                   (format t "Delta-v needed to reach target vector: ~a~%"
                             delta-v)
                   (format t "Delta-v needed to get down to moon:~a~%"
                             (+ (hm-dv-1 *mass-moon* (+ *rad-moon* 100d3) to-moon) (hm-dv-2 *mass-moon*  (+ *rad-moon* 100d3) to-moon)))
                   (setf (body-v rocket ) 
                         target-v) )
                  )
              (when (and  not-yet-boost
                          (> (vdistance (body-pos rocket)
                                      (pos-earth curtime))
                          (* 3 *orb-rad-moon* ))
                         (< (abs (v. (vunit (v- (body-pos rocket)
                                                (pos-earth curtime)))
                                     (vunit (body-v rocket))))
                            0.01))
                (setf not-yet-boost nil)
                (format t "Current v: ~A~%" (vlength (body-v rocket)))
                (format t "Target v: ~A~%" (hm-v-2 *mass-earth* 
                                  (+ *rad-earth* 200d3)
                                  (vlength (v- (body-pos rocket)
                                               (pos-earth curtime)))))
                (format t "Current v: ~A~%" (vlength (v- (body-v rocket)
                                                         )))
                (setf (body-v rocket)
                      (v* (vc (vunit (v- (vec 0 0 0)
                                         (body-pos rocket )
                                     ))
                              (vec 0 0 1))
                         (hm-v-2 (+ *mass-moon* *mass-earth*) 
                                  (+ *rad-earth* 3000d3)
                                  (vlength (v- (body-pos rocket)
                                               (pos-earth curtime)))))))
             
              (when (and not-yet-encountered-earth 
                         (< (vdistance (body-pos rocket)
                                      (pos-earth curtime))
                           (+ *rad-earth* 100000d3))
                         (< (abs (v. (vunit (v- (body-pos rocket)
                                                (pos-earth curtime)))
                                     (vunit (body-v rocket))))
                            0.005))
                (setf not-yet-encountered-earth nil)
                (format t "current delta-v ~A~%"
                        (vlength
                         (v-   (v* (vc (vunit (v- (pos-earth curtime)
                                         (body-pos rocket )
                                     ))
                              (vec 0 0 1))
                         (orbital-v  *mass-earth* (vdistance (body-pos rocket)
                                      (pos-earth curtime))))

                             (body-v rocket))))
               (format t "Delta-v needed to get down to LEO~A~%"
                             (+ (hm-dv-1 *mass-earth* (+ *rad-earth* 200d3) (vlength (v- (body-pos rocket) (pos-earth curtime)) )) 
                                (hm-dv-2 *mass-earth* (+ *rad-earth* 200d3) (vlength (v- (body-pos rocket) (pos-earth curtime))))))
                (setf (body-v rocket)
                      (v* (vc (vunit (v- (pos-earth curtime)
                                         (body-pos rocket )
                                     ))
                              (vec 0 0 1))
                         (orbital-v  *mass-earth* (vdistance (body-pos rocket)
                                      (pos-earth curtime)))))

               
                )

              (when (or (< (vdistance (body-pos rocket)
                                      (pos-moon curtime))
                           *rad-moon*)
                        (< (vdistance (body-pos rocket)
                                      (pos-earth curtime))
                           *rad-earth*))
                    (progn (format t "COLLISION!!!~%")
                           (iter:finish)))))
      
    (if (or (eq dt -50) (eq dt -5))  
        (progn 

        (push r-y rockets)
        (push r-x rockets)
        (push timevec timevecs)
        (vgplot:close-all-plots)
        
        (3d-plot 
                 e-x e-y timevec "earth" 
                 m-x m-y timevec "moon" 
                 (pop rockets ) (pop rockets) (pop timevecs) "LEO to L2" 
                 ;(subseq (pop rockets ) 0 3000) (subseq (pop rockets) 0 3000) (subseq (pop timevecs) 0 800)  "L2 to Lunar Ballistic Capture" 
                 
                 )
         
         (format-plot nil "set key left bottom")
         (format-plot nil "set key font \", 15\"")
         (format-plot nil "set title font \", 15\"")
         (format-plot nil "set tics font \", 8\"")
         (format-plot nil "set title \"\"")
         (format-plot nil "set title bottom")
         (format-plot nil "set view equal xy") 
         
         (format-plot nil "set size ratio -1")
         (format-plot nil "set size square")
         )
        )
    (progn  
      
        (push r-y rockets)
        (push r-x rockets)
        (push timevec timevecs)
)))))


 
