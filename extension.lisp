(in-package #:comp-phys)

(defstruct (swpvec
             (:constructor %make-swpvec)
             (:copier %copy-swpvec))
    (current (vector) :type simple-vector) 
    (next (vector) :type simple-vector))

(defun make-swpvec (&optional (vec1 (vector) )
                              (vec2 (vector) vec2-p))
    (when (not vec2-p)
          (setf vec2 vec1))
    (assert (equal (length vec1)
                   (length vec2)))
    (%make-swpvec :current vec1 :next vec2))

(defun copy-array (array &key
                   (element-type (array-element-type array))
                   (fill-pointer (and (array-has-fill-pointer-p array)
                                      (fill-pointer array)))
                   (adjustable (adjustable-array-p array)))
    "Returns an undisplaced copy of ARRAY, with same fill-pointer and
    adjustability (if any) as the original, unless overridden by the keyword
    arguments."
    (let* ((dimensions (array-dimensions array))
           (new-array (make-array dimensions
                                  :element-type element-type
                                  :adjustable adjustable
                                  :fill-pointer fill-pointer)))
      (dotimes (i (array-total-size array))
        (setf (row-major-aref new-array i)
              (row-major-aref array i)))
      new-array))

(defun equarray (arr1 arr2)
    (assert (arrayp arr1))
    (assert (arrayp arr2))
    (if (eql (length arr1)
             (length arr2))
        (dotimes (i (length arr1) T)
            (when (not (equal (aref arr1 i)
                              (aref arr2 i)))
                  (return nil)))
        nil))

(defun eqswpvec (obj1 obj2)
    (and (equarray (swpvec-current obj1)
                   (swpvec-current obj2))
         (equarray (swpvec-next obj1)
                   (swpvec-next obj2))))

(defun copy-swpvec (obj)
    (make-swpvec (copy-array (swpvec-current obj))
                 (copy-array (swpvec-next obj))))

(defun copy-celestial (vec)
    (map-into (make-array (length vec))
                           #'copy-body
                           vec))

(defun copy-celestial-swpvec (vec)
    (make-swpvec (map-into (make-array (length (swpvec-current vec)))
                           #'copy-body
                           (swpvec-current vec))
                 (map-into (make-array (length (swpvec-next vec)))
                           #'copy-body
                           (swpvec-next vec))))

(defun test-swaps ()
  (let* ((a (make-swpvec (vector 1 2 3 4 5 6 7 8 9) (vector 11 12 13 14 15 16 17 18 19))) 
         (b (copy-swpvec a))) 
    (print a)
    (print b)
    (print (eqswpvec a b))
    (setf (aref (swpvec-current b) 0) 9999) 
    (print a) 
    (print b)
    (print (eqswpvec a b))
    "And that's it"
    ))

(defun swpvec-swap (vec)
    (let ((temp (swpvec-current vec))) 
      (setf (swpvec-current vec)
            (swpvec-next vec))
      (setf (swpvec-next vec)
            temp)))


(defun shift-to-zero-net-momentum (celestvec)
    (let* (
          (total-mass (reduce 
                        (lambda (bod1 bod2)
                          (+ (if (eq (type-of bod1) 
                                     'body) 
                                 (body-mass bod1)
                                 bod1)
                             (body-mass bod2)))
                        celestvec))
          (net-x (/ (reduce 
                      (lambda (bod1 bod2) 
                        (+ (if (body-p bod1) 
                             (* (body-mass bod1) 
                                (vx (body-v bod1)))
                              bod1) 
                           (* (body-mass bod2) 
                              (vx (body-v bod2))))) 
                      celestvec)
                    total-mass))
          (net-y (/ (reduce 
                      (lambda (bod1 bod2) 
                        (+ (if (body-p bod1) 
                             (* (body-mass bod1) 
                                (vy (body-v bod1)))
                              bod1) 
                           (* (body-mass bod2) 
                              (vy (body-v bod2))))) 
                      celestvec)
                    total-mass))
          (net-z (/ (reduce 
                      (lambda (bod1 bod2) 
                        (+ (if (body-p bod1) 
                             (* (body-mass bod1) 
                                (vz (body-v bod1)))
                              bod1) 
                           (* (body-mass bod2) 
                              (vz (body-v bod2))))) 
                      celestvec)
                    total-mass))
          (copy (copy-celestial celestvec)))
          
          (map-into copy
                    (lambda (bod)
                      (setf (vx (body-v bod))
                            (- (vx (body-v bod))
                               net-x))
                      (setf (vy (body-v bod))
                            (- (vy (body-v bod))
                               net-y))
                      (setf (vz (body-v bod))
                            (- (vz (body-v bod))
                               net-z))
                      bod)
                    copy)
          copy))

(defparameter *celestials-initial*
    (vector 
      (make-body (VEC3 4671759.51188055d0 0.0d0 0.0d0)
                 (VEC3 -0.0d0 9293.99482022799d0 0.0d0)
                 *mass-earth*
                 "Earth")
      (make-body (VEC3 -3.797282404881194d8 4.650329742596318d-8 0.0d0)
                 (VEC3 -1.2624583775320132d-13 -1030.8754968460996d0 0)
                 *mass-moon*
                 "Moon")))
(defparameter *celes-sim-init* 
    (make-swpvec *celestials-initial*
                 *celestials-initial*))

(defparameter *earth* 
  (make-body (VEC3 *orb-rad-earth* 0.0d0 0.0d0)
                   (VEC3 0 0 0)
                   *mass-earth*
                   "Earth"))

(defparameter *mass-sun* 1.9885d30)
(defparameter *mass-jupiter* 1.8982d27)

(defparameter *dead-inital*
    (vector 
        *earth*
        (make-body (VEC3 (* -1 *orb-rad-moon*) 0.0d0 0.0d0)
                   (VEC3 0 500 0)
                   *mass-moon*
                   "Moon"
                   )
       ; (make-body (VEC3 149.6d9 0.0d0 0.0d0)
       ;            (VEC3 0 30000 0)
       ;            *mass-sun*
       ;            "Sun")
       ; (make-body (vec3 (+ 149.6d9 778.57d9) 0 0)
       ;            (vec3 0 0 0)
       ;            *mass-jupiter*
       ;            "Jupiter")
        (make-body
                   (VEC3 (* -1 (+ *orb-rad-moon* *lagrangedist*)) 0.0d0 0.0d0)
                   (vec 0 (* 1022
                             (/ (+ *orb-rad-moon* *lagrangedist*)
                                *orb-rad-moon*))  
                        0)
                   0
                   "Small rock"
                   nil)))

(defparameter *dead-sim-init*
    (make-swpvec *dead-inital*
                 *dead-inital*))

(defun grav-acceleration (from-pos to-bod)
    (declare (optimize (speed 3) (safety 0)))
    (declare (type vec3 from-pos)
             (type body to-bod)
             (type double-float *big-g* ))
    ;(format t "this is a test~%")
    (let* ((dirvec (v- (body-pos to-bod)
                       from-pos))
          (dir (vunit dirvec))
          (magsqr (expt (vlength dirvec) 2)))
      ;(format t "dir is ~A~%" dir)
      (v* dir
          (/ (* *big-g* 
               (body-mass to-bod))
            magsqr))))

(defun total-grav-acceleration (pos all-bodies 
                          &optional (excluded-indices (vector)))
    (declare (optimize (speed 3) (safety 0))
             (type simple-vector excluded-indices)
             (type simple-vector all-bodies)
             (inline grav-acceleration))
    (let ((total (vec 0d0 0d0 0d0)))
          (dotimes (i (length all-bodies))
                   (when (and (not (find i excluded-indices))
                              (body-grav-significant 
                                (aref all-bodies i))) 
                         (nv+ total
                              (grav-acceleration pos
                                                 (aref all-bodies
                                                       i)))))
          total))

(defun iterate-n-body-r4k (body-old body-new all-bodies dt 
                                &optional (excluded-indices (vector))
                                          (other-accel (vec 0 0 0)))
    (%iterate-n-body-r4k body-old body-new all-bodies (round dt 1) 
                          excluded-indices other-accel))

(defun %iterate-n-body-r4k (body-old body-new all-bodies dt 
                            &optional (excluded-indices (vector))
                                      (other-accel (vec 0 0 0)))
    (declare (optimize (safety 0) (speed 3))
             (type integer dt))
    (let* ((sixth-dt (/ (round dt 1) 6))
           (half-dt (/ (round dt 1) 2))
           (orig-v (body-v body-old))
           (orig-a (v+ 
                     (total-grav-acceleration (body-pos body-old)
                                              all-bodies
                                              excluded-indices)
                     other-accel))

           (z1 (v+ (body-pos body-old) 
                    (v* (body-v body-old) 
                         half-dt)))
           (z1-dot (v+ (body-v body-old) 
                        (v* orig-a
                             half-dt)))
           (z1-a (v+ 
                   (total-grav-acceleration z1
                                            all-bodies
                                            excluded-indices)
                   other-accel))
           (z2 (v+ (body-pos body-old) 
                   (v* z1-dot 
                       half-dt)))
           (z2-dot (v+ (body-v body-old) 
                       (v* z1-a
                           half-dt)))
           (z2-a (v+ (total-grav-acceleration z2
                                              all-bodies
                                              excluded-indices)
                     other-accel))
           (z3 (v+ (body-pos body-old) 
                   (v* z2-dot 
                        dt)))
           (z3-dot (v+ (body-v body-old) 
                       (v* z2-a
                            dt)))
           (z3-a (v+ (total-grav-acceleration z3
                                              all-bodies
                                              excluded-indices)
                     other-accel))
           )
      (setf (body-v body-new)
            (v+ (body-v body-old)
               (v* sixth-dt
                  (v+ orig-a
                     (v* 2 z1-a)
                     (v* 2 z2-a)
                     z3-a))))
      (setf (body-pos body-new)
            (v+ (body-pos body-old)
               (v* sixth-dt
                  (v+ orig-v
                     (v* 2 z1-dot)
                     (v* 2 z2-dot)
                     z3-dot))))))

(defun iterate-n-bodies-r4k (body-swpvec dt)
    (declare (optimize (speed 3) (safety 0)))
    (dotimes (i (length (swpvec-current body-swpvec)))
        (iterate-n-body-r4k (aref (swpvec-current body-swpvec) i) 
                            (aref (swpvec-next  body-swpvec) i)
                            (swpvec-current body-swpvec)
                            dt
                            (vector i)))
    (swpvec-swap body-swpvec))

(defmacro simulating-n-bodies-while-doing ((start-state &optional (duration 3600) (dt 20) (?debug? nil)) &body body)
    `(let* ((sim-steps (round (/ ,duration ,dt)))
           (bodies (copy-celestial-swpvec ,start-state)))
      (progn 
        (dotimes (i sim-steps)
          ,@body
          (iterate-n-bodies-r4k bodies ,dt)
          (when ,?debug? (format t "~A~%" (aref (swpvec-current bodies) 0)))))
      (assert (not (eqswpvec ,start-state bodies)))))

(defun sim-while-plotting (initial-state &key (duration (days 1)) (dt 20) (bodies-to-plot '() bodies-to-plot-p ) (n-plot-points 1000))
    (let* ((copied-state (copy-celestial-swpvec initial-state))
           (bodies-to-plot (if bodies-to-plot-p
                               (coerce bodies-to-plot 'vector)
                               (range 
                                 (length 
                                   (swpvec-current copied-state)))))
           (sim-steps (round (/ duration dt)))
           (plot-data (make-array (list (* 3 (length (swpvec-current copied-state)))
                                         n-plot-points)))
           (timevec (range 0 duration (/ duration n-plot-points)))
           (plot-ratio (/ sim-steps n-plot-points))
           (plot-coords
             (nreverse
               (cons `(,(1+ sim-steps) ,(1+ n-plot-points))
                 (nreverse (iter 
                             (for i from 0 to (1- n-plot-points))
                             (collecting (list 
                                           (min 
                                             (round (* i plot-ratio))
                                             sim-steps)
                                           i))))))))
      (format t "~A~%" plot-ratio)
      (print (swpvec-current copied-state))
      
      
      (simulating-n-bodies-while-doing (copied-state duration dt)
          (when (= i (caar plot-coords))
            (let* ((sim-index i)
                   (plot-index (cadr (pop plot-coords))))
              ;(format t "~A~%" (body-pos (aref (swpvec-current bodies) 0)))
              (iter (for body-index in-vector bodies-to-plot)
                    (setf (aref plot-data (+ 0 (* 3 body-index)) plot-index ) 
                          (vx (body-pos 
                                (aref (swpvec-current bodies)
                                      body-index))))
                    (setf (aref plot-data (+ 1 (* 3 body-index)) plot-index ) 
                          (vy (body-pos 
                                (aref (swpvec-current bodies)
                                      body-index))))
                    (setf (aref plot-data (+ 2 (* 3 body-index)) plot-index ) 
                          (aref timevec plot-index))))))
      (close-all-plots)
      (format-plot nil "set zlabel \"Time elapsed (s)\"")
      (format-plot nil "set view equal xy")
      ;(format-plot nil "set size ratio -1")
      (format-plot nil "set size square")
      (eval `(3d-plot 
              ,@(iter (for body-index in-vector bodies-to-plot)
                      (collect (make-array n-plot-points
                                           :displaced-to plot-data 
                                           :displaced-index-offset (* (+ 0 (* 3 body-index))
                                                                      n-plot-points)))
                      (collect (make-array n-plot-points
                                           :displaced-to plot-data 
                                           :displaced-index-offset (* (+ 1 (* 3 body-index))
                                                                      n-plot-points)))
                      ;(collect (map (type-of timevec) #'cos timevec))
                      ;(collect (map (type-of timevec) #'sin timevec))
                      (collect timevec)
                      (collect (body-name (aref (swpvec-current copied-state) body-index))))))
      
      ))

(defun run ()
    (sim-while-plotting *dead-sim-init* 
                        :duration (days 30) 
                        :n-plot-points 1000
                        :bodies-to-plot '(0 1 )

                        ))


;--------------------------------------------------------

(defparameter *e* 2.718281828459045235360287471352d0)

(defun veloc (pos E)
    (* (sqrt 2) (sqrt (+ E (/ 1 (* pos pos))))))

(defun doit () 
  (close-all-plots)
  (let* ((e -2)
         (upto (if (>= e 0) 
                        40
                        (sqrt (/ 1 (* -1 e)))))
         (veccy (range 0.1 upto 0.001))) 
    (plot 
    veccy
    (map 'vector 
         (lambda (x) (veloc x e))
         veccy))
       
    (format-plot nil "set tics font \", 10\"")
    (format-plot nil "set title font \", 18\"")
    (format-plot nil "set xlabel font \", 18\"")
    (format-plot nil "set ylabel font \", 18\"")
    (format-plot nil "set title \"Velocity magnitude vs. position for a gravitationally bound particle, E<0\"")
    (format-plot nil "set xlabel \"Distance from attracting mass\"")
    (format-plot nil "set ylabel \"Magnitude of velocity\"")
      )

  )

(defun doit-2 () 
  (close-all-plots)
  (let* ((e -2)
         (upto (if (>= e 0) 
                        40
                        (sqrt (/ 1 (* -1 e)))))
         (veccy (range 0.1 upto 0.001))) 
    (plot 
     (vector 0) (vector 0)
    )
       
    (format-plot nil "set tics font \", 10\"")
    (format-plot nil "set title font \", 18\"")
    (format-plot nil "set xlabel font \", 18\"")
    (format-plot nil "set ylabel font \", 18\"")
    (format-plot nil "set title \"Allowed positions for a particle of E<0\"")
    (format-plot nil "set xlabel \"x\"")
    (format-plot nil "set ylabel \"y\"")
     
      
      )

  )

(defun doit-3 () 
  (close-all-plots)
  (let* ((e -2)
         (upto (if (>= e 0) 
                        40
                        (sqrt (/ 1 (* -1 e)))))
         (veccy (range 0.1 upto 0.001))) 
    (plot 
     (vector 0) (vector 0)
    )
       
    (format-plot nil "set tics font \", 10\"")
    (format-plot nil "set title font \", 18\"")
    (format-plot nil "set xlabel font \", 18\"")
    (format-plot nil "set ylabel font \", 18\"")
    (format-plot nil "set title \"Allowed positions for a particle\"")
    (format-plot nil "set xlabel \"x\"")
    (format-plot nil "set ylabel \"y\"")
     
      
      )

  )

(defun jacobi (x y u1 u2 n)
    (+ 
      (* (* n n)
         (+ (* x x)
            (* y y)))
      (* 2
         (+ (/ u1 (sqrt (+ (expt (+ x u1) 2)
                           (expt y 2))))
            (/ u2 (sqrt (+ (expt (- x u2) 2)
                           (expt y 2))))))))

(defun jacobi-go (minval)
    (let* ((u1 (- 1 1/80))
           (u2 1/80)
           (xes (range -3 3 0.05))
           (yes (range -3 3 0.05))
           (points nil))
      (iter (for i from 0 to (1- (length xes)))
            (format t "~%")
            (iter (for j from 0 to (1- (length yes)))
               (format t "~A "
                       (if 
                        (and (>  (jacobi (aref xes j)
                               (aref yes i)
                               u1
                               u2
                               (/ (* 2 pi)
                                  (days 30)))
                            3.1)
                             (<  (jacobi (aref xes j)
                               (aref yes i)
                               u1
                               u2
                               (/ (* 2 pi)
                                  (days 30)))
                            3.3))
                        0
                        (jacobi (aref xes i)
                               (aref yes j)
                               u1
                               u2
                               (/ (* 2 pi)
                                  (days 30)))))
                      
               
               ))))

