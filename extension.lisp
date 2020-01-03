(in-package #:comp-phys)

(defstruct (swpvec
             (:constructor %make-swpvec))
    (current (vector) :type simple-vector) 
    (next (vector) :type simple-vector))

(defun make-swpvec (&optional (vec1 (vector) )
                              (vec2 (vector) vec2-p))
    (when (not vec2-p)
          (setf vec2 vec1))
    (assert (equal (length vec1)
                   (length vec2)))
    (%make-swpvec :current vec1 :next vec2))

(defun swpvec-swap (vec)
    (let ((temp (swpvec-current vec))) 
      (setf (swpvec-current vec)
            (swpvec-next vec))
      (setf (swpvec-next vec)
            temp)))

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

(defun grav-acceleration (from-pos to-bod)
    (declare (optimize (speed 3) (safety 0)))
    (declare (type vec3 from-pos)
             (type body to-bod)
             (type double-float *big-g* ))
    (let* ((dirvec (v- (body-pos to-bod)
                       from-pos))
          (dir (vunit dirvec))
          (magsqr (expt (vlength dirvec) 2)))
      (v* dir
          (/ (* *big-g* 
               (body-mass to-bod))
            magsqr))))

(defun total-grav-acceleration (pos all-bodies 
                          &optional (excluded-index -1))
    (declare (optimize (speed 3) (safety 0))
             (type integer excluded-index)
             (type simple-vector all-bodies)
             (inline grav-acceleration))
    (let ((total (vec 0d0 0d0 0d0)))
          (dotimes (i (length all-bodies))
                   (when (not (equal excluded-index i)) 
                         (nv+ total
                              (grav-acceleration pos
                                                 (aref all-bodies
                                                       i)))))
          total))

(defun iterate-n-body-r4k (body-old body-new all-bodies dt 
                                &optional (excluded-index -1)
                                          (other-accel (vec 0 0 0)))
    (%iterate-n-body-r4k body-old body-new all-bodies (round dt 1) 
                          excluded-index other-accel))

(defun %iterate-n-body-r4k (body-old body-new all-bodies dt 
                            &optional (excluded-index -1)
                                      (other-accel (vec 0 0 0)))
    (declare (optimize (safety 1) (speed 3))
             (type integer dt))
    (let* ((sixth-dt (/ (round dt 1) 6))
           (half-dt (/ (round dt 1) 2))
           (orig-v (body-v body-old))
           (orig-a (v+ 
                     (total-grav-acceleration (body-pos body-old)
                                              all-bodies
                                              excluded-index)
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
                                            excluded-index)
                   other-accel))
           (z2 (v+ (body-pos body-old) 
                   (v* z1-dot 
                       half-dt)))
           (z2-dot (v+ (body-v body-old) 
                       (v* z1-a
                           half-dt)))
           (z2-a (v+ (total-grav-acceleration z2
                                              all-bodies
                                              excluded-index)
                     other-accel))
           (z3 (v+ (body-pos body-old) 
                   (v* z2-dot 
                        dt)))
           (z3-dot (v+ (body-v body-old) 
                       (v* z2-a
                            dt)))
           (z3-a (v+ (total-grav-acceleration z3
                                              all-bodies
                                              excluded-index)
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
    (dotimes (i (length (swpvec-current body-swpvec)))
        (iterate-n-body-r4k (aref (swpvec-current body-swpvec) i) 
                            (aref (swpvec-next  body-swpvec) i)
                            (swpvec-current body-swpvec)
                            dt
                            i))
    (swpvec-swap body-swpvec))

(defmacro simulating-n-bodies-while-doing ((start-state &optional (duration 3600) (dt 20) (?debug? nil)) &body body)
    `(let* ((sim-steps (round (/ ,duration ,dt)))
           (bodies (copy-swpvec ,start-state)))
      (progn 
        (dotimes (i sim-steps)
          ,@body
          (iterate-n-bodies-r4k bodies ,dt)
          (when ,?debug? (format t "~A~%" (aref (swpvec-current bodies) 0)))))))

(defun sim-while-plotting (start-state &key (duration 3600) (dt 20) (bodies-to-plot '() bodies-to-plot-p))

    (let* ((bodies-to-plot (if bodies-to-plot-p
                               bodies-to-plot
                               (coerce 
                                 (range 
                                   (length 
                                     (swpvec-current start-state))) 
                                 'list))))
      (close-all-plots)
      (simulating-n-bodies-while-doing (start-state duration dt))
      (format t "~A~%" bodies-to-plot)))


