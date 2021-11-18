(in-package #:comp-phys)

(defmacro simulating-3-restricted-while-doing ((start-state &optional (duration 3600) (dt 20) (?debug? nil)) &body body)
    `(let* ((sim-steps (round (/ ,duration ,dt)))
           (bodies (copy-celestial-swpvec ,start-state)))
      (progn 
        (dotimes (i sim-steps)
          ,@body
          (setf (body-pos (aref (swpvec-next bodies) 0)) (pos-earth (* i ,dt)))
          (setf (body-v (aref (swpvec-next bodies) 0)) (v-earth (* i ,dt)))
          (setf (body-pos (aref (swpvec-next bodies) 1)) (pos-moon (* i ,dt)))
          (setf (body-v (aref (swpvec-next bodies) 1)) (v-moon (* i ,dt)))
          (iterate-n-body-r4k (aref (swpvec-current bodies) 2) 
                              (aref (swpvec-next bodies) 2)
                              (swpvec-current bodies)
                              dt
                              (vector 2))
          (swpvec-swap bodies)
          (setf (body-pos (aref (swpvec-next bodies) 2)) (pos-earth (* i ,dt)))
          (when ,?debug? (format t "~A~%" (aref (swpvec-current bodies) 0)))))
      (assert (not (eqswpvec ,start-state bodies)))))

(defun sim-3-restricted-while-plotting (initial-state &key (duration (days 1)) (dt 20) (bodies-to-plot '() bodies-to-plot-p ) (n-plot-points 1000))
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
      
      (simulating-3-restricted-while-doing (copied-state duration dt)
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
                          (aref timevec plot-index)
                          ;(vz (body-pos 
                                ;(aref (swpvec-current bodies)
                                      ;body-index)))
                          
                          )))))
      (close-all-plots)
      ;(format-plot nil "set zlabel \"Time elapsed (s)\"")
      (format-plot nil "set view equal xy")
      ;(format-plot nil "set size ratio -1")
      (format-plot nil "set size square")
      (format-plot nil  "unset surfaces")
      (format-plot nil  "set isosamples 100,100")
      (format-plot nil  "set pm3d implicit at st")

      (eval `(3d-plot 
              ,@(iter (for body-index in-vector bodies-to-plot)
                      (collect 
                        (make-array n-plot-points
                                    :displaced-to plot-data 
                                    :displaced-index-offset 
                                    (* (+ 0 (* 3 body-index))
                                       n-plot-points)))
                      (collect 
                        (make-array n-plot-points
                                    :displaced-to plot-data 
                                    :displaced-index-offset 
                                    (* (+ 1 (* 3 body-index))
                                       n-plot-points)))
                      ;(collect (map (type-of timevec) #'cos timevec))
                      ;(collect (map (type-of timevec) #'sin timevec))
                      (collect 
                        (make-array n-plot-points
                                    :displaced-to plot-data 
                                    :displaced-index-offset 
                                    (* (+ 2 (* 3 body-index))
                                       n-plot-points)))
                      ;(collect timevec)
                      (collect (body-name (aref (swpvec-current copied-state) body-index))))

              ))
      
      ))

(defparameter *puppet-em-initial*
    (vector (make-body (pos-earth 0)
                       (v-earth 0)
                       *mass-earth*
                       "Earth")
            (make-body (pos-moon 0)
                       (v-moon 0)
                       *mass-moon*
                       "Moon")
            (make-body (vec 2d8 0 0)
                       (vec 0 1000 0)
                       10
                       "Rocket")))
(defparameter *puppet-em*
    (make-swpvec *puppet-em-initial*
                 *puppet-em-initial*))

(defun pr-run ()
    (sim-3-restricted-while-plotting *puppet-em*
                                     :duration (days 1) 
                                     :n-plot-points 1000
                        ;:bodies-to-plot '(0 1)

                        ))
