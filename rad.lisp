(in-package #:comp-phys)


(defstruct rad-sample
           mean-t
           N)

(defparameter *sample-t-nil*
              (make-rad-sample :mean-t 100
                               :N 1000))

(defun rad-decay-rate (sample)
    (* -1.0 (/ (rad-sample-n sample)
               (rad-sample-mean-t sample))))

(defun rad-next-n (sample &optional (timestep 0.01))
    (+ (rad-sample-n sample) 
       (* (rad-decay-rate sample) 
          timestep)))


    
(declaim (inline rad-decay))
(defun rad-decay (sample-t0  
                  time-length
                  &optional (timestep 0.01d0))
    (let (
          (time-length (coerce time-length 'long-float))
          (timestep (coerce timestep 'long-float)))
      (simulate-rad-decay sample-t0 time-length timestep)))

(declaim (inline simulate-rad-decay))
(defun simulate-rad-decay (sample-t0  
                           time-length
                           &optional (timestep 0.01l0))
    (declare (optimize (speed 3) (safety 0)))
    (declare (type double-float time-length timestep))
    (declare (inline * ) )
    (declare (inline +))
    (declare (inline aref))
    (declare (inline rad-next-n) (inline rad-sample-n))
    
    (let* ((steps (round (/ time-length timestep)))
          (results (make-array steps :element-type 'double-float))
          (time-array (make-array steps :element-type 'double-float))
          (current-sample (make-rad-sample 
                            :mean-t (rad-sample-mean-t sample-t0)
                            :n (rad-sample-n sample-t0))))
      
      (let ((j 0l0))
        (dotimes (i (length time-array))
          (setf (aref time-array i) (* j timestep))
          (setf j (+ j 1l0))))

      (dotimes (i (length results)) 
        (setf (rad-sample-n current-sample)
              (rad-next-n current-sample timestep))
        (setf (aref results i) (rad-sample-n current-sample)))
      (list time-array results)))

(defun rad-timestep (sample &optional (dt 0.01))
    (make-rad-sample :mean-t (rad-sample-mean-t sample)
                     :N (+ (rad-sample-n sample)
                           (* (rad-decay-rate sample)
                              dt))))
