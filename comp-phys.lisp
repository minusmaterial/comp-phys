;;;; comp-phys.lisp

(in-package #:comp-phys)

(defun expand-lists (lists)
    (if (eq nil lists)
        nil
        `(,@(car lists) ,@(expand-lists (cdr lists))))) 

(defun derivative (f x &optional (dx 0.1) (max-diff 0.01))
    (let ((der (derivative-once f x dx)))
      (print der)
      (if (< (abs (- der
                     (derivative-once f x (* dx 0.5))))
             max-diff)
          der
          (derivative f x (* dx 0.25) max-diff))
      ))

(defun derivative-once (f x dx)
    (/ (- (funcall f (+ x dx)) 
                     (funcall f x))
                   dx)) 

(defun test-plot ()
    (plot (range 0 20) 
          (map 'vector #'sin 
               #(0 1 2 3 4 5 6 7 8 9 10 11 12))
          ";curve;")
    (title "Curve"))

(defmacro plot-functions (&rest fns)
    `(3d-plot ,@(expand-lists (mapcar #'eval fns))))

(defun helix-plot (&optional (o 0))
    (let ((z (range 0 10 0.1))) 
      (list (map 'vector #'sin (map 'vector (lambda (x) (+ x o)) z))
            (map 'vector #'cos (map 'vector (lambda (x) (+ x o)) z))
            z)))


