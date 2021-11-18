(in-package #:comp-phys)

(defun pos-to-rotating (pos rate curtime)
    (vec3 (+ (* (vx pos) (cos (* -1 rate curtime)))
             (* -1 (vy pos) (sin (* -1 rate curtime))))
          (+ (* (vx pos) (sin (* -1 rate curtime)))
             (* (vy pos) (cos (* -1 rate curtime))))
          (vz pos)))

(defun rotating-to-pos (pos rate curtime)
(vec3 (+ (* (vx pos) (cos (* -1 rate curtime)))
             (* -1 (vy pos) (sin (* -1 rate curtime))))
          (+ (* (vx pos) (sin (* -1 rate curtime)))
             (* (vy pos) (cos (* -1 rate curtime))))
          (vz pos))
  )

(defun vel-to-rotating (pos vel rate curtime)
    (v- vel
        (vc (vec 0 0 rate)
            (pos-to-rotating pos rate curtime))))

(defun rotating-to-vel (pos vel rate )
    (v+ vel
        (vc (vec 0 0 rate)
            pos)))

(defun orbital-v (mass distance)
    (sqrt (/ (* *big-g* mass)
             distance)))

(defun rocket-energy (rocket curtime)
    (+ (* (expt (vlength (body-v rocket)) 2)
          (body-mass rocket)
          1/2)
       (* -1
          *big-g* 
          (body-mass rocket)
          *mass-earth*
          (/ 1
             (vlength (v- (body-v rocket)
                          (pos-earth curtime)))))
       (* -1
          *big-g* 
          (body-mass rocket)
          *mass-moon*
          (/ 1
             (vlength (v- (body-v rocket)
                          (pos-moon curtime)))))))

(defun hm-v-1 (mass r1 r2)
    (sqrt (* *big-g*
             mass
             (- (/ 2 r1)
                (/ 2 (+ r1 r2))))))

(defun hm-dv-1 (mass r1 r2)
         (- (hm-v-1 mass r1 r2)
            (orbital-v mass r1)))

(defun hm-v-2 (mass r1 r2)
    (sqrt (* *big-g*
             mass
             (- (/ 2 r2)
                (/ 2 (+ r1 r2))))))

(defun hm-dv-2 (mass r1 r2)
    (- (orbital-v mass r2)
       (hm-v-2 mass r1 r2)))

(defun hm-transfer-period (mass r1 r2)
    (* pi
       (sqrt (/ (expt (/ (+ r1 r2) 2) 
                      3)
                (* *big-g* mass)))))

(defun mu (mass-1 mass-2)
    (/ mass-2
       (+ mass-1 mass-2)))

(defun k (rad-1 rad-2 pos)
    (/ pos
       (+ rad-1 rad-2)))

(defun rho (mu k)
    (+(* mu 
         (expt (abs (+ k -1 mu)) -3))
      (* (- 1 mu)
         (expt (abs (+ k mu)) -3))))

(defun a (mu k)
  (+ (* 2 (rho mu k))
     1))

(defun b (mu k)
  (- (rho mu k)
     1))


