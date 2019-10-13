;;;; comp-phys.asd

(asdf:defsystem #:comp-phys
  :description "Describe comp-phys here"
  :author "Your Name <your.name@example.com>"
  :license  "Specify license here"
  :version "0.0.1"
  :serial t
  :depends-on (:vgplot :rtg-math :iterate)
  :components ((:file "package")
  			   (:file "rad")
			   (:file "rocket-to-moon")
               (:file "comp-phys")))
