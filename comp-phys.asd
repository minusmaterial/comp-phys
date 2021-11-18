;;;; comp-phys.asd

(asdf:defsystem #:comp-phys
  :description "Describe comp-phys here"
  :author "Your Name <your.name@example.com>"
  :license  "Specify license here"
  :version "0.0.1"
  :serial t
  :depends-on (:magicl :vgplot :3d-vectors :3d-matrices :iterate)
  :components ((:file "package")
  			   (:file "tools" )
  			   (:file "rad")
			   (:file "rocket-to-moon")
			   (:file "extension")
               (:file "comp-phys")
			   (:file "pr3bp")
			   ))
