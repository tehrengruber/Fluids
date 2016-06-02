#|
 This file is a part of Flint
 Author: Nicolas Hafner <shinmera@tymoon.eu>
|#

(in-package #:cl-user)
(asdf:defsystem fluids
  :version "1.0.0"
  :license "Artistic"
  :author "Till Ehrengruber <till@ehrengruber.ch>"
  :maintainer "Till Ehrengruber <till@ehrengruber.ch>"
  :description "Simple 3d fluid simulation"
  :homepage "https://github.com/tehrengruber/Fluids"
  :serial T
  :components ((:file "package")
               (:file "simulation"))
  :depends-on (:3d-vectors))
