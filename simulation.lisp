;; Load 3rd party packages
(ql:quickload :3d-vectors)
(use-package :3d-vectors)

;;
;; Parameter
;;
(defclass fluid-simulation-parameter ()
  ((smoothing-length
    :initform 0.0457
    :accessor smoothing-length)
   (dt
    :initform 0.01
    :accessor dt)
   (particle-mass
    :initform 0.02
    :accessor particle-mass)
   (gas-stiffness
    :initform 5
    :accessor gas-stiffness)
   (rest-density
    :initform 998.29
    :accessor rest-density)
   (rcut
    :initform 0.01
    :accessor rcut) ; usually equal to the smoothing length
   (viscosity
    :initform 3.5
    :accessor viscosity)
   (surface-threshold
    :initform 7.065
    :accessor surface-threshold)
   (surface-tension
    :initform 0.0728
    :accessor surface-tension)
   (radius
    :initform 0.01
    :accessor radius)
   (gravity-acceleration
    :initform (vec 0 0 -9.81)
    :accessor gravity-acceleration)))

;;
;; Particle
;;
(defclass particle ()
  ((location
    :initarg :location
    :initform (vec 0 0 0)
    :accessor location)
   (velocity
    :initarg :velocity
    :initform (vec 0 0 0)
    :accessor velocity)
   (acceleration
    :initform (vec 0 0 0)
    :accessor acceleration)
   (normal
    :initform (vec 0 0 0)
    :accessor normal)
   (density
    :initform 0
    :accessor density)
   (pressure
    :initform 0
    :accessor pressure)
   (is-on-surface
    :initform NIL
    :accessor is-on-surface)))

(defmethod print-object ((particle particle) stream)
  (print-unreadable-object (particle stream :type T)
    (format stream "~a" (location particle))))

;;
;; Generate particle configuration
;;
(defun generate-particle-configuration (nx ny nz &key (box (vec 1 1 1)))
  (let ((conf (make-array (* nx ny nz) :initial-element nil))
        (dx (/ (vx box) (1- nx)))
        (dy (/ (vy box) (1- ny)))
        (dz (/ (vz box) (1- nz)))
        (pos 0))
    (dotimes (i nx conf)
      (dotimes (j ny)
        (dotimes (k nz)
          (print ())
          (setf (aref conf pos) (make-instance 'particle :location (vec (* i dx) (* j dy) (* k dz))))
          (incf pos))))))

;; calculate squared norm of the vector v
(defun squared-norm (v)
  (+ (* (vx v) (vx v))
     (* (vy v) (vy v))
     (* (vz v) (vz v))))

;;
;; Simulation
;;
(defvar *params*)
(defvar *particles*)

(defun simulation (steps)
  ;; parameters
  (setf *params* (make-instance 'fluid-simulation-parameter))
  ;; setup particles
  (setf *particles* (generate-particle-configuration 10 10 10))
  ;;
  ;; update density and pressure
  ;;
  (loop for particle across *particles* do
        ;; set density to zero
        (setf (density particle) 0)
        ;; iterate over all neighbours
        (loop for neighbour across *particles* do
              (when (< (squared-norm (v- (location particle) (location neighbour))) (expt (smoothing-length *params*) 2))
                ;; calculate density from all neighbours (including the particle itself)
                (incf (density particle)
                      (wpoly6 (squared-norm (v- (location particle) (location neighbour))) (smoothing-length *params*)))))
        ;;(incf (density particle) (wpoly6 0 (smoothing-length *params*)))
        (setf (density particle) (* (particle-mass *params*) (density particle)))
        (setf (pressure particle) (* (gas-stiffness *params*) (- (density particle) (rest-density *params*)))))
  ;;
  ;; compute forces
  ;;
  (let ((smoothing-length (smoothing-length *params*)))
    (loop for particle across *particles* do
          (let ((force-pressure (vec 0 0 0))
                (force-viscosity (vec 0 0 0))
                (force-surface (vec 0 0 0))
                (color-field-normal (vec 0 0 0))
                (color-field-laplacian 0))
            (loop for neighbour across *particles* do
                  (let* ((diff (v- (location particle) (location neighbour)))
                         (radius-squared (squared-norm diff)))
                    (unless (eql particle neighbour)
                      ;; pressure force
                      (nv+ force-pressure (v* (wspiky-gradient diff radius-squared smoothing-length)
                                              (+ (/ (pressure particle) (expt (density particle) 2))
                                                 (/ (pressure neighbour) (expt (density particle) 2)))))
                      ;; viscous force
                      (nv+ force-viscosity (v* (/ (wviscosity-laplacian radius-squared smoothing-length) (density neighbour))
                                               (v- (velocity neighbour) (velocity particle)))))
                    ;; compute colour field
                    (nv+ (v/ (wpoly6-gradient diff radius-squared smoothing-length)
                             (density neighbour)))
                    (incf color-field-laplacian (/ (wpoly6-laplacian radius-squared smoothing-length) (density neighbour)))
                    ))
            ;; collect forces
            (nv* force-pressure (* -1 (particle-mass *params*) (density particle)))
            (nv* force-viscosity (* (viscosity *params*) (particle-mass *params*)))
            (nv* color-field-normal (particle-mass *params*))
            (setf (normal particle) (v- color-field-normal))
            (setf color-field-laplacian (* (particle-mass *params*) color-field-laplacian))

            ;; surface tension force
            (let ((color-field-normal-magnitude (vlength color-field-normal)))
              (when (> color-field-normal-magnitude (surface-threshold *params*))
                (setf (is-on-surface particle) T)
                (setf force-surface (v* (* (/ (- (surface-tension *params*))
                                              color-field-normal-magnitude)
                                           color-field-laplacian)
                                        color-field-normal))))
            ;; else flag = false

            ;; add sph forces
            (setf (acceleration particle) (v+ (v/ (density particle) (v+ force-pressure force-viscosity force-surface))
                                              (gravity-acceleration *params*))))
          (let* ((dt (dt *params*)) (new-location (v+ (location particle) (v* (velocity particle) dt (v* dt dt (acceleration particle)))))
                 (new-velocity (v/ (v- new-location (location particle)) dt)))
            (setf (location particle) new-location)
            (setf (velocity particle) new-velocity)))))


;;
;; Smoothing kernel
;;
(defun wpoly6 (radius-squared smoothing-length)
  (assert (>= radius-squared 0))
  (assert (<= radius-squared (expt smoothing-length 2)))
  (* (/ 315 (* 64 pi (expt smoothing-length 9)))
     (expt (- (expt smoothing-length 2) radius-squared) 3)))

(defun wpoly6-gradient (diff radius-squared smoothing-length)
  (assert (>= radius-squared 0))
  (assert (<= radius-squared (expt smoothing-length 2)))
  (v* (* (/ -945 (* 32 pi (expt smoothing-length 9)))
         (expt (- (expt smoothing-length 2) radius-squared) 2))
      diff))

(defun wpoly6-laplacian (radius-squared smoothing-length)
  (assert (>= radius-squared 0))
  (assert (<= radius-squared (expt smoothing-length 2)))
  (* (/ -945 (* 32 pi (expt smoothing-length 9)))) (- (expt smoothing-length 2) radius-squared) (- (* 3 (expt smoothing-length 2)) (* 7 radius-squared)))

(defun wspiky-gradient (diff radius-squared smoothing-length)
  (assert (>= radius-squared 0))
  (assert (<= radius-squared (expt smoothing-length 2)))
  (let ((radius (sqrt radius-squared)))
    (v* (* (/ -45 (* pi (expt smoothing-length 6))) (expt (- smoothing-length radius) 2))
        (v/ diff radius))))

(defun wviscosity-laplacian (radius-squared smoothing-length)
  (* (/ 45 (* pi (expt smoothing-length 6)))
     (- smoothing-length (sqrt radius-squared))))
