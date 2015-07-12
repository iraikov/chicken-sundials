;;
;; Example of Izhikevich fast spiking neuron
;;


(use sundials srfi-4)


(define  reltol 1.0e-6)
(define  abstol 1.0e-6)


(define  t0 0.0)
(define  tf 100.0)

(define k     1.0)
(define Vpeak 25.0)
(define Vt    -55.0)
(define Vr    -40.0)
(define Vb    -55.0)
(define Cm    20.0)
(define Isyn  0.0)

;; State initial values 
(define  V      -65.0)
(define  tspike   0.0)
(define  spike     #f)

(define  FS_a 0.2)
(define  FS_b 0.025)
(define  FS_c -45.0)
(define  FS_U (* FS_b V))

(define Iext 100.)

(define (UU v)  (if (< v Vb) 0. (* FS_b (- v Vb) (- v Vb) (- v Vb))))

(define (subthreshold t yy)
  (let ((v (f64vector-ref yy 0))
	(u (f64vector-ref yy 1)))
    (let (
	  (dv (/ (+ (* k (- v Vr) (- v Vt)) (- u) Iext) Cm))
	  (du (* FS_a (- (UU v) u)))
	  )
      (f64vector dv du)
      )))


(define (vreset t h yy)
  (let ((v (f64vector-ref yy 0))
	(u (f64vector-ref yy 1)))
  (let ((v1 FS_c)
	(u1 (+ u (* h FS_a (- (* FS_b (- FS_c (* FS_b (- v Vb) (- v Vb) (- v Vb)) )) u)))))

    (f64vector v1 u1)
    )))
    


(define (spike-detect t yy)
  (let ((v (f64vector-ref yy 0)))
    (f64vector (floor (- v Vpeak)))))
	     

(define (main)
  
  (let ((yy (f64vector V FS_U)))

    ;; CVODE initialization 
    (let ((h 1e-1)
	  (solver (cvode-create-solver
		   t0 yy subthreshold  
		   tstop: tf
		   abstol: abstol
		   reltol: reltol
		   events: (s32vector 0)
		   event-fn: spike-detect
		   )))
      
      ;; In loop, call CVodeSolve, print results, and test for error. 
      
      (let recur ((tnext (+ t0 h)) (iout 1))

	(let ((flag  (cvode-solve solver tnext)))

	  (if (negative? flag) (error 'main "CVODE solver error" flag))

	  (if (= 1 flag)
	      (begin
		;; discrete event: re-initialize solver to integrate over the discontinuity
		(print-results solver)
		(cvode-reinit-solver solver (+ tnext h) (vreset (cvode-t solver) h (cvode-yy solver)) )
		(recur tnext (+ iout 1)))

	      (if (< (cvode-t solver) tnext) 
		  (recur tnext (+ 1 iout))
		  (begin

		    (print-results solver)
		    (if (< tnext tf)
			(recur (if (= 1 flag) (+ (+ tnext h) h) (+ tnext h)) (+ 1 iout)))
		    )))
	  ))
      ))
  )

(define (print-results solver)
  (let ((yy (cvode-yy solver))
	(t (cvode-t solver)))
    (let ((spike (>= (f64vector-ref (spike-detect t yy) 0) 0)))
    (printf "~A ~A ~A ~A ~A ~A ~A ~A~%" 
	    t (or (and spike 1) 0) (or (and spike t) 0.0)
	    (f64vector-ref yy 0)
	    (f64vector-ref yy 1)
	    (cvode-get-last-order solver)
	    (cvode-get-num-steps solver)
	    (cvode-get-last-step solver)
	    ))
    ))

      
(main)
