;;
;; Example of Izhikevich regular spiking neuron
;;


(use sundials srfi-4)


(define  reltol 1.0e-6)
(define  abstol 1.0e-6)


(define  t0 0.0)
(define  tf 240.0)

(define k1     0.04)
(define k2     5.0)
(define k3     140.0)
(define theta 30.0)

;; State initial values 
(define  V      -65.0)

(define  RS_a 0.02)
(define  RS_b 0.2)
(define  RS_c -65.0)
(define  RS_d 8.0)

(define  U      (* RS_b V))

(define Iext 10.)

(define (subthreshold t yy)
  (let ((V (f64vector-ref yy 0))
	(U (f64vector-ref yy 1)))
    (let (
	  (dv (+ (* k1 V V) (* k2 V) k3 (- U) Iext))
	  (du (* RS_a (- (* RS_b V) U)))
	  )
      (f64vector dv du)
      )))


(define (vreset t h yy)
  (let ((V (f64vector-ref yy 0))
	(U (f64vector-ref yy 1)))
  (let ((v1 RS_c)
	(u1 (+ U RS_d)))
    (f64vector v1 u1)
    )))
    


(define (spike-detect t yy)
  (let ((v (f64vector-ref yy 0)))
    (f64vector (floor (- v theta)))))
	     

(define (main)
  
  (let ((yy (f64vector V U)))

    ;; CVODE initialization 
    (let ((h 0.1)
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
