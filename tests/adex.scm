
;;
;; Adaptive Exponential Integrate-and-Fire model.
;;

(use sundials srfi-4)


(define  reltol 1.0e-4)
(define  abstol 1.0e-4)


(define  t0 0.0)
(define  tf 100.0)


(define C       200.0 )
(define gL      10.0 )
(define EL      -58.0 )
(define VT      -50.0 )
(define Delta   2.0 )
(define theta   0.0 )
(define trefractory  0.25)

(define a       2.0 )
(define tau_w  120.0 )
(define b       100.0 )
(define Vr     -46.0 )

(define Isyn     210.0 )

(define (subthreshold t yy)
  (let ((V (f64vector-ref yy 0))
        (W (f64vector-ref yy 1)))
    (f64vector (/ (+ (* (- gL ) (- V EL))
                     (* gL Delta (exp (/ (- V VT) Delta)))
                     (- W) Isyn) 
                  C)
               (/ (- (* a (- V EL)) W) tau_w)
               )))

(define (ressc t yy yp)
  (let ((dd (subthreshold t yy)))
    (let ((v (- (f64vector-ref yp 0) (f64vector-ref dd 0) ))
          (w (- (f64vector-ref yp 1) (f64vector-ref dd 1) )))
      (f64vector v w))))


(define (spike-detect t yy)
  (let ((V (f64vector-ref yy 0)))
    (f64vector (- V theta))
    ))

(define (residual-spike-detect t yy yp)
  (spike-detect t yy))


(define (main)
  
  (let ((yy (f64vector -70.0 0.0)))

    ;; CVODE initialization 
    (let ((dt 1e-1)
	  (solver (cvode-create-solver
		   t0 yy subthreshold  
		   tstop: tf
		   abstol: abstol
		   reltol: reltol
		   events: (s32vector 0)
		   event-fn: spike-detect
		   )))
      
      ;; In loop, call CVodeSolve, print results, and test for error. 
      
      (let recur ((tnext (+ t0 dt)) (iout 1))

	(let ((flag  (cvode-solve solver tnext)))

	  (if (negative? flag) (error 'main "CVODE solver error" flag))

	  (if (= 1 flag)
              (let ((yy (cvode-yy solver)))
                ;; discrete event: re-initialize solver to integrate over the discontinuity
                (cvode-reinit-solver solver 
                                     (+ tnext trefractory) 
                                     (f64vector Vr (+ (f64vector-ref yy 1) b))
                                     )))
	  
	  (print-results solver)

	  (if (< tnext tf)
	      (recur (if (= 1 flag) (+ (+ tnext trefractory) dt) (+ tnext dt)) (+ 1 iout)))
	  ))
      
      )))


(define (ida-main)
  (let* ((yy (f64vector -70.0 0.0))
	 (yp (subthreshold 0.0 yy))
	 
	 ;; Integration limits 
	 (t0  0.0)
	 (dt  1e-2))
    
    ;; IDA initialization 
    (let ((solver (ida-create-solver t0 yy yp ressc  
                                     tstop: tf
				     abstol: 1e-6
				     reltol: 1e-6
                                     events: (s32vector 0)
                                     residual-event: residual-spike-detect
                                     )))

      ;; In loop, call IDASolve, print results, and test for error. 
      
      (let recur ((tnext (+ t0 dt)) (iout 1))

	(let ((flag  (ida-solve solver tnext)))
	  (if (negative? flag) (error 'main "IDA solver error" flag))

	  (print-results/ida solver tnext)

	  (if (= 1 flag)
              (let ((yy (ida-yy solver)))
                ;; discrete event: re-initialize solver to integrate over the discontinuity
                (ida-reinit-solver solver 
                                   tnext 
                                   (f64vector Vr (+ (f64vector-ref yy 1) b))
                                   (f64vector 0.0 0.0)
                                   )))

	  (if (< tnext tf)
	      (recur (+ tnext dt) (+ 1 iout)))
	  ))
      
      (ida-destroy-solver solver)
      
      ))
  )


(define (print-results solver)
  (let ((yy (cvode-yy solver)))
    (printf "~A ~A ~A ~A ~A ~A~%" 
	    (cvode-t solver)
	    (f64vector-ref yy 0)
	    (f64vector-ref yy 1)
	    (cvode-get-last-order solver)
	    (cvode-get-num-steps solver)
	    (cvode-get-last-step solver)
	    )))


(define (print-results/ida solver t)
  (let ((yy (ida-yy solver)))
    (printf "~A ~A ~A~%" 
	    t
	    (f64vector-ref yy 0)
	    (f64vector-ref yy 1)
	    )))

      
(ida-main)
