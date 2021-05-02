
;;
;; Simple integrate-and-fire example to illustrate integrating over
;; discontinuities.
;;

(import sundials srfi-4  (chicken format))


(define  reltol 1.0e-3)
(define  abstol 1.0e-4)


(define  t0 0.0)
(define  tf 100.0)



(define gL        0.2 )
(define vL      -70.0 )
(define Isyn     20.0 )
(define C         1.0 )
(define theta    25.0 )
(define vreset  -65.0 )
(define trefractory  5.0 )


(define (subthreshold t yy)
  (let ((v (f64vector-ref yy 0)))
    (f64vector (/ (+ (* (- gL) (- v vL)) Isyn) C))))

(define (spike-detect t yy . rest)
  (let ((v (f64vector-ref yy 0)))
    (f64vector (- v theta))
    ))
  
(define (ressc t yy yp)
  (let ((dd (subthreshold t yy)))
    (f64vector (- (f64vector-ref yp 0) (f64vector-ref dd 0) ))))


(define (main)
  
  (let ((yy (f64vector -65.0)))

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
	      ;; discrete event: re-initialize solver to integrate over the discontinuity
	      (cvode-reinit-solver solver (+ tnext trefractory) (f64vector vreset)) )
	  
	  (print-results solver)

	  (if (< tnext tf)
	      (recur (if (= 1 flag) (+ (+ tnext trefractory) dt) (+ tnext dt)) (+ 1 iout)))
	  ))
      
      )))


(define (print-results solver)
  (let ((yy (cvode-yy solver)))
    (printf "~A ~A ~A ~A ~A~%" 
	    (cvode-t solver)
	    (f64vector-ref yy 0)
	    (cvode-get-last-order solver)
	    (cvode-get-num-steps solver)
	    (cvode-get-last-step solver)
	    )))


(define (main-ida)
  
  (let* ((yy (f64vector -65.0))
         (yp (subthreshold 0.0 yy))
         (events (s32vector 0)))

    ;; IDA initialization 
    (let ((dt 1e-1)
          (solver (ida-create-solver t0 yy yp ressc  
                                     tstop: tf
				     abstol: 1e-6
				     reltol: 1e-7
                                     events: events
                                     residual-event: spike-detect
                                     alg-or-diff: (s32vector 0)
                                     )))
      
      ;; In loop, call CVodeSolve, print results, and test for error. 
      
      (let recur ((tnext (+ t0 dt)) (iout 1))

	(let ((flag  (ida-solve solver tnext)))

	  (if (negative? flag) (error 'main "IDA solver error" flag))

	  (if (= 1 flag)
              (begin
                ;; discrete event: re-initialize solver to integrate over the discontinuity
                (ida-reinit-solver solver 
                                   (+ trefractory tnext )
                                   (f64vector vreset)
                                   (f64vector 0.0)
                                   ))
              )

	  (print-results-ida solver)

	  (if (< tnext tf)
	      (recur (if (= 1 flag) (+ (+ tnext trefractory) dt) (+ tnext dt)) (+ 1 iout)))
	  ))
      
      )))

(define (print-results-ida solver)
  (let ((yy (ida-yy solver))
        (t (ida-t solver)))
    (printf "~A ~A~%" 
	    t
	    (f64vector-ref yy 0)
	    )))

      
(main)
(main-ida)
