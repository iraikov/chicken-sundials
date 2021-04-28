;;
;; Compute Euler's number
;;


(import mathh sundials srfi-4 (chicken format) (chicken memory))


;; Problem Constants 

(define NEQ 1)

(define TEND  1.0)


(define (ressc t yy yp rr data)
  (let ((v (- (f64vector-ref yp 0) (f64vector-ref yy 0))))
    ;(print "yy = " yy " yp = " yp " v = " v)
    (f64vector v)))


(define (ressc/unsafe t yy yp rr data)
  (let ((v (- (pointer-f64-ref yp) (pointer-f64-ref yy))))
    (pointer-f64-set! rr v)
    ))


(define (main)
  
  (let ((yy (f64vector 1.0))
	(yp (f64vector 1.0))

	 ;; Integration limits 
	 (t0  0.0)
	 (tf  TEND)
	 (dt  1e-2))
    
    ;; IDA initialization 
    (let ((solver (ida-create-solver t0 yy yp ressc  
                                     tstop: tf
				     abstol: 1e-14
				     reltol: 1e-14)))

      ;; In loop, call IDASolve, print results, and test for error. 
      
      (let recur ((tnext (+ t0 dt)) (iout 1))

	(let ((flag  (ida-solve solver tnext)))
	  (if (negative? flag) (error 'main "IDA solver error" flag))

	  (print-results solver)

	  (if (< tnext tf)
	      (recur (+ tnext dt) (+ 1 iout)))
	  ))

      (let ((yy (ida-yy solver)))
        (print yy)
	(assert (< (abs (- 2.71828182846 (f64vector-ref yy 0) )) 1e-12)) )
      
      (ida-destroy-solver solver)
      
      )))


(define (main/unsafe)
  
  (let ((yy (f64vector 1.0))
	(yp (f64vector 1.0))

	 ;; Integration limits 
	 (t0  0.0)
	 (tf  TEND)
	 (dt  1e-2))
    
    ;; IDA initialization 
    (let ((solver (ida-create-solver/unsafe t0 yy yp ressc/unsafe
                                            tstop: tf
					    abstol: 1e-14
					    reltol: 1e-14)))

      ;; In loop, call IDASolve, print results, and test for error. 
      
      (let recur ((tnext (+ t0 dt)) (iout 1))

	(let ((flag  (ida-solve solver tnext)))
	  (if (negative? flag) (error 'main "IDA solver error" flag))

	  (print-results solver)

	  (if (< tnext tf)
	      (recur (+ tnext dt) (+ 1 iout)))
	  ))

      (let ((yy (ida-yy solver)))
	(assert (< (abs (- 2.71828182846 (f64vector-ref yy 0) )) 1e-12)) )
      
      (ida-destroy-solver solver)
      
      )))


(define (print-results solver)
  (let ((yy (ida-yy solver))
	(t (ida-t solver)))
    (printf "~A ~A ~A ~A~%" 
	    t
	    (f64vector-ref yy 0)
	    (ida-get-last-order solver)
	    (ida-get-num-steps solver)
	    (ida-get-last-step solver)
	    )))
      
      
;(main)
(main/unsafe)
