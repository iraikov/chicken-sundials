;;
;; Morris-Lecar model
;;


(import mathh sundials srfi-4 (chicken format) (chicken memory))


(define TEND  150.0)

  	                   
;; Model parameters

(define Istim  50.0)
(define vl     -50)
(define vk     -70)
(define vca    100)
(define gl     2.0)
(define gk     8.0)
(define gca    4.0)
(define c      20.0)
(define v1     -1.0)
(define v2     15)
(define v3     10)
(define v4     14.5)
(define phi    0.0667)


;; State functions

(define (minf v) (* 0.5 (+ 1 (tanh (/ (- v v1) v2)))))
(define (winf v) (* 0.5 (+ 1 (tanh (/ (- v v3) v4)))))
(define (lamw v) (* phi (cosh (/ (- v v3) (* 2 v4)))))
  	                   
;; Model equations

(define (rhs t yy)
  (let ((v (f64vector-ref yy 0))
	(w (f64vector-ref yy 1)))

  (let ((ica (* gca (* (minf v)  (- vca v))))
	(ik  (* gk  (* w (- vk v )))))
    
    (let ((dv (/ (+ Istim (* gl (- vl v)) ica ik) c))
	  (dw (* (lamw v) (- (winf v) w))))
    
      (f64vector dv dw)

      ))
  ))

(define (rhs/unsafe t yy yp data)
  
  (let ((v (pointer-f64-ref yy))
	(w (pointer-f64-ref (pointer+f64 yy 1))))

  (let ((ica (* gca (* (minf v)  (- vca v))))
	(ik  (* gk  (* w (- vk v )))))
    
    (let ((dv (/ (+ Istim (* gl (- vl v)) ica ik) c))
	  (dw (* (lamw v) (- (winf v) w))))
    
      (pointer-f64-set! yp dv)
      (pointer-f64-set! (pointer+f64 yp 1) dw)

      0
      ))
  ))


(define (ida-main)

  (define (ressc t yy yp)
    (let ((dd (rhs t yy)))

      (let ((v (- (f64vector-ref yp 0) (f64vector-ref dd 0) ))
	    (w (- (f64vector-ref yp 1) (f64vector-ref dd 1) )))

	(f64vector v w))))
  
  (let* ((yy (f64vector -60.899 0.0149)) ;; v w
	 (yp (rhs 0.0 yy))
	 
	 ;; Integration limits 
	 (t0  0.0)
	 (tf  TEND)
	 (dt  1e-1))
    

    ;; IDA initialization 
    (let ((solver (ida-create-solver t0 yy yp ressc  
                                     tstop: tf
				     abstol: 1e-10
				     reltol: 1e-10)))

      ;; In loop, call IDASolve, print results, and test for error. 
      
      (let recur ((tnext (+ t0 dt)) (iout 1))

	(let ((flag  (ida-solve solver tnext)))
	  (if (negative? flag) (error 'main "IDA solver error" flag))

	  (print-results/ida solver tnext)

	  (if (< tnext tf)
	      (recur (+ tnext dt) (+ 1 iout)))
	  ))
      
      (ida-destroy-solver solver)
      
      )))



(define (cvode-main)
  
  (let ((yy (f64vector -60.899 0.0149));; v w

	;; Integration limits 
	(t0  0.0)
	(tf  TEND)
	(dt  0.1))
    
    ;; CVODE initialization 
    (let ((solver (cvode-create-solver
		   t0 yy rhs
                   tstop: tf
		   abstol: 1e-6
		   reltol: 1e-6)))

      ;; In loop, call CVodeSolve, print results, and test for error. 

      (let recur ((tnext (+ t0 dt)) (iout 1))

	(let ((flag  (cvode-solve solver tnext)))
	  (if (negative? flag) (error 'main "CVODE solver error" flag))

	  (print-results/cvode solver tnext)

	  (if (< tnext tf)
	      (recur (+ tnext dt) (+ 1 iout)))
	  ))
      
      (let ((yy (cvode-yy solver)))
	(let ((v (f64vector-ref yy 0))
	      (w (f64vector-ref yy 1)))
	  (printf "v = ~A w = ~A~%" v w)
	  ))

      (cvode-destroy-solver solver)
      
      )))

(define (cvode-main/unsafe)
  
  (let ((yy (f64vector -60.899 0.0149)) ;; v w

	;; Integration limits 
	(t0  0.0)
	(tf  TEND)
	(dt  0.1))
    
    ;; CVODE initialization 
    (let ((solver (cvode-create-solver/unsafe
		   t0 yy rhs/unsafe
                   tstop: tf
		   abstol: 1e-4
		   reltol: 1e-4)))

      ;; In loop, call CVodeSolve, print results, and test for error. 

      (let recur ((tnext (+ t0 dt)) (iout 1))

	(let ((flag  (cvode-solve solver tnext)))
	  (if (negative? flag) (error 'main "CVODE solver error" flag))

	  (print-results/cvode solver tnext)

	  (if (< tnext tf)
	      (recur (+ tnext dt) (+ 1 iout)))
	  ))
      
      (let ((yy (cvode-yy solver)))
	(let ((v (f64vector-ref yy 0))
	      (w (f64vector-ref yy 1)))
	  (printf "v = ~A w = ~A~%" v w)
	  ))

      (cvode-destroy-solver solver)
      
      )))


(define (print-results/cvode solver t)
  (let ((yy (cvode-yy solver)))
    (printf "~A ~A ~A ~A ~A ~A~%" 
	    t 
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
      
(cvode-main/unsafe)
(cvode-main)
(ida-main)      
