;;
;; Hodgkin-Huxley pulse propagation model
;;

(import mathh sundials srfi-4 (chicken format))

(define neg -)
(define pow expt)

(define TEND  250.0)

  	                   
;; Model parameters

(define (I_stim t) 10)
(define C_m       1)
(define E_Na      50)
(define E_K       -77)
(define E_L       -54.4)
(define gbar_Na   120)
(define gbar_K    36)
(define g_L       0.3)

;; Rate functions

(define (amf v)   (* 0.1    (/ (+ v 40)  (- 1.0 (exp (/ (neg (+ v 40)) 10))))))
(define (bmf v)   (* 4.0    (exp (/ (neg (+ v 65)) 18))))
(define (ahf v)   (* 0.07   (exp (/ (neg (+ v 65)) 20))))
(define (bhf v)   (/ 1.0    (+ 1.0 (exp (/ (neg (+ v 35)) 10)))))
(define (anf v)   (* 0.01   (/ (+ v 55) (- 1 (exp (/ (neg (+ v 55)) 10))))))
(define (bnf v)   (* 0.125  (exp (/ (neg (+ v 65)) 80))))
  	                   
;; Model equations

(define (rhs t yy)

  (let ((v (f64vector-ref yy 0))
	(m (f64vector-ref yy 1))
	(h (f64vector-ref yy 2))
	(n (f64vector-ref yy 3)))

    ;; transition rates at current step
    (let ((am  (amf v))
	  (an  (anf v))
	  (ah  (ahf v))
	  (bm  (bmf v))
	  (bn  (bnf v))
	  (bh  (bhf v))

	  (g_Na (* gbar_Na  (* h (pow m 3))))
	  (g_K  (* gbar_K   (pow n 4))))
      
      (let (

	    ;; currents
	    (I_Na   (* (- v E_Na) g_Na))
	    (I_K    (* (- v E_K)  g_K))
	    (I_L    (* g_L  (- v E_L))))
		  
	(let (
	      ;; state equations
	      (dm (- (* am (- 1 m))  (* bm m)))
	      (dh (- (* ah (- 1 h))  (* bh h)))
	      (dn (- (* an (- 1 n))  (* bn n)))
	      (dv (/ (- (I_stim t) (+ I_L I_Na I_K)) C_m))
	      )
    
	  (f64vector dv dm dh dn)
	  
	  )))
    ))


(define (rhs/unsafe t yy yp _)
  (let ((v (pointer-f64-ref yy))
	(m (pointer-f64-ref (pointer+-f64 yy 1)))
	(h (pointer-f64-ref (pointer+-f64 yy 2)))
	(n (pointer-f64-ref (pointer+-f64 yy 3))))


    ;; transition rates at current step
    (let ((am  (amf v))
	  (an  (anf v))
	  (ah  (ahf v))
	  (bm  (bmf v))
	  (bn  (bnf v))
	  (bh  (bhf v))

	  (g_Na (* gbar_Na  (* h (pow m 3))))
	  (g_K  (* gbar_K   (pow n 4))))
      
      (let (

	    ;; currents
	    (I_Na   (* (- v E_Na) g_Na))
	    (I_K    (* (- v E_K)  g_K))
	    (I_L    (* g_L  (- v E_L))))
		  
	(let (
	      ;; state equations
	      (dm (- (* am (- 1 m))  (* bm m)))
	      (dh (- (* ah (- 1 h))  (* bh h)))
	      (dn (- (* an (- 1 n))  (* bn n)))
	      (dv (/ (- (I_stim t) I_L I_Na I_K) C_m))
	      )
    
	  (pointer-f64-set! yp dv)
	  (pointer-f64-set! (pointer+-f64 yp 1) dm)
	  (pointer-f64-set! (pointer+-f64 yp 2) dh)
	  (pointer-f64-set! (pointer+-f64 yp 3) dn)
	  
	  )))
    ))



(define (main)
  
  (let ((yy (f64vector -65  0.052 0.596 0.317)) ;; v m h n

	;; Integration limits 
	(t0  0.0)
	(tf  TEND)
	(dt  1e-1))
    
    ;; CVODE initialization 
    (let ((solver (cvode-create-solver
		   t0 yy rhs  
                   tstop: tf
		   abstol: 1e-4
		   reltol: 1e-4)))

      ;; In loop, call CVodeSolve, print results, and test for error. 
      
      (let recur ((tnext (+ t0 dt)) (iout 1))

	(let ((flag  (cvode-solve solver tnext)))
	  (if (negative? flag) (error 'main "CVODE solver error" flag))

	  ;; (print-results/cvode solver)

	  (if (< tnext tf)
	      (recur (+ tnext dt) (+ 1 iout)))
	  ))
      

      (let ((yy (cvode-yy solver)))
	(let ((v (f64vector-ref yy 0))
	      (m (f64vector-ref yy 1))
	      (h (f64vector-ref yy 2))
	      (n (f64vector-ref yy 3)))
	  (printf "v = ~A m = ~A h = ~A n = ~A ~%" v m h n)
	  ))
      
      (cvode-destroy-solver solver)
      
      )))

(define (main/unsafe)
  
  (let ((yy (f64vector -65  0.052 0.596 0.317)) ;; v m h n

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

	  (if (< tnext tf)
	      (recur (+ tnext dt) (+ 1 iout)))
	  ))
      

      (let ((yy (cvode-yy solver)))
	(let ((v (f64vector-ref yy 0))
	      (m (f64vector-ref yy 1))
	      (h (f64vector-ref yy 2))
	      (n (f64vector-ref yy 3)))
          (printf "v = ~A m = ~A h = ~A n = ~A ~%" v m h n)
	  ))
      
      (cvode-destroy-solver solver)
      
      )))


(define (ida-main)

  (define (ressc t yy yp)
    (let ((dd (rhs t yy)))

      (let ((dv (- (f64vector-ref yp 0) (f64vector-ref dd 0) ))
	    (dm (- (f64vector-ref yp 1) (f64vector-ref dd 1) ))
	    (dh (- (f64vector-ref yp 2) (f64vector-ref dd 2) ))
	    (dn (- (f64vector-ref yp 3) (f64vector-ref dd 3) ))
            )
	(f64vector dv dm dh dn))))
  
  (let* ((yy (f64vector -65 0.052 0.596 0.317)) ;; v m h n
	 (yp (rhs 0.0 yy))
	 
	 ;; Integration limits 
	 (t0  0.0)
	 (tf  TEND)
	 (dt  1e-1))
    

    ;; IDA initialization 
    (let ((solver (ida-create-solver t0 yy yp ressc  
                                     tstop: tf
				     abstol: 1e-4
				     reltol: 1e-4)))

      ;; In loop, call IDASolve, print results, and test for error. 
      
      (let recur ((tnext (+ t0 dt)) (iout 1))

	(let ((flag  (ida-solve solver tnext)))
	  (if (negative? flag) (error 'main "IDA solver error" flag))

	  (print-results/ida solver)

	  (if (< tnext tf)
	      (recur (+ tnext dt) (+ 1 iout)))
	  ))
      (print-results/ida solver)
      (ida-destroy-solver solver)
      
      )))



(define (print-results/cvode solver )
  (let ((yy (cvode-yy solver))
	(t  (cvode-t solver)))
    (printf "~A ~A ~A ~A ~A~%" 
	    t 
	    (f64vector-ref yy 0)
	    (f64vector-ref yy 1)
	    (f64vector-ref yy 2)
	    (f64vector-ref yy 3)
	    )))
      

(define (print-results/ida solver )
  (let ((yy (ida-yy solver))
	(t  (ida-t solver)))
    (printf "~A ~A ~A ~A ~A~%" 
	    t 
	    (f64vector-ref yy 0)
	    (f64vector-ref yy 1)
	    (f64vector-ref yy 2)
	    (f64vector-ref yy 3)
	    )))
      
      
(ida-main)
;(main)
;(main/unsafe)
