;;
;; Morris-Lecar model
;;

(use mathh sundials srfi-4)


(define TEND  150.0)

  	                   
;; Model parameters

(define con  -0.4)
(define t0 0.0)
(define y0 1.75)
(define TEND 5.0)

(define (exact t) (* y0 (exp (* con (- t t0)))))

;; Solve the test problem dy/dt = -t^3 y^3 

;;   which has the exact solution y = 1/sqrt(C + t^4/2)
;;   where C = 1/y_0^2 - t_0^4/2



;; State functions

(define (rhs t yy)
  (let ((y (f64vector-ref yy 0)))
    (let ((dy (* con y)))
    
      (f64vector dy)

      ))
  )


(define (cvode-main tol)
  
  (let ((yy (f64vector y0))

	;; Integration limits 
	(t0  0.0)
	(tf  TEND)
	(dt  0.1))
    
    ;; CVODE initialization 
    (let ((solver (cvode-create-solver
		   t0 yy rhs
                   tstop: tf
		   abstol: tol
		   reltol: tol)))

      ;; In loop, call CVodeSolve, print results, and test for error. 

      (let recur ((tnext (+ t0 dt)) (iout 1))

	(let ((flag  (cvode-solve solver tnext)))
	  (if (negative? flag) (error 'main "CVODE solver error" flag))


	  (if (< tnext tf)
	      (recur (+ tnext dt) (+ 1 iout))
              (print-results/cvode solver)
              )
	  ))
      
      (cvode-destroy-solver solver)
      
      )))


(define (print-results/cvode solver)
  (let ((yy (cvode-yy solver))
        (t (cvode-t solver)))
    (printf "~A ~A ~A ~A~%" 
            t
	    (f64vector-ref yy 0)
	    (- (exact t) (f64vector-ref yy 0))
            (cvode-get-last-step solver)
	    )))
      
(cvode-main 1e-6)
(cvode-main 1e-8)
(cvode-main 1e-10)

