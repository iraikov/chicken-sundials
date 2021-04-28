
;;
;; Simple integrate-and-fire example to illustrate integrating over
;; discontinuities.
;;


(import mathh sundials srfi-4 (chicken format))



(define  t0 0.0)
(define  tf 100.0)


(define w         13.77)
(define R         1.0 )
(define tau      20.0 )
(define theta    25.0 )
(define vreset   10.0 )
(define trefractory  2.0 )
(define tausyn   0.1 )
(define e   2.7238 )


(define (subthreshold t yy)
  (let ((v (f64vector-ref yy 0))
        (a (f64vector-ref yy 1))
        (b (f64vector-ref yy 2)))
    (f64vector (- (* R a) v) 
               (/ (+ (- a) b) tausyn)
               (/ (- b) tausyn)
               )))

  
(define (ressc t yy yp)
  (let ((dd (subthreshold t yy)))
    (let (
          (v (- (f64vector-ref yp 0) (f64vector-ref dd 0) ))
          (a (- (f64vector-ref yp 1) (f64vector-ref dd 1) ))
          (b (- (f64vector-ref yp 2) (f64vector-ref dd 2) ))
          )
      (f64vector v a b))))


(define (event-detect t yy yp)
  (let ((V (f64vector-ref yy 0)))
    (f64vector (- V theta) (- t 5.0))
    ))
    

(define (main)
  
  (let* ((yy (f64vector 0.0 0.0 0.0))
         (yp (subthreshold 0.0 yy))
         (events (s32vector 0 0)))

    ;; CVODE initialization 
    (let ((dt 1e-1)
          (solver (ida-create-solver t0 yy yp ressc  
                                     tstop: tf
				     abstol: 1e-6
				     reltol: 1e-7
                                     events: events
                                     residual-event: event-detect
                                     alg-or-diff: (s32vector 1 0)
                                     )))
      
      ;; In loop, call CVodeSolve, print results, and test for error. 
      
      (let recur ((tnext (+ t0 dt)) (iout 1))

	(let ((flag  (ida-solve solver tnext)))

	  (if (negative? flag) (error 'main "IDA solver error" flag))

	  (if (= 1 flag)
              (begin
                ;; discrete event: re-initialize solver to integrate over the discontinuity
                (cond ((> (s32vector-ref events 0) 0)
                       (ida-reinit-solver solver 
                                          (+ trefractory tnext )
                                          (f64vector vreset (f64vector-ref yy 1) (f64vector-ref yy 2))
                                          (f64vector 0.0 (f64vector-ref yp 1) (f64vector-ref yp 2))
                                          ))
                      ((> (s32vector-ref events 1) 0)
                       (ida-reinit-solver solver 
                                          (ida-t solver)
                                          (f64vector (f64vector-ref yy 0) (f64vector-ref yy 1) (+ w (f64vector-ref yy 2)))
                                          (f64vector 0.0 0.0 0.0)
                                          ))
                      )
                ))

	  (print-results solver)

	  (if (< tnext tf)
	      (recur (if (= 1 flag) (+ (+ tnext trefractory) dt) (+ tnext dt)) (+ 1 iout)))
	  ))
      
      )))


(define (print-results solver)
  (let ((yy (ida-yy solver))
        (t (ida-t solver)))
    (printf "~A ~A ~A ~A~%" 
	    t
	    (f64vector-ref yy 0)
	    (f64vector-ref yy 1)
	    (f64vector-ref yy 2)
	    )))

      
(main)
