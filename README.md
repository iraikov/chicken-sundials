# sundials

## Description

The Chicken `sundials` library provides bindings to the solvers from
the SUNDIALS library (http://computation.llnl.gov/casc/sundials/)
. SUNDIALS (SUite of Nonlinear and DIfferential/ALgebraic equation
Solvers) is a collection of solvers for systems of ordinary
differential equations and differential-algebraic equations.

The Chicken `sundials` library provides interfaces to the CVODE and
IDA solvers and has been tested with SUNDIALS versions 2.4.0 and 2.5.0.


## Library procedures

### IDA solver interface


<procedure>(ida-create-solver TSTART TSTOP VARIABLES DERIVATIVES RESIDUAL-MAIN [RESIDUAL-INIT] [RESIDUAL-EVENT] [EVENTS] [ALG-OR-DIFF]   [SUPPRESS] [IC] [USER-DATA] [RELTOL] [ABSTOL]) => IDA-SOLVER</procedure>

Creates and initializes an object representing a problem to be solved
with the IDA solver.  

Arguments `TSTART` and `TSTOP` must be real numbers that represent
the beginning and end of the independent variable range.

Arguments `VARIABLES` and `DERIVATIVES` must be SRFI-4
`f64vector` objects that hold respectively the initial values and
derivatives of the system variables.

Argument `RESIDUAL-MAIN` is used to compute the residual function
`F` and must be a procedure of the following form:

 (LAMBDA T YY YP DATA)

or

 (LAMBDA T YY YP)

depending on whether the `USER-DATA` optional argument is set, where 

; `T`  :  real-valued independent variable
; `YY` : SRFI-4 `f64vector` with current variable values
; `YP` : SRFI-4 `f64vector` with current variable derivatives
; `DATA` : is a user data object (if set)

This procedure must return a SRFI-4 `f64vector` containing the
residual vector.

Optional keyword argument `RESIDUAL-EVENT` must be a procedure of
the same form as `RESIDUAL-MAIN`, which computes a rootfinding
problem to be solved during the integration of the system. It is set
only if argument `EVENTS` is given.

Optional keyword argument `EVENTS` is an SRFI-4 `s32vector` that
is used for storage of root finding solutions. It must be given if
`RESIDUAL-EVENT` is given.

Optional keyword argument `ALG-OR-DIFF` must be an SRFI-4
`s32vector` which indicates the algebraic and differential variables
in the system. A value of 1 indiciates differential variable, and a
value of 0 indicates an algebraic one. This is required if the
`SUPPRESS` argument is given and true.

Optional keyword argument `SUPPRESS` is a boolean flag that
indicates whether algebraic variables must be suppressed in the local
error test. If it is true (suppress), then the argument
`ALG-OR-DIFF` must be given.

Optional keyword argument `IC` is a boolean flag that indicates
whether the solver must calculate consistent initial conditions, or
whether it must use the initial conditions given by `VARIABLES`.

Optional keyword argument `USER-DATA` is an object that will be
passed as an additional argument to the residual functions.

Optional keyword arguments `RELTOL` and `ABSTOL` specify relative
and absolute error tolerance, respectively. These both default to
1e-4.


<procedure>(ida-reinit-solver IDA-SOLVER T0 Y0 YP0)</procedure>

Re-initializes IDA for the solution of a problem.

<procedure>(ida-destroy-solver IDA-SOLVER)</procedure>

Deallocates the memory associated with the given solver.

<procedure>(ida-solve IDA-SOLVER T)</procedure>

Integrates the system over an interval in the independent
variable. This procedure returns either when the given `T` is
reached, or when a root is found.

<procedure>(ida-yy IDA-SOLVER)</procedure>

Returns the vector of current state values of the system.

<procedure>(ida-yp IDA-SOLVER)</procedure>

Returns the vector of current state derivative values of the system.

<procedure>(ida-get-last-order IDA-SOLVER)</procedure>

Returns the order used during the last solver step.

<procedure>(ida-get-last-step IDA-SOLVER)</procedure>

Returns the steps size used during the last solver step.

<procedure>(ida-get-num-steps IDA-SOLVER)</procedure>

Returns the cumulative number of steps taken by the solver.


### CVODE solver interface


<procedure>(cvode-create-solver TSTART TSTOP VARIABLES RHS-FN [LMM] [ITER] [EWT-FN] [EVENT-FN] [EVENTS] [USER-DATA] [RELTOL] [ABSTOL]) => CVODE-SOLVER</procedure>

Creates and initializes an object representing a problem to be solved
with the CVODE solver.  

Arguments `TSTART` and `TSTOP` must be real numbers that represent
the beginning and end of the independent variable range.

Arguments `VARIABLES` must be a SRFI-4 `f64vector` object that
holds the initial values of the system variables.

Argument `RHS-FN` is used to compute the right-hand side of the
equations, and must be a procedure of the following form:

 (LAMBDA T YY DATA)

or

 (LAMBDA T YY)

depending on whether the `USER-DATA` optional argument is set, where 

; `T`  :  real-valued independent variable
; `YY` : SRFI-4 `f64vector` with current variable values
; `DATA` : is a user data object (if set)

This procedure must return a SRFI-4 `f64vector` containing the
residual vector.

Optional keyword argument `EWT-FN` must be a procedure of the same
form as `(LAMBDA YY)`, which computes error weights for the system
variables, and which can be used in place of relative and absolute
error tolerance.

Optional keyword argument `EVENT-FN` must be a procedure of
the same form as `RHS-FN`, which computes a rootfinding
problem to be solved during the integration of the system. It is set
only if argument `EVENTS` is given.

Optional keyword argument `EVENTS` is an SRFI-4 `s32vector` that
is used for storage of root finding solutions. It must be given if
`EVENT-FN` is given.

Optional keyword argument `LMM` specifies the linear multistep
method to be used and can be one of `cvode-lmm/adams` (default) or
`cvode-lmm/bdf`. `cvode-lmm/bdf` is recommended for stiff
problems.

Optional keyword argument `ITER` specifies the iteration type to be
used and can be one of `cvode-iter/functional` (default) or
`cvode-iter/newton`. `cvode-iter/newton` is recommended for stiff
problems.

Optional keyword argument `USER-DATA` is an object that will be
passed as an additional argument to the residual functions.

Optional keyword arguments `RELTOL` and `ABSTOL` specify relative
and absolute error tolerance, respectively. These both default to
1e-4. They are only set of `EWT-FN` is not specified.


<procedure>(cvode-reinit-solver CVODE-SOLVER T0 Y0 YP0)</procedure>

Re-initializes CVODE for the solution of a problem.

<procedure>(cvode-destroy-solver CVODE-SOLVER)</procedure>

Deallocates the memory associated with the given solver.

<procedure>(cvode-solve CVODE-SOLVER T)</procedure>

Integrates the system over an interval in the independent
variable. This procedure returns either when the given `T` is
reached, or when a root is found.

<procedure>(cvode-yy CVODE-SOLVER)</procedure>

Returns the vector of current state values of the system.


## Example

```scheme

 ;;
 ;; Hodgkin-Huxley model
 ;;
 
 (use mathh sundials srfi-4)
 
 (define neg -)
 (define pow expt)
 
 (define TEND  500.0)
 
   	                   
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
 
 ;; State functions
 
 (define (minf v) (* 0.5 (+ 1 (tanh (/ (- v v1) v2)))))
 (define (winf v) (* 0.5 (+ 1 (tanh (/ (- v v3) v4)))))
 (define (lamw v) (* phi (cosh (/ (- v v3) (* 2 v4)))))
   	                   
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
 	      (dv (/ (- (I_stim t) I_L I_Na I_K) C_m))
 	      )
    
 	  (f64vector dv dm dh dn)
 	  
 	  )))
     ))
   
  (let ((yy (f64vector -65  0.052 0.596 0.317)) ;; v m h n
 
 	;; Integration limits 
 	(t0  0.0)
 	(tf  TEND)
 	(dt  1e-2))
    
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
  
          (print-results solver tnext)
 
 	  (if (< tnext tf)
 	      (recur (+ tnext dt) (+ 1 iout)))
 	  ))
       
 
      (cvode-destroy-solver solver)
       
 (define (print-results solver t)
   (let ((yy (cvode-yy solver)))
     (printf "~A ~A ~A ~A ~A ~A ~A ~A~%" 
 	    t 
 	     (f64vector-ref yy 0)
	     (f64vector-ref yy 1)
	     (f64vector-ref yy 2)
	     (f64vector-ref yy 3)
	     (cvode-get-last-order solver)
	     (cvode-get-num-steps solver)
	     (cvode-get-last-step solver)
	     )))
```

## License


  Copyright 2011-2016 Ivan Raikov.
  All rights reserved.
  
  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:
  
  Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
  
  Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.
  
  Neither the name of the author nor the names of its contributors may
  be used to endorse or promote products derived from this software
  without specific prior written permission.
  
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
  COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
  OF THE POSSIBILITY OF SUCH DAMAGE.
