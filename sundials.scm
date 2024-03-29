;;
;;
;; Chicken Scheme bindings for SUNDIALS (SUite of Nonlinear and
;; DIfferential/ALgebraic equation Solvers).
;;
;;  Copyright 2011-2021 Ivan Raikov.
;;
;; 
;;  Redistribution and use in source and binary forms, with or without
;;  modification, are permitted provided that the following conditions
;;  are met:
;; 
;;  - Redistributions of source code must retain the above copyright
;;  notice, this list of conditions and the following disclaimer.
;; 
;;  - Redistributions in binary form must reproduce the above
;;  copyright notice, this list of conditions and the following
;;  disclaimer in the documentation and/or other materials provided
;;  with the distribution.
;; 
;;  - Neither name of the copyright holders nor the names of its
;;  contributors may be used to endorse or promote products derived
;;  from this software without specific prior written permission.
;; 
;;  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND THE
;;  CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
;;  INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
;;  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
;;  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR THE
;;  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
;;  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
;;  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
;;  USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
;;  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
;;  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
;;  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
;;  POSSIBILITY OF SUCH DAMAGE.
;;
;;  
;;


(module sundials

	(
	 ida-create-solver ida-create-solver/unsafe
         ida-reinit-solver ida-destroy-solver ida-solve ida-yy ida-yp ida-t
	 ida-get-last-order ida-get-num-steps ida-get-last-step

	 cvode-lmm/adams cvode-lmm/bdf
	 cvode-iter/functional cvode-iter/newton
	 cvode-create-solver cvode-create-solver/unsafe
	 cvode-reinit-solver cvode-destroy-solver cvode-solve cvode-yy cvode-t
	 cvode-get-last-order cvode-get-num-steps cvode-get-last-step

         pointer+f64
	 )

	(import scheme (chicken base) (chicken foreign) (chicken fixnum) srfi-4 srfi-69
                (only (chicken memory) move-memory! 
                      pointer+ pointer-f64-ref pointer-f64-set!))



#>

#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <ida/ida.h>
#include <ida/ida_direct.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <cvode/cvode.h>           
#include <cvode/cvode_diag.h>      
#include <nvector/nvector_serial.h>

#define C_f64vector_i_ref(b, i)      (((double *)C_data_pointer(C_block_item((b), 1)))[ i ])

static void chicken_panic (C_char *) C_noret;
static void chicken_panic (C_char *msg)
{
  C_word *a = C_alloc (C_SIZEOF_STRING (strlen (msg)));
  C_word scmmsg = C_string2 (&a, msg);
  C_halt (scmmsg);
  exit (5); /* should never get here */
}

static void chicken_throw_exception(C_word value) C_noret;
static void chicken_throw_exception(C_word value)
{
  char *aborthook = C_text("\003syserror-hook");

  C_word *a = C_alloc(C_SIZEOF_STRING(strlen(aborthook)));
  C_word abort = C_intern2(&a, aborthook);

  abort = C_block_item(abort, 0);
  if (C_immediatep(abort))
    chicken_panic(C_text("`##sys#error-hook' is not defined"));

#if defined(C_BINARY_VERSION) && (C_BINARY_VERSION >= 8)
  C_word rval[3] = { abort, C_SCHEME_UNDEFINED, value };
  C_do_apply(3, rval);
#else
  C_save(value);
  C_do_apply(1, abort, C_SCHEME_UNDEFINED);
#endif

}

void chicken_error (char *msg, C_word obj) 
{
  size_t msglen;
  C_word *a;
  C_word scmmsg;
  C_word list;

  msglen = strlen (msg);
  a = C_alloc (C_SIZEOF_STRING (msglen) + C_SIZEOF_LIST(2));
  scmmsg = C_string2 (&a, (char *) msg);
  list = C_list(&a, 2, scmmsg, obj);
  chicken_throw_exception(list);
}

typedef struct Solver_User_Data_struct
{
    unsigned int data_index;
} Solver_User_Data;

typedef struct IDA_Solver_Handle_struct
{
    void* ida_mem;
    double tret;
    C_word syy;
    C_word syp;
    N_Vector yy;
    N_Vector yp;
    N_Vector id;
    int event_number;
    int* events;
    double abstol;
    double reltol;
    SUNMatrix A;
    SUNLinearSolver LS;
    Solver_User_Data *user_data;
    
} IDA_Solver_Handle;


typedef struct CVODE_Solver_Handle_struct
{
    void* cvode_mem;
    double tret;
    C_word syy;
    N_Vector yy;
    N_Vector id;
    int event_number;
    int* events;
    double abstol;
    double reltol;
    Solver_User_Data *user_data;

} CVODE_Solver_Handle;

<#

(define (pointer+f64 x n) (pointer+ x (* 8 n)))
(define-foreign-type IDASolverHandle "IDA_Solver_Handle")
(define-foreign-type CVODESolverHandle "CVODE_Solver_Handle")


(define cvode-lmm/adams (foreign-value "CV_ADAMS" int))
(define cvode-lmm/bdf (foreign-value "CV_BDF" int))

(define cvode-iter/functional (foreign-value "CV_FUNCTIONAL" int))
(define cvode-iter/newton (foreign-value "CV_NEWTON" int))

(define ida-residual-init-global  (make-hash-table))
(define ida-residual-main-global  (make-hash-table))
(define ida-residual-event-global (make-hash-table))

(define ida-data-global (make-hash-table))

(define cvode-rhs-global    (make-hash-table))
(define cvode-event-global  (make-hash-table))
(define cvode-ewt-global    (make-hash-table))

(define cvode-data-global (make-hash-table))



(define c_nvector_data_pointer 
  (foreign-safe-lambda* (nonnull-c-pointer double) ((nonnull-c-pointer x))
#<<EOF
   N_Vector v = (N_Vector)x;
   C_return (N_VGetArrayPointer(v));
EOF
))

(define c_user_data_index 
  (foreign-safe-lambda* unsigned-int ((nonnull-c-pointer x))
#<<EOF
   C_return (((Solver_User_Data *)x)->data_index);
EOF
))


(define-external (ida_residual_main_cb (double t) 
				   (c-pointer yy) 
				   (c-pointer yp) 
				   (c-pointer rr)
				   (c-pointer data)) int
   (let ((data-index (c_user_data_index data)))
     (let ((v (and ((hash-table-ref ida-residual-main-global data-index)
                    t (c_nvector_data_pointer yy) 
                    (c_nvector_data_pointer yp) (c_nvector_data_pointer rr )
                    (hash-table-ref/default ida-data-global data-index #f)) 0)))
       v))
   )


(define-external (ida_residual_init_cb (double t) 
				   (c-pointer yy) 
				   (c-pointer yp) 
				   (c-pointer rr)
				   (c-pointer data)) int 
   (let ((data-index (c_user_data_index data)))
     (let ((v ((hash-table-ref ida-residual-init-global data-index)
               t (c_nvector_data_pointer yy) 
               (c_nvector_data_pointer yp) (c_nvector_data_pointer rr )
               (hash-table-ref/default ida-data-global data-index #f))))
       v))
   )



(define-external (ida_residual_event_cb (double t) 
                                        (c-pointer yy) 
                                        (c-pointer yp) 
                                        (c-pointer rr)
                                        (c-pointer data)) int
   (let ((data-index (c_user_data_index data)))
     (let ((v (and ((hash-table-ref ida-residual-event-global data-index)
                    t (c_nvector_data_pointer yy) 
                    (c_nvector_data_pointer yp) rr
                    (hash-table-ref/default ida-data-global data-index #f)) 0)))
       v))
   )


(define-external (cvode_rhs_cb (double t) 
			       (c-pointer yy) 
			       (c-pointer yp) 
			       (c-pointer data)) int
   (let ((data-index (c_user_data_index data)))
     (let ((v (and ((hash-table-ref cvode-rhs-global data-index)
                    t (c_nvector_data_pointer yy) (c_nvector_data_pointer yp) 
                    (hash-table-ref/default cvode-data-global data-index #f)) 0)))
       v))
   )


(define-external (cvode_event_cb (double t) 
				 (c-pointer yy) 
				 (c-pointer gout) 
				 (c-pointer data)) int
   (let ((data-index (c_user_data_index data)))
     (let ((v (and ((hash-table-ref cvode-event-global data-index)
                    t (c_nvector_data_pointer yy) gout
                    (hash-table-ref/default cvode-data-global data-index #f)) 0)))
       v))
   )


(define-external (cvode_ewt_cb (c-pointer yy) 
			       (c-pointer ewt) 
			       (c-pointer data)) int
   (let ((data-index (c_user_data_index data)))
     (let ((v (and ((hash-table-ref cvode-ewt-global data-index)
                    (c_nvector_data_pointer yy) (c_nvector_data_pointer ewt) 
                    (hash-table-ref/default cvode-data-global data-index #f)) 0)))
       v))
   )


#>

extern int  ida_residual_init_cb(double,void *,void *,void *,void *);
extern int ida_residual_main_cb(double, void *, void *, void *, void *);

extern  int  ida_residual_event_cb(double,void *,void *,void *,void *);

extern  int  cvode_rhs_cb(double,void *,void *,void *);
extern  int  cvode_event_cb(double,void *,void *,void *);
extern  int  cvode_ewt_cb(void *,void *,void *);


void adjust_zero_crossings (N_Vector v, double abstol)
{
    int i;
    for (i = 0; i < NV_LENGTH_S(v); i++)
    {
        if (fabs(NV_Ith_S(v,i)) < abstol) NV_Ith_S(v,i) = 0;
    }
    return;
}


IDA_Solver_Handle* ida_create_solver
( double time_start,
  double time_stop,

  int variable_number,
  C_word variables, 
  C_word derivatives, 
  int ic,
  int* alg_or_diff,
  int suppress,

  int event_number, 
  int* events, 

  unsigned int data_index,

  double abstol,
  double reltol
)
{
    IDA_Solver_Handle* solver_handle;
    assert ((solver_handle = malloc (sizeof(struct IDA_Solver_Handle_struct))) != NULL);
    assert((solver_handle->user_data = malloc (sizeof(struct Solver_User_Data_struct))) != NULL);

    solver_handle->user_data->data_index = data_index;
    solver_handle->tret = 0.0;
    solver_handle->event_number = event_number;
    solver_handle->events = events;

    int flag = 0;
    int i = 0;

    solver_handle->syy = variables;
    solver_handle->syp = derivatives;

    solver_handle->yy = N_VNew_Serial(variable_number);
    solver_handle->yp = N_VNew_Serial(variable_number);
    solver_handle->id = N_VNew_Serial(variable_number);

    for (i = 0; i < variable_number; i++)
    {
        NV_Ith_S(solver_handle->id,i) = (alg_or_diff[i] == 0) ? 0.0 : 1.0;
        NV_Ith_S(solver_handle->yy,i) = C_f64vector_i_ref(variables, i);
        NV_Ith_S(solver_handle->yp,i) = C_f64vector_i_ref(derivatives, i);
    }

    solver_handle->abstol = abstol;
    solver_handle->reltol = reltol;

    solver_handle->ida_mem = IDACreate();
    if (solver_handle->ida_mem == NULL) 
    {
       chicken_error("could not allocate memory with IDACreate", C_SCHEME_UNDEFINED);
    }

    flag = IDASetUserData(solver_handle->ida_mem, (void *)solver_handle->user_data);
    if (flag != IDA_SUCCESS) 
    {
       chicken_error("could not set user data with IDASetUserData", C_fix(flag));
    }

    flag = IDAInit (solver_handle->ida_mem, 
		    (IDAResFn)ida_residual_main_cb, 
		    time_start, 
		    solver_handle->yy, 
		    solver_handle->yp);
    if (flag != IDA_SUCCESS) 
        chicken_error("could not initialize solver with IDAInit", C_fix(flag)) ;	

    flag = IDASStolerances(solver_handle->ida_mem, solver_handle->reltol, solver_handle->abstol);
    if (flag != IDA_SUCCESS) 
        chicken_error("could not set tolerances with IDASStolerances", C_fix(flag)) ;	
    if (flag != IDA_SUCCESS) 
        chicken_error("could not initialize solver with IDAMalloc", C_fix(flag)) ;	

    flag = IDASetId(solver_handle->ida_mem, solver_handle->id);
    if (flag != IDA_SUCCESS) 
        chicken_error("could not set algebraic variables with IDASetId", C_fix(flag)) ;	
	
    flag = IDASetSuppressAlg(solver_handle->ida_mem, suppress);
    if (flag != IDA_SUCCESS) 
        chicken_error("could not set error suppression flag with IDASetSuppressAlg", C_fix(flag)) ;	

    /* Create dense SUNMatrix for use in linear solves */
    solver_handle->A = SUNDenseMatrix(variable_number, variable_number);
    if ((void *)solver_handle->A == NULL)
        chicken_error("could not initialize SUNDenseMatrix object", C_fix(1));


    /* Create dense SUNLinearSolver object */
    solver_handle->LS = SUNDenseLinearSolver(solver_handle->yy, solver_handle->A);
    if ((void *)solver_handle->LS == NULL)
        chicken_error("could not initialize SUNDenseLinearSolver object", C_fix(1));
    
    /* Attach the matrix and linear solver */
    flag = IDADlsSetLinearSolver(solver_handle->ida_mem, solver_handle->LS, solver_handle->A);
    if (flag != IDA_SUCCESS) 
        chicken_error("could not attach matrix and linear solver with IDADlsSetLinearSolver", C_fix(flag));	
  
    if (ic)
    {
       flag = IDACalcIC(solver_handle->ida_mem, IDA_YA_YDP_INIT, time_start + solver_handle->abstol);
      if (flag != IDA_SUCCESS) 
        chicken_error("could not compute initial conditions with IDACalcIC", C_fix(flag)) ;	
    }

    if (event_number > 0)
    { 
       flag = IDARootInit(solver_handle->ida_mem, event_number, (IDARootFn)ida_residual_event_cb);
       if (flag != IDA_SUCCESS) 
          chicken_error("could not initialize event variables with IDARootInit", C_fix(flag)) ;	

    }

    if (time_stop > 0.0)
    {
      flag = IDASetStopTime(solver_handle->ida_mem, time_stop);
    }
    if (flag != IDA_SUCCESS) 
         chicken_error("could not set stop time with IDASetStopTime", C_fix(flag)) ;	

    return solver_handle;
}


void ida_reinit_solver (IDA_Solver_Handle* solver_handle, double t0, C_word y0, C_word yp0)
{
    int flag; N_Vector ytmp, yptmp;
    size_t variable_number = NV_LENGTH_S(solver_handle->yy);

    ytmp = N_VNew_Serial(variable_number);
    yptmp = N_VNew_Serial(variable_number);

    for (size_t i = 0; i < variable_number; i++)
    {
        NV_Ith_S(ytmp,i) = C_f64vector_i_ref(y0, i);
        NV_Ith_S(yptmp,i) = C_f64vector_i_ref(yp0, i);
    }

    N_VScale(1.0, ytmp, solver_handle->yy);
    N_VScale(1.0, yptmp, solver_handle->yp);

    
    flag = IDAReInit(solver_handle->ida_mem, t0, solver_handle->yy, solver_handle->yp);

    if (flag != CV_SUCCESS) 
         chicken_error("could not set reinitialize solver time with IDAReInit", C_fix(flag)) ;	

    N_VDestroy_Serial(ytmp);
    N_VDestroy_Serial(yptmp);
    
    return;
}


void ida_destroy_solver (IDA_Solver_Handle* solver_handle)
{
    IDAFree(&(solver_handle->ida_mem));
    N_VDestroy_Serial(solver_handle->yy);
    N_VDestroy_Serial(solver_handle->yp);
    N_VDestroy_Serial(solver_handle->id);
    SUNLinSolFree(solver_handle->LS);
    SUNMatDestroy(solver_handle->A);
    free(solver_handle->user_data);
    free(solver_handle);
    return;
}


int ida_solve (IDA_Solver_Handle* solver_handle, double tout)
{
    int flag;

    if (solver_handle->event_number > 0)
    {
        long int ngevals;
        flag = IDAGetNumGEvals(solver_handle->ida_mem, &ngevals);

        if (flag != IDA_SUCCESS) 
          chicken_error("error in IDAGetNumGEvals", C_fix(flag)) ;	


        if (ngevals == 0)
        {
	 flag = IDASolve (solver_handle->ida_mem, 
			  tout,
			  &(solver_handle->tret),
			  solver_handle->yy,
			  solver_handle->yp, 
			  IDA_ONE_STEP);

            switch (flag)
            {
                case IDA_SUCCESS:      return 0;
                case IDA_ROOT_RETURN:  return 0;
                case IDA_TSTOP_RETURN: return 2;
                default:               chicken_error("unknown status code returned by IDASolve", C_fix(flag)) ;	

            }
        }
    }

    flag = IDASolve (solver_handle->ida_mem,
		     tout,
		     &(solver_handle->tret),
		     solver_handle->yy,
		     solver_handle->yp, 
		     IDA_NORMAL);

    switch (flag)
    {
     case IDA_SUCCESS:
      return 0;
      
      case IDA_ROOT_RETURN:
      flag = IDAGetRootInfo(solver_handle->ida_mem, solver_handle->events);
      if (flag != IDA_SUCCESS) 
         chicken_error("could not obtain even information with IDAGetRootInfo", C_fix(flag)) ;	
      adjust_zero_crossings(solver_handle->yy, solver_handle->abstol);
      adjust_zero_crossings(solver_handle->yp, solver_handle->abstol);
      return 1;
      
      case IDA_TSTOP_RETURN:
      return 2;
      
      default: 
      chicken_error("unknown status code returned by IDASolve", C_fix(flag)) ;	

      }
    }



CVODE_Solver_Handle *cvode_create_solver
( int lmm,
  int iter, 
  int order, 

  double time_start,
  double time_stop,

  int variable_number,
  C_word variables, 

  int ewt,

  int event_number, 
  int* events, 

  unsigned int data_index,

  double abstol,
  double reltol
)
{
    CVODE_Solver_Handle* solver_handle;
    assert ((solver_handle = malloc (sizeof(struct CVODE_Solver_Handle_struct))) != NULL);
    assert((solver_handle->user_data = malloc (sizeof(struct Solver_User_Data_struct))) != NULL);

    solver_handle->user_data->data_index = data_index;
    solver_handle->tret = 0.0;
    solver_handle->event_number = event_number;
    solver_handle->events = events;

    int flag = 0;
    int i = 0;

    solver_handle->syy = variables;
    solver_handle->yy = N_VNew_Serial(variable_number);
    for (i = 0; i < variable_number; i++)
    {
        NV_Ith_S(solver_handle->yy,i) = C_f64vector_i_ref(variables, i);
    }


    solver_handle->abstol = abstol;
    solver_handle->reltol = reltol;

    solver_handle->cvode_mem = CVodeCreate(lmm, iter);
    if (solver_handle->cvode_mem == NULL) 
    {
       chicken_error("could not allocate memory with CVodeCreate", C_SCHEME_UNDEFINED);
    }

    if (lmm == CV_BDF)
    {
       flag = CVodeSetStabLimDet (solver_handle->cvode_mem, 1);
       if (flag != CV_SUCCESS) 
       {
          chicken_error("could not set stability limit detection with CVodeSetStabLimtDet", C_fix(flag));
       }
    }

    flag = CVodeSetMaxOrd (solver_handle->cvode_mem,order);
    if (flag != CV_SUCCESS) 
    {
        chicken_error("could not set maximum order with CVodeSetMaxOrd", C_fix(flag));
    }


    flag = CVodeSetUserData(solver_handle->cvode_mem, (void *)solver_handle->user_data);

    if (flag != CV_SUCCESS) 
    {
       chicken_error("could not set user data with CVodeSetUserData", C_SCHEME_UNDEFINED);
    }

    flag = CVodeInit (solver_handle->cvode_mem, 
		    (CVRhsFn)cvode_rhs_cb, 
		    time_start, 
		    solver_handle->yy);
    if (flag != CV_SUCCESS) 
        chicken_error("could not initialize solver with CVodeInit", C_fix(flag)) ;	

    if (ewt > 0)
    {
       flag = CVodeWFtolerances(solver_handle->cvode_mem, (CVEwtFn)cvode_ewt_cb);
       if (flag != CV_SUCCESS) 
          chicken_error("could not set error weight function with CVodeWFtolerances", C_fix(flag)) ;	
    } else
    {
       flag = CVodeSStolerances(solver_handle->cvode_mem, solver_handle->reltol, solver_handle->abstol);
       if (flag != CV_SUCCESS) 
          chicken_error("could not set tolerances with CVodeSStolerances", C_fix(flag)) ;	
    }

    if (iter == CV_NEWTON)
    {
       flag = CVDiag(solver_handle->cvode_mem);

       if (flag != CV_SUCCESS) 
          chicken_error("could not initialize linear solver with CVDiag", C_fix(flag)) ;	
       
    }

    if (event_number > 0)
    { 
       flag = CVodeRootInit(solver_handle->cvode_mem, event_number, (CVRootFn)cvode_event_cb);
       if (flag != CV_SUCCESS) 
          chicken_error("could not initialize event variables with CVodeRootInit", C_fix(flag)) ;	

    }


    if (time_stop > 0.0)
    {
     flag = CVodeSetStopTime(solver_handle->cvode_mem, time_stop);
     if (flag != CV_SUCCESS) 
         chicken_error("could not set stop time with CVodeSetStopTime", C_fix(flag)) ;	
    }

    return solver_handle;
}


void cvode_reinit_solver (CVODE_Solver_Handle* solver_handle, double t0, C_word y0)
{
    int flag; N_Vector ytmp;
    size_t variable_number = NV_LENGTH_S(solver_handle->yy);
        
    ytmp = N_VNew_Serial(variable_number);
    for (size_t i = 0; i < variable_number; i++)
    {
        NV_Ith_S(ytmp,i) = C_f64vector_i_ref(y0, i);
    }

    N_VScale(1.0, ytmp, solver_handle->yy);
    flag = CVodeReInit(solver_handle->cvode_mem, t0, solver_handle->yy);
    if (flag != CV_SUCCESS) 
         chicken_error("could not set reinitialize solver time with CVodeReInit", C_fix(flag)) ;	

    N_VDestroy_Serial(ytmp);
    
    return;
}


void cvode_destroy_solver (CVODE_Solver_Handle* solver_handle)
{
    CVodeFree(&(solver_handle->cvode_mem));
    N_VDestroy_Serial(solver_handle->yy);
    free(solver_handle->user_data);
    free(solver_handle);
    return;
}


int cvode_solve (CVODE_Solver_Handle* solver_handle, double tout)
{
    int flag;

    if (solver_handle->event_number > 0)
    {
        long int nevals;
        flag = CVodeGetNumRhsEvals(solver_handle->cvode_mem, &nevals);
        if (flag != CV_SUCCESS) 
           chicken_error("error in CVodeGetNumRhsEvals", C_fix(flag)) ;	

        if (nevals == 0)
        {
	 flag = CVode (solver_handle->cvode_mem, 
		       tout,
		       solver_handle->yy,
		       &(solver_handle->tret),
		       CV_ONE_STEP);

            
            switch (flag)
            {
                case CV_SUCCESS:      return 0;
                case CV_ROOT_RETURN:  return 1;
                case CV_TSTOP_RETURN: return 2;
                default:              chicken_error("unknown status code returned by CVode", C_fix(flag)) ;	

            }
        }
    }

    flag = CVode (solver_handle->cvode_mem, 
		  tout,
		  solver_handle->yy,
		  &(solver_handle->tret),
		  CV_NORMAL);

    switch (flag)
    {
     case CV_SUCCESS:
      return 0;
      
      case CV_ROOT_RETURN:
      flag = CVodeGetRootInfo(solver_handle->cvode_mem, solver_handle->events);
      if (flag != CV_SUCCESS) 
         chicken_error("could not obtain event information wtih CVodeGetRootInfo", C_fix(flag)) ;	
      adjust_zero_crossings(solver_handle->yy, solver_handle->abstol);
      return 1;
      
      case CV_TSTOP_RETURN:
      return 2;
      
      default: 
      chicken_error("unknown status code returned by CVode", C_fix(flag)) ;	

      }
    }

<#

(define c-cvode-create-solver
  (foreign-safe-lambda (nonnull-c-pointer CVODESolverHandle) 
		  "cvode_create_solver" 
		  int     ;; linear multistep method
		  int     ;; iter
		  int     ;; max-order

		  double  ;; time_start
		  double  ;; time_stop
		  
		  int ;; variable_number
		  scheme-object ;; variables
		  
		  int       ;; ewt
		  
		  int ;; event_number
		  s32vector ;; events
		  
		  unsigned-int ;; user data
		  
		  double ;; abstol
		  double ;; reltol
		  ))

  
(define c-ida-create-solver
  (foreign-safe-lambda (nonnull-c-pointer IDASolverHandle) 
		  "ida_create_solver" 
		  double    ;; time_start
		  double    ;; time_stop
		  
		  int ;; variable_number
		  scheme-object ;; variables
		  scheme-object ;; derivatives
		  bool ;; calculate initial conditions
		  s32vector ;; alg_or_diff
		  bool
		  
		  int ;; event_number
		  s32vector ;; events

		  unsigned-int ;; user data

		  double ;; abstol
		  double ;; reltol
		  ))


(define (ida-create-solver
	 tstart variables derivatives 
	 residual-main #!key 
	 (tstop 0.0)
	 (residual-init #f)
	 (residual-event #f)
	 (alg-or-diff (make-s32vector (f64vector-length variables) 1))
	 (suppress #f)
	 (ic #f)
	 (user-data #f)
	 (events (make-s32vector 0))
	 (reltol  1.0e-6)
	 (abstol  1.0e-6)
	 )
  
  (assert (= (f64vector-length variables) (f64vector-length derivatives)))
  
  (let (
	(n (f64vector-length variables))
	(nevents (s32vector-length events))
	(data-index (hash-table-size ida-data-global))
	)

    (if user-data (hash-table-set! ida-data-global data-index user-data))

    (let ((residual-main1
	   (if user-data
	       (lambda (t yy yp rr data)  
		 (let ((yy1 (make-f64vector n))
		       (yp1 (make-f64vector n)))
		   (move-memory! yy yy1 (fx* 8 n))
		   (move-memory! yp yp1 (fx* 8 n))
		   (let ((v (residual-main t yy1 yp1 data)))
		     (move-memory! v rr (fx* 8 n))
		     )))
	       (lambda (t yy yp rr data)  
		 (let ((yy1 (make-f64vector n))
		       (yp1 (make-f64vector n)))
		   (move-memory! yy yy1 (fx* 8 n))
		   (move-memory! yp yp1 (fx* 8 n))
		   (let ((v (residual-main t yy1 yp1)))
		     (move-memory! v rr (* 8 n))
		     )))
	       ))

	  (residual-init1
	   (and residual-init
		(if user-data
		    (lambda (t yy yp rr data)  
		      (let ((yy1 (make-f64vector n))
			    (yp1 (make-f64vector n)))
			(move-memory! yy yy1 (fx* 8 n))
			(move-memory! yp yp1 (fx* 8 n))
			(let ((v (residual-init t yy1 yp1 data)))
			  (move-memory! v rr (fx* 8 n))
			  )))
		    (lambda (t yy yp rr data)  
		      (let ((yy1 (make-f64vector n))
			    (yp1 (make-f64vector n)))
			(move-memory! yy yy1 (fx* 8 n))
			(move-memory! yp yp1 (fx* 8 n))
			(let ((v (residual-init t yy1 yp1)))
			  (move-memory! v rr (* 8 n))
			  )))
		    )))

	  (residual-event1 
	   (and residual-event
		(if user-data
		    (lambda (t yy yp rr data)  
		      (let ((yy1 (make-f64vector n))
			    (yp1 (make-f64vector n)))
			(move-memory! yy yy1 (fx* 8 n))
			(move-memory! yp yp1 (fx* 8 n))
			(let ((v (residual-event t yy1 yp1 data)))
			  (move-memory! v rr (fx* 8 nevents))
			  )))
		    (lambda (t yy yp rr data)  
		      (let ((yy1 (make-f64vector n))
			    (yp1 (make-f64vector n)))
			(move-memory! yy yy1 (fx* 8 n))
			(move-memory! yp yp1 (fx* 8 n))
			(let ((v (residual-event t yy1 yp1)))
			  (move-memory! v rr (* 8 nevents))
			  )))
		    )))

	  )

      (hash-table-set! ida-residual-main-global data-index residual-main1)
      (if residual-init (hash-table-set! ida-residual-init-global data-index residual-init1))
      (if residual-event (hash-table-set! ida-residual-event-global data-index residual-event1))
    
      (c-ida-create-solver
       tstart tstop  
       n variables derivatives
       ic alg-or-diff suppress (s32vector-length events) events
       data-index abstol reltol)

      )))


(define (ida-create-solver/unsafe
	 tstart variables derivatives 
	 residual-main #!key 
	 (tstop 0.0)
	 (residual-init #f)
	 (residual-event #f)
	 (alg-or-diff (make-s32vector (f64vector-length variables) 1))
	 (suppress #f)
	 (ic #f)
	 (user-data #f)
	 (events (make-s32vector 0))
	 (reltol  1.0e-6)
	 (abstol  1.0e-6)
	 )

  (assert (= (f64vector-length variables) (f64vector-length derivatives)))
  
  (let ((data-index (hash-table-size ida-data-global)))

    (hash-table-set! ida-residual-main-global data-index residual-main)

    (if user-data (hash-table-set! ida-data-global data-index user-data))
    (if residual-init (hash-table-set! ida-residual-init-global data-index residual-init))
    (if residual-event (hash-table-set! ida-residual-event-global data-index residual-event))
    
    (c-ida-create-solver 
     tstart tstop  
     (f64vector-length variables) variables derivatives
     ic alg-or-diff suppress (s32vector-length events) events
     data-index abstol reltol)

    ))
		       
			   

(define c-ida-reinit-solver
  (foreign-safe-lambda void "ida_reinit_solver" (nonnull-c-pointer IDASolverHandle) double scheme-object scheme-object ))

(define (ida-reinit-solver solver t0 y0 yp0)
  (assert (f64vector? y0))
  (c-ida-reinit-solver solver t0 y0 yp0))

(define c-ida-destroy-solver
  (foreign-safe-lambda void "ida_destroy_solver" (nonnull-c-pointer IDASolverHandle) ))

(define (ida-destroy-solver solver)
  (c-ida-destroy-solver solver))

(define ida-solve
  (foreign-safe-lambda int "ida_solve" (nonnull-c-pointer IDASolverHandle) double ))

(define ida-get-last-order
  (foreign-safe-lambda* int (((nonnull-c-pointer IDASolverHandle) handle))
#<<EOF
    int result;
    IDAGetLastOrder (handle->ida_mem, &result);
    C_return (result);
EOF
))

(define ida-get-num-steps
  (foreign-safe-lambda* long (((nonnull-c-pointer IDASolverHandle) handle))
#<<EOF
    long result;
    IDAGetNumSteps (handle->ida_mem, &result);
    C_return (result);
EOF
))

(define ida-get-last-step
  (foreign-safe-lambda* double (((nonnull-c-pointer IDASolverHandle) handle))
#<<EOF
    double result;
    IDAGetLastStep (handle->ida_mem, &result);
    C_return (result);
EOF
))

(define ida-syy 
  (foreign-safe-lambda* scheme-object (((nonnull-c-pointer IDASolverHandle) handle))
#<<EOF
    C_return (handle->syy);
EOF
))

(define ida-syp
  (foreign-safe-lambda* scheme-object (((nonnull-c-pointer IDASolverHandle) handle))
#<<EOF
    C_return (handle->syp);
EOF
))

(define c_ida_yy 
  (foreign-safe-lambda* void (((nonnull-c-pointer IDASolverHandle) handle) (f64vector result))
#<<EOF
   memcpy(result,
	  NV_DATA_S(handle->yy),
	  (NV_LENGTH_S(handle->yy))*sizeof(double));
EOF
))

(define c_ida_yy_length
  (foreign-safe-lambda* unsigned-int (((nonnull-c-pointer IDASolverHandle) handle))
#<<EOF
   C_return (NV_LENGTH_S(handle->yy));
EOF
))

(define (ida-yy handle)
  (let ((v (make-f64vector (c_ida_yy_length handle) 0.0)))
    (c_ida_yy handle v)
    v))

(define c_ida_yp 
  (foreign-safe-lambda* void (((nonnull-c-pointer IDASolverHandle) handle) (f64vector result))
#<<EOF
   memcpy(result,
	  NV_DATA_S(handle->yp),
	  (NV_LENGTH_S(handle->yp))*sizeof(double));
EOF
))

(define c_ida_yp_length
  (foreign-safe-lambda* unsigned-int (((nonnull-c-pointer IDASolverHandle) handle))
#<<EOF
   C_return (NV_LENGTH_S(handle->yp));
EOF
))

(define (ida-yp handle)
  (let ((v (make-f64vector (c_ida_yp_length handle) 0.0)))
    (c_ida_yp handle v)
    v))
  

(define c_ida_t
  (foreign-safe-lambda* double (((nonnull-c-pointer IDASolverHandle) handle))
#<<EOF
   double result;
   result = handle->tret;
   C_return (result);
EOF
))

(define (ida-t handle)
    (c_ida_t handle))
  

(define (cvode-create-solver
	 tstart variables  
	 rhs-fn #!key 
	 (tstop 0.0)
	 (lmm cvode-lmm/adams)
	 (iter cvode-iter/functional)
	 (max-order 2)
	 (ewt-fn #f)
	 (event-fn #f)
	 (user-data #f)
	 (events (make-s32vector 0))
	 (reltol  1.0e-6)
	 (abstol  1.0e-6)
	 )

  (let (
	(n (f64vector-length variables))
	(nevents (s32vector-length events))
	(data-index (hash-table-size cvode-data-global))
	)

    (if user-data (hash-table-set! cvode-data-global data-index user-data))

    (let ((rhs-fn1
	   (if user-data
	       (lambda (t yy yp data)  
		 (let ((yy1 (make-f64vector n)))
		   (move-memory! yy yy1 (fx* 8 n))
		   (let ((v (rhs-fn t yy1 data)))
		     (move-memory! v yp (fx* 8 n))
		     )))
	       (lambda (t yy yp data)  
		 (let ((yy1 (make-f64vector n)))
		   (move-memory! yy yy1 (fx* 8 n))
		   (let ((v (rhs-fn t yy1)))
		     (move-memory! v yp (fx* 8 n))
		     )))
	       ))
	  (ewt-fn1
	   (if user-data
	       (lambda (yy ewt data)  
		 (let ((yy1 (make-f64vector n)))
		   (move-memory! yy yy1 (fx* 8 n))
		   (let ((v (ewt-fn yy1 data)))
		     (move-memory! v ewt (fx* 8 n))
		     )))
	       (lambda (yy ewt data)  
		 (let ((yy1 (make-f64vector n)))
		   (move-memory! yy yy1 (fx* 8 n))
		   (let ((v (ewt-fn yy1)))
		     (move-memory! v ewt (fx* 8 n))
		     )))
	       ))
	  (event-fn1
	   (and event-fn
		(if user-data
		    (lambda (t yy gout data)  
		      (let ((yy1 (make-f64vector n)))
			(move-memory! yy yy1 (fx* 8 n))
			(let ((v (event-fn t yy1 data)))
			  (move-memory! v gout (fx* 8 nevents))
			  )))
		    (lambda (t yy gout data)  
		      (let ((yy1 (make-f64vector n)))
			(move-memory! yy yy1 (fx* 8 n))
			(let ((v (event-fn t yy1)))
			  (move-memory! v gout (fx* 8 nevents))
			  )))
	       )))
	  )
    
      (hash-table-set! cvode-rhs-global data-index rhs-fn1)
      (if ewt-fn (hash-table-set! cvode-ewt-global data-index ewt-fn1))
      (if event-fn (hash-table-set! cvode-event-global data-index event-fn1))
    
    (c-cvode-create-solver
     lmm iter max-order tstart tstop  
     n variables
     (if ewt-fn 1 0) (s32vector-length events) events
     data-index abstol reltol)

    )))

(define (cvode-create-solver/unsafe
	 tstart variables  
	 rhs-fn #!key 
	 (tstop 0.0)
	 (lmm cvode-lmm/adams)
	 (iter cvode-iter/functional)
	 (max-order 2)
	 (ewt-fn #f)
	 (event-fn #f)
	 (user-data #f)
	 (events (make-s32vector 0))
	 (reltol  1.0e-6)
	 (abstol  1.0e-6)
	 )
  (let ((data-index (hash-table-size cvode-data-global)))

    (if user-data (hash-table-set! cvode-data-global data-index user-data))

    (hash-table-set! cvode-rhs-global data-index rhs-fn)
    (if ewt-fn (hash-table-set! cvode-ewt-global data-index ewt-fn))
    (if event-fn (hash-table-set! cvode-event-global data-index event-fn))
    
    (c-cvode-create-solver
     lmm iter max-order tstart tstop  
     (f64vector-length variables)  variables
     (if ewt-fn 1 0) (s32vector-length events) events
     data-index abstol reltol)
    ))


(define c-cvode-reinit-solver
  (foreign-safe-lambda void "cvode_reinit_solver" (nonnull-c-pointer CVODESolverHandle) double scheme-object ))

(define (cvode-reinit-solver solver t0 y0)
  (assert (f64vector? y0))
  (c-cvode-reinit-solver solver t0 y0))

		       
(define c-cvode-destroy-solver
  (foreign-safe-lambda void "cvode_destroy_solver" (nonnull-c-pointer CVODESolverHandle) ))

(define (cvode-destroy-solver solver)
  (c-cvode-destroy-solver solver))



(define cvode-solve
  (foreign-safe-lambda int "cvode_solve" (nonnull-c-pointer CVODESolverHandle) double ))
			   
(define cvode-syy 
  (foreign-safe-lambda* scheme-object (((nonnull-c-pointer CVODESolverHandle) handle))
#<<EOF
    C_return (handle->syy);
EOF
))


(define c_cvode_t
  (foreign-safe-lambda* double (((nonnull-c-pointer CVODESolverHandle) handle))
#<<EOF
   double result;
   result = handle->tret;
   C_return (result);
EOF
))

(define c_cvode_yy 
  (foreign-safe-lambda* void (((nonnull-c-pointer CVODESolverHandle) handle) (f64vector result))
#<<EOF
   memcpy(result,
	  NV_DATA_S(handle->yy),
	  (NV_LENGTH_S(handle->yy))*sizeof(double));
EOF
))

(define c_cvode_yy_length
  (foreign-safe-lambda* unsigned-int (((nonnull-c-pointer CVODESolverHandle) handle))
#<<EOF
   C_return (NV_LENGTH_S(handle->yy));
EOF
))

(define (cvode-yy handle)
  (let ((v (make-f64vector (c_cvode_yy_length handle) 0.0)))
    (c_cvode_yy handle v)
    v))

(define (cvode-t handle)
  (c_cvode_t handle))

(define cvode-get-last-order
  (foreign-safe-lambda* int (((nonnull-c-pointer CVODESolverHandle) handle))
#<<EOF
    int result;
    CVodeGetLastOrder (handle->cvode_mem, &result);
    C_return (result);
EOF
))

(define cvode-get-num-steps
  (foreign-safe-lambda* long (((nonnull-c-pointer CVODESolverHandle) handle))
#<<EOF
    long result;
    CVodeGetNumSteps (handle->cvode_mem, &result);
    C_return (result);
EOF
))

(define cvode-get-last-step
  (foreign-safe-lambda* double (((nonnull-c-pointer CVODESolverHandle) handle))
#<<EOF
    double result;
    CVodeGetLastStep (handle->cvode_mem, &result);
    C_return (result);
EOF
))




)
