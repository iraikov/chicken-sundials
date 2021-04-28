;; -*- Hen -*-

(import (chicken base) (chicken string) (chicken format) (chicken process)
        (chicken process-context) (chicken port) srfi-13 compile-file)
(define args (command-line-arguments))

(define (sundials-try-compile header ldflags cppflags)
  (and (try-compile 
	(string-append "#include <stdlib.h>\n"
		       "#include <stdio.h>\n"
	 	       "#include <math.h>\n"
		       header "\n" 
		       "int main(int argc, char **argv) { void *x; x = NULL; CVodeFree(&x); IDAFree(&x); return 0; }\n")
	ldflags: ldflags
	cflags: cppflags
        verbose: #t)
       (cons ldflags cppflags)
       ))

(define-syntax sundials-test 
  (syntax-rules ()
    ((_ (flags ...))
     (condition-case (sundials-try-compile flags ...)
		     (t ()    #f)))))

(define sundials-headers 
#<<EOF
#include <sundials/sundials_config.h>
#include <ida/ida.h>
#include <ida/ida_direct.h>
#include <cvode/cvode.h>           
#include <cvode/cvode_diag.h>      
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
EOF
)


(define ld+cpp-options
  (or (sundials-test (sundials-headers "-lsundials_ida -lsundials_cvode -lsundials_nvecserial -lblas -llapack " 
				  "-no-pie"))
      (sundials-test (sundials-headers "-lsundials_ida -lsundials_cvode -lsundials_nvecserial -lblas -llapack" 
				  "-I/usr/include/ida -I/usr/include/cvode"))
      (sundials-test (sundials-headers "-lsundials_ida -lsundials_cvode -lsundials_nvecserial -lblas -llapack" 
				  "-I/usr/include/sundials"))
      (sundials-test (sundials-headers "-lsundials_ida -lsundials_cvode -lsundials_nvecserial -lblas -llapack" 
				  "-I/opt/local/include"))
      (error "unable to figure out location of SUNDIALS")))



(define (detect-sundials-version ldflags)
  (let ((version-env (get-environment-variable "SUNDIALS_VERSION")))
    (cond (version-env 
           (string-split version-env "."))
          (else
           (begin
             (compile-file "detect-sundials-version.scm"
                           options:  `("-d2" ,(sprintf "-L \"~A\"" ldflags))
                           output-file: "detect-sundials-version"
                           load: #f
                           verbose: #t)
             (let ((version-str (condition-case 
                                 (with-output-to-string
                                   (lambda () (load "./detect-sundials-version" )))
                                 [var () #f])))
               (string-split version-str ".")))
           ))
    ))


(define sundials-version (detect-sundials-version (car ld+cpp-options)))

(print "SUNDIALS version is " sundials-version)

(define cmd (intersperse (append args (list (sprintf "-L \"~A\"" (car ld+cpp-options)) 
                                            (sprintf "-C \"~A\"" (cdr ld+cpp-options))
                                            (sprintf "-C \"-D SUNDIALS_VERSION_MAJOR=~A\"" (car sundials-version))
                                            (sprintf "-C \"-D SUNDIALS_VERSION_MINOR=~A\"" (cadr sundials-version)))
                                            ) " "))
(system (string-concatenate cmd))
