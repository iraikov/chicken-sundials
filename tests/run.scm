
(import scheme (chicken base) (chicken process-context) (chicken pathname))
(define prefix (pathname-directory (program-name)))

(load (make-pathname prefix "e.scm"))
(load (make-pathname prefix "ml.scm"))
(load (make-pathname prefix "hh.scm"))
(load (make-pathname prefix "iaf.scm"))
(load (make-pathname prefix "izhfs.scm"))
