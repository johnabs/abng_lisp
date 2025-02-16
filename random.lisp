(ql:quickload :distributions)
(ql:quickload :lla)
(ql:quickload :alexandria)
(ql:quickload :array-operations/all)


(defun ones (length)
  (make-array length :initial-element 1))

(defun ones-m (side-length)
  (make-array (list side-length side-length) :initial-element 1))

(defparameter *rv-uniform* (distributions:r-uniform 0 1))

(defun rand-pmf (len)
  (let* (
         (start (aops:generate (lambda () (random 1.0)) (list len)))
         (sum (lla:dot (ones len) start))
         )
    (aops:each (alexandria:rcurry #'/ sum) start)
    )
  )

(defun cumsum (array)
  (let* ((len (array-dimension array 0))
         (result (make-array len :element-type (array-element-type array))))
    (setf (aref result 0) (aref array 0))
    (loop for i from 1 below len do
      (setf (aref result i)
            (+ (aref result (1- i)) (aref array i))))
    result))

(defun vector-sample-n (weights n)
  (let ((u01 (distributions:r-uniform 0 1))
        (probs (cumsum (aops:flatten weights)))
        )
    (coerce  (loop for i in (alexandria:iota n)
                   collect (position-if (lambda (x) (> x (distributions:draw u01))) probs :from-end nil)
                   ) 'vector)
    )
  )
