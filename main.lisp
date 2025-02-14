;; Enable maximum debug settings, can change to space and speed for better performance later
(declaim (optimize (speed 0) (space 0) (debug 3)))

(ql:quickload :cl-ana)
(ql:quickload :distributions)
(ql:quickload :alexandria)
(ql:quickload :lla)
                                        ;(ql:quickload :aops)
(ql:quickload :array-operations/all)

(load "graphs.lisp")
(load "actions.lisp")


;; Action Functions



;; Macros to modify actions to accept thresholds.

;; Testing code
(defparameter *test*  '((0 1 1 0)
                        (1 0 1 0)
                        (1 1 0 0)
                        (0 0 0 0)))

(defparameter *a1t* (action1 *test*))
(defun scal (x y) (aops:each (alexandria:curry #'* x) y))
(defparameter *actions* #(action1 action2 action3 action4 action5 action6 action7))


(defparameter *res1* (map-array *actions* (lambda (x) (funcall x *test*)) ))

(print (cumsum (aops:each #'scal (rand-pmf 7) *res1*)))

(sum (lla:dot (ones len) start))

(defun rand-pmf (len)
  (let* (
         (start (aops:generate (lambda () (random 1.0)) (list len)))
         (sum (lla:dot (ones len) start))
         )
    (aops:each (alexandria:rcurry #'/ sum) start)
    )
  )

(print (aops:generate (lambda () (random 1.0)) '(7)))

(print (rand-pmf 7))


(print (scal 0.2 *a1t*))
(print *a1t*)
(print (action1 *test*))
(print (action2 *test*))
(print (action3 *test*))
(print (action4 *test*))
(print (action5 *test*))
(print (action6 *test*))
(print (action7 *test*))

(print (vector-sample-n (action 6 *test*) '(0.8 0.2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0) 1000))
;; Usage:
(defmethod map-array (array function
                      &optional (retval (make-array (array-dimensions array))))
  "Apply FUNCTION to each element of ARRAY
Return a new array, or write into the optional 3rd argument."
  (dotimes (i (array-total-size array) retval)
    (setf (row-major-aref retval i)
          (funcall function (row-major-aref array i)))))


;; Write the code to pre-filter the actions by which ones can't be feasible based on the attempted fixed-point theorem
;; and action thresholding.

;; After ruling out infeasible actions and finding thresholds for feasible ones, run the optimization on the remainder
;; with the constraint that the starting graph must still be a fixed point of the generator. Then use other standard
;; optimization techniques as before.

;; Leverage LP/NLP/MINLP problem properties on coerced generators to see how many feasible optimal solutions there are and if our other fixed points fall in that region.
