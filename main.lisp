;; Enable maximum debug settings, can change to space and speed for better performance later
(declaim (optimize (speed 0) (space 0) (debug 3)))

;;(ql:quickload :trivial-time)

(load "random.lisp")
(load "graphs.lisp")
(load "actions.lisp")

;; Load target graph and pre-compute its inverse for fixed-point checks.
(defparameter *target* #2A((0 1 1 0)
                           (1 0 1 0)
                           (1 1 0 0)
                           (0 0 0 0)))

(defparameter *inverse-g-star* (aops:each #'- (ones-m (array-dimension *target* 0)) *target*))

(defun fixed-point-p (action graph)
  (equalp (lla:asum (aops:each #'* *inverse-g-star* (funcall action graph))) 0)
  )




;; Macros to modify actions to accept thresholds.
(defmacro thresholds ())


;; Testing code
(print (aops:vectorize-reduce #'+ (*inverse-g-star*) *inverse-g-star*))
(print (lla:asum *inverse-g-star*))



(defparameter *a1t* (action1 *target*))
(defparameter *actions* #(action1 action2 action3 action4 action5 action6 action7))


(defparameter *res1* (map-array *actions* (lambda (x) (funcall x *target*))))

(print (cumsum (aops:each #'scal (rand-pmf 7) *res1*)))

(sum (lla:dot (ones len) start))


(print (rand-pmf 7))

(print (action1 *target*))
(print (action2 *target*))
(print (action3 *target*))
(print (action4 *target*))
(print (action5 *target*))
(print (action6 *target*))
(print (action7 *target*))

(trivial-time:benchmark (10000) (action7 *target*))

(print (vector-sample-n (action4 *target*) 100))
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
