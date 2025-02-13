;; Enable maximum debug settings, can change to space and speed for better performance later
(declaim (optimize (speed 0) (space 0) (debug 3)))

;; Load relevant dependencies
                                        ;(ql:quickload :cl-ana)
                                        ;(ql:quickload :distributions)
                                        ;(ql:quickload :alexandria)


(load "graphs.lisp")


;; Action Functions
(defun action1 (adjacency-matrix)
  "Compute the degree-distribution of each node in a graph from its adjacency matrix.

   Parameters:
   adjacency-matrix - a square matrix representing the adjacency matrix of the graph
   adjacency-matrix - a square matrix representing the adjacency matrix of the graph
   adjacency-matrix - a square matrix representing the adjacency matrix of the graph

   Returns:
   a vector where the index corresponds to the node and the value is the degree of that node divided by the total degree of all nodes"
  (let*
      (

       (nodes (length adjacency-matrix))
       (deg (cl-ana.linear-algebra:euclidean-dot  (ones nodes) adjacency-matrix))
       (d2 (cl-ana.tensor:tensor-map 'nan-div (cl-ana.linear-algebra:euclidean-dot adjacency-matrix deg) deg))
       )
    (cl-ana.tensor:tensor-map #'nan-div d2 (make-list nodes :initial-element (cl-ana.linear-algebra:euclidean-dot d2 (ones nodes))))
    )

  )


(defun action2 (adjacency-matrix)
  "Compute the degree-distribution of each node in a graph from its adjacency matrix.

   Parameters:
   adjacency-matrix - a square matrix representing the adjacency matrix of the graph
   adjacency-matrix - a square matrix representing the adjacency matrix of the graph
   adjacency-matrix - a square matrix representing the adjacency matrix of the graph

   Returns:
   a vector where the index corresponds to the node and the value is the degree of that node divided by the total degree of all nodes"
  (let*
      (

       (nodes (length adjacency-matrix))
       (deg (cl-ana.linear-algebra:euclidean-dot  (ones nodes) adjacency-matrix))
       )
    (cl-ana.tensor:tensor-map #'nan-div deg (make-list nodes :initial-element (cl-ana.linear-algebra:euclidean-dot deg (ones nodes))))
    )

  )

(defun action3 (adjacency-matrix)
  "Compute the probability of forming an edge based on neighbors' page rank scores.

   Parameters:
   adjacency-matrix - a square matrix representing the adjacency matrix of the graph
   adjacency-matrix - a square matrix representing the adjacency matrix of the graph
   adjacency-matrix - a square matrix representing the adjacency matrix of the graph

   Returns:
   a vector where the index corresponds to the node and the value is the degree of that node divided by the total degree of all nodes"
  (let*
      (

       (nodes (length adjacency-matrix))
       (pr (page-rank adjacency-matrix))
       (norm (cl-ana.tensor:tensor-map #'nan-div pr  (make-list nodes :initial-element (vecsum pr)) ))
       )
    norm
    )

  )

(defun action4 (g)
  "Compute normalized betweenness-based matrix from adjacency matrix g"
  (let* ((n (length g))
         ;; Convert adjacency matrix to adjacency list for brandes
         (bc-values (brandes g)) ; get betweenness centrality
         (bc-sum (reduce #'+ bc-values))
         (norm (mapcar #'(lambda (x) (nan-div x bc-sum))
                       bc-values))
         ;; Create matrix of repeated normalized values
         (m (make-array (list n n)
                        :initial-contents
                        (loop repeat n collect norm))))

    ;; Set diagonal to 0
    (dotimes (i n)
      (setf (aref m i i) 0))

    ;; Normalize entire matrix
    (let ((matrix-sum (loop for i below n
                            sum (loop for j below n
                                      sum (aref m i j)))))
      (unless (zerop matrix-sum)
        (dotimes (i n)
          (dotimes (j n)
            (setf (aref m i j)
                  (/ (aref m i j) matrix-sum))))))

    m))

(defun action5 (adj-matrix &key (cumulative t))
  "Process adjacency matrix to find normalized second-order connections"
  (let* ((n (length adj-matrix))
         ;; Initialize A3 with zeros
         ;; Find positions where A2 - A1 > 0 and set to 1
         (diff (cl-ana.tensor:tensor-map #'g0p (cl-ana.tensor:tensor-- (cl-ana.linear-algebra:matrix-mult *test2* *test2*) *test2*)) ))
    ;; Clear diagonal
    (dotimes (i n)
      (setf (cl-ana.tensor:tensor-ref (cl-ana.tensor:tensor-ref diff i) i) 0))
    ;; Calculate row sums
    (let* (
           (sums (row-sums diff))
           ;; Normalize by row sums
           (normalized (divide-by-row-sums diff sums)))
      ;; Apply cumulative sum if requested
      (if cumulative
          (cumulative-sum-rows normalized)
          normalized))))

(defun action6 (adj-tensor)
  "Process cl-ana tensor adjacency matrix using neighbor intersections"
  (let* ((n (first (cl-ana.tensor:tensor-dimensions adj-tensor)))
         (result-tensor (cl-ana.tensor:make-tensor (list n n)
                                                   :initial-element 0.0d0))
         (neighbor-lists (make-array n)))

    ;; Precompute neighbor lists
    (dotimes (i n)
      (setf (aref neighbor-lists i) (neighbors adj-tensor i)))

    ;; Calculate inverse log values
    (dotimes (i n)
      (dotimes (j n)
        (setf (cl-ana.tensor:tref result-tensor i j)
              (inv-log adj-tensor i (aref neighbor-lists j)))))


    ;; Clear diagonal
    (dotimes (i n)
      (setf (cl-ana.tensor:tref result-tensor i i) 0.0d0))


    ;; Normalize matrix
    ;; (dotimes (i n)
    ;;   (let ((row-sum (loop for j below n
    ;;                        sum (cl-ana.tensor:tref result-tensor i j))))
    ;;     (unless (zerop row-sum)
    ;;       (dotimes (j n)
    ;;         (setf (cl-ana.tensor:tref result-tensor i j)
    ;;               (/ (cl-ana.tensor:tref result-tensor i j) row-sum))))))

    ;; Apply cumulative sum if requested
    (normalize-tensor (tensor-to-list result-tensor))))



;;(defun action6 (adj-list &key (cumulative t))
;;  "Process adjacency list using neighbor intersections and inverse log values"
;;  (let* ((n (length adj-list))
;;         (result-matrix (make-array (list n n) :element-type 'double-float
;;                                               :initial-element 0.0d0))
;;         (neighbor-lists (loop for i below n collect (neighbors adj-list i))))
;;
;;    ;; Calculate inverse log values
;;    (dotimes (i n)
;;      (dotimes (j n)
;;        (setf (aref result-matrix i j)
;;              (inv-log adj-list (list i (nth j neighbor-lists))))))
;;
;;    ;; Clear diagonal
;;    (dotimes (i n)
;;      (setf (aref result-matrix i i) 0.0d0))
;;
;;    ;; Normalize rows and handle NaN
;;    (dotimes (i n)
;;      (let ((row-sum (loop for j below n sum (aref result-matrix i j))))
;;        (unless (zerop row-sum)
;;          (dotimes (j n)
;;            (setf (aref result-matrix i j)
;;                  (/ (aref result-matrix i j) row-sum))))))
;;
;;    ;; Apply cumulative sum if requested
;;    (when cumulative
;;      (dotimes (i n)
;;        (let ((cumsum 0.0d0))
;;          (dotimes (j n)
;;            (incf cumsum (aref result-matrix i j))
;;            (setf (aref result-matrix i j) cumsum)))))
;;
;;    ;; Convert result back to nested list format
;;    (array-to-list result-matrix)))


(defun action7 (adj-list )
  "Process adjacency list using Jaccard index"
  (let* ((n (length adj-list))
         (result-matrix (make-array (list n n) :element-type 'double-float
                                               :initial-element 0.0d0))
         ;; (neighbor-lists (loop for i below n
         ;;                       collect (neighbors adj-list i)))
         )

    ;; Calculate Jaccard indices
    (dotimes (i n)
      (dotimes (j n)
        (setf (aref result-matrix i j)
              (jacc (list (neighbors adj-list i)
                          (neighbors adj-list j))))))

    ;; Clear diagonal
    (dotimes (i n)
      (setf (aref result-matrix i i) 0.0d0))

    ;; Normalize rows and handle NaN
    (dotimes (i n)
      (let ((row-sum (loop for j below n sum (aref result-matrix i j))))
        (unless (zerop row-sum)
          (dotimes (j n)
            (setf (aref result-matrix i j)
                  (/ (aref result-matrix i j) (* row-sum n))
                  )))))

    ;; Apply cumulative sum if requested
    ;;(when cumulative
    ;;  (dotimes (i n)
    ;;    (let ((cumsum 0.0d0))
    ;;      (dotimes (j n)
    ;;        (incf cumsum (aref result-matrix i j))
    ;;        (setf (aref result-matrix i j) cumsum)))))

    ;; Convert result to nested list format
    ;;(loop for i below n
    ;;      collect (loop for j below n
    ;;                    collect (coerce (aref result-matrix i j) 'double-float)))
    (list-to-tensor (cl-ana.linear-algebra:lisp-2d-array->tensor result-matrix))
    ))



;; Macros to modify actions to accept thresholds.

;; Testing code
(defparameter *test*  '((0 1 1 1)
                        (1 0 1 0)
                        (1 1 0 0)
                        (1 0 0 0)))

(defparameter *test2* (list-to-tensor *test2*))
(print *test2*)

(defparameter *test*  '((0 1 0 )
                        (1 0 1 )
                        (0 1 0 )))
(print (adjacency-list-to-matrix (matrix-to-adjacency-list *test*)))
(defparameter *temp* (cl-ana.tensor:tensor-flatten *test*))
(print *temp*)
(print (position-if (lambda (x) (> x 0)) *temp* :from-end nil))
(print (brandes *test2*))
(print (action1 *test*))
(print (action2 *test*))
(print (action3 *test*))
;; Implement action 4 betweenness centrality
(print (action4 *test2*))
(print (action5 *test* :cumulative nil))

;; Apparently action6 has an error somewhere, this needs to be addressed.
(print (action6 *test* :cumulative nil))

(print (action7 *test*))


;; Usage:



;; Write the code to pre-filter the actions by which ones can't be feasible based on the attempted fixed-point theorem
;; and action thresholding.

;; After ruling out infeasible actions and finding thresholds for feasible ones, run the optimization on the remainder
;; with the constraint that the starting graph must still be a fixed point of the generator. Then use other standard
;; optimization techniques as before.

;; Leverage LP/NLP/MINLP problem properties on coerced generators to see how many feasible optimal solutions there are and if our other fixed points fall in that region.
