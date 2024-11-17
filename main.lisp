;; Enable maximum debug settings, can change to space and speed for better performance later
(declaim (optimize (speed 0) (space 0) (debug 3)))

;; Load relevant dependencies
                                        ;(ql:quickload :cl-graph)
(ql:quickload :cl-ana)
(ql:quickload :distributions)
(ql:quickload :alexandria)


;; Helper Functions
(defun ones (l) (cl-ana.linear-algebra:make-vector l :initial-element 1))
(defun nan-div (x y) (if (equalp 0 y) 0 (/ x y))
  )
;; Start with a simple reimplementation of ABNG code
;;



;; Action Statistical Properties Functions
(defun vecsum (mat)  (cl-ana.linear-algebra:euclidean-dot mat (ones (length mat))))
(defun degree (adjacency-matrix)
  "Compute the degree of each node in a graph from its adjacency matrix.

   Parameters:
   adjacency-matrix - a square matrix representing the adjacency matrix of the graph

   Returns:
   a vector where the index corresponds to the node and the value is the degree of that node"
  (cl-ana.linear-algebra:euclidean-dot  (ones (length adjacency-matrix)) adjacency-matrix)
  )

(defun normalize-matrix (matrix)
  "Normalize the columns of the adjacency matrix"
  (let* ((n (length matrix))
         (col-sums (make-array n :initial-element 0))
         (normalized (make-array (list n n) :initial-element 0)))

    ;; Calculate column sums
    (loop for i from 0 below n do
      (loop for j from 0 below n do
        (incf (aref col-sums j) (nth i (nth j matrix)))))

    ;; Normalize each column
    (loop for i from 0 below n do
      (loop for j from 0 below n do
        (if (> (aref col-sums j) 0)
            (setf (aref normalized i j)
                  (/ (float (nth i (nth j matrix)))
                     (float (aref col-sums j)))))))
    normalized))

(defun page-rank (adjacency-matrix &key (damping 0.85) (iterations 100) (epsilon 1e-8))
  "Compute PageRank values for nodes in a graph.

   Parameters:
   adjacency-matrix - adjacency matrix representation of the graph
   damping - damping factor (default: 0.85)
   iterations - maximum number of iterations (default: 100)
   epsilon - convergence threshold (default: 1e-8)

   Returns:
   Vector of PageRank values for each node"
  (let* ((n (length adjacency-matrix))
         (normalized-matrix (normalize-matrix adjacency-matrix))
         (rank-vector (make-array n :initial-element (/ 1.0 n)))
         (teleport-vector (make-array n :initial-element (/ 1.0 n))))

    ;; Power iteration method
    (dotimes (iter iterations)
      (let ((new-rank (make-array n :initial-element 0.0))
            (diff 0.0))

        ;; Calculate new rank values
        (loop for i from 0 below n do
          (setf (aref new-rank i)
                (+ (* (- 1 damping) (aref teleport-vector i))
                   (* damping
                      (loop for j from 0 below n
                            sum (* (aref normalized-matrix i j)
                                   (aref rank-vector j))))))

          ;; Calculate difference for convergence check
          (incf diff (abs (- (aref new-rank i) (aref rank-vector i)))))

        ;; Update rank vector
        (loop for i from 0 below n do
          (setf (aref rank-vector i) (aref new-rank i)))

        ;; Check for convergence
        (when (< diff epsilon)
          (return-from page-rank rank-vector))))

    rank-vector))

;; Example usage
(let ((adjacency-matrix '((0 1 0 0 0 0)
                          (1 0 1 0 0 0 )
                          (0 1 0 1 0 0 )
                          (0 0 1 0 1 0 )
                          (0 0 0 1 0 1 )
                          (0 0 0 0 1 0 )
                          )
                        ))
  (action3 adjacency-matrix))



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



;; Macros to modify actions to accept thresholds.

;; Testing code
(defparameter *test*  '((0 1 0 1)
                        (1 0 0 0)
                        (0 0 0 0)
                        (1 0 0 0)))
(print (action1 *test*))
(print (action2 *test*))


;; Write the code to pre-filter the actions by which ones can't be feasible based on the attempted fixed-point theorem
;; and action thresholding.


;; After ruling out infeasible actions and finding thresholds for feasible ones, run the optimization on the remainder
;; with the constraint that the starting graph must still be a fixed point of the generator. Then use other standard
;; optimization techniques as before.


;; Question to ask Mario: feedback between


;; Confirm


;; Leverage LP/NLP/MINLP problem properties on coerced generators to see how many feasible optimal solutions there are and if our other fixed points fall in that region.
