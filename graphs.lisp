
(ql:quickload :cl-ana)
(ql:quickload :distributions)
(ql:quickload :alexandria)

(defun ones (l) (cl-ana.linear-algebra:make-vector l :initial-element 1))
(defun nan-div (x y) (if (equalp 0 y) 0 (/ x y)
                         )
  )

;; Start with a simple reimplementation of ABNG code
(defun list-to-tensor (nested-list)
  (let* ((n (length nested-list))
         (tensor (cl-ana.tensor:make-tensor (list n n) :initial-element 0)))
    (dotimes (i n)
      (dotimes (j n)
        (setf (cl-ana.tensor:tref tensor i j)
              (nth j (nth i nested-list)))))
    tensor))

(defun tensor-to-list (tensor)
  (let* ((dims (cl-ana.tensor:tensor-dimensions tensor))
         (rows (first dims))
         (cols (second dims)))
    (loop for i below rows
          collect (loop for j below cols
                        collect (cl-ana.tensor:tref tensor i j)))))
(defun list-to-array (lst)
  "Convert nested list to 2D array"
  (let* ((n (length lst))
         (arr (make-array (list n n))))
    (dotimes (i n)
      (dotimes (j n)
        (setf (aref arr i j) (nth j (nth i lst)))))
    arr))

(defun col-sums (matrix)
  "Calculate sum of each row in matrix"
  (cl-ana.linear-algebra:euclidean-dot matrix (ones (length matrix)))
  )

(defun row-sums (matrix)
  "Calculate sum of each row in matrix"
  (cl-ana.linear-algebra:euclidean-dot (ones (length matrix)) (cl-ana.linear-algebra:matrix-transpose matrix))
  )

(defun eachrow (matrix)
  "Calculate sum of each row in matrix"
  (mapcar (lambda (x) (cl-ana.tensor:tref matrix x)) (alexandria:iota (length matrix)))
  )

(defun eachcol (matrix)
  "Calculate sum of each row in matrix"
  (eachrow (cl-ana.linear-algebra:matrix-transpose matrix))
  )

(defun array-to-list (arr)
  "Convert 2D array to nested list"
  (let ((n (array-dimension arr 0)))
    (loop for i below n
          collect (loop for j below n
                        collect (coerce (aref arr i j) 'double-float)))))

(defun degree (adj-list vertex)
  "Calculate degree of a vertex in the graph"
  (loop for x in (nth vertex adj-list)
        sum x))

;; (defun neighbors (adj-list vertex)
;;   "Get list of neighbors for a vertex"
;;   (loop for i from 0
;;         for val in (nth vertex adj-list)
;;         when (= 1 val)
;;           collect i))

(defun ith-degree (adj-tensor vertex)
  "Calculate degree of the ith vertex in the tensor graph"
  (let ((n (first (cl-ana.tensor:tensor-dimensions adj-tensor))))
    (loop for i below n
          sum (cl-ana.tensor:tref adj-tensor vertex i))))

(defun vector-intersection (vec1 vec2)
  "Calculate intersection of two vectors"
  (let ((result '()))
    (dotimes (i (length vec1))
      (when (find (aref vec1 i) vec2)
        (push (aref vec1 i) result)))
    (coerce (nreverse result) 'vector)))

(defun neighbors (adj-tensor vertex)
  "Get neighbors for a vertex from tensor adjacency matrix"
  (let* ((n (first (cl-ana.tensor:tensor-dimensions adj-tensor)))
         (neighbors-list '()))
    (dotimes (i n)
      (when (= 1 (cl-ana.tensor:tref adj-tensor vertex i))
        (push i neighbors-list)))
    (coerce (nreverse neighbors-list) 'vector)))


(defun g0p (x)  (if (> x 0) 1 0))

(defun subtract-matrices (mat1 mat2)
  "Subtract two matrices represented as nested lists"
  (let* ((n (length mat1))
         (result (make-list n :initial-element
                            (make-list n :initial-element 0))))
    (dotimes (i n)
      (dotimes (j n)
        (setf (nth j (nth i result))
              (- (nth j (nth i mat1))
                 (nth j (nth i mat2))))))
    result))

(defun divide-by-row-sums (matrix row-sums)
  "Divide each matrix element by its row sum, replace NaN with 0"
  (let* ((n (length matrix))
         (result (list-to-tensor (make-list n :initial-element
                                            (make-list n :initial-element 0.0)))))
    (dotimes (i n)
      (let ((row-sum (aref row-sums i)))
        (dotimes (j n)
          (setf (aref (aref result i) j)
                (if (zerop row-sum)
                    0.0
                    (/ (aref (aref matrix i) j)
                       row-sum))))))
    result))

(defun cumulative-sum-rows (matrix)
  "Calculate cumulative sum along each row"
  (let* ((n (length matrix))
         (result (make-list n :initial-element
                            (make-list n :initial-element 0.0))))
    (dotimes (i n)
      (let ((sum 0.0))
        (dotimes (j n)
          (incf sum (nth j (nth i matrix)))
          (setf (nth j (nth i result)) sum))))
    result))
(defun matrix-multiply (mat1 mat2)
  "Multiply two matrices represented as nested lists"
  (let* ((n (length mat1))
         (result (make-list n :initial-element
                            (make-list n :initial-element 0))))
    (dotimes (i n)
      (dotimes (j n)
        (setf (nth j (nth i result))
              (loop for k below n
                    sum (* (nth k (nth i mat1))
                           (nth j (nth k mat2)))))))
    result))


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

(defun inv-log (adj-tensor vertex neighbors-vec)
  "Calculate inverse log based on intersection of neighbor sets"
  (let* ((vertex-neighbors (neighbors adj-tensor vertex))
         (common (vector-intersection vertex-neighbors neighbors-vec)))
    (if (plusp (length common))
        (loop for i below (length common)
              sum (/ 1.0d0 (log (1+ (ith-degree adj-tensor (aref common i))))))
        0.0d0)))


;;(defun inv-log (adj-list pair)
;;  "Calculate inverse log based on intersection of neighbor sets"
;;  (destructuring-bind (vertices neighbors-list) pair
;;    (let* ((vertex-neighbors (neighbors adj-list vertices))
;;           (common (intersection vertex-neighbors neighbors-list)))
;;      (if common
;;          (loop for vertex in common
;;                sum (/ 1.0d0 (log (1+ (aref (degree adj-list) vertex)))))
;;          0.0d0))))


(defun vector-union (vec1 vec2)
  "Calculate union of two sorted vectors"
  (let ((result '())
        (i 0)
        (j 0))
    (loop while (and (< i (length vec1))
                     (< j (length vec2)))
          do (cond ((= (aref vec1 i) (aref vec2 j))
                    (push (aref vec1 i) result)
                    (incf i)
                    (incf j))
                   ((< (aref vec1 i) (aref vec2 j))
                    (push (aref vec1 i) result)
                    (incf i))
                   (t
                    (push (aref vec2 j) result)
                    (incf j))))
    ;; Add remaining elements
    (loop while (< i (length vec1))
          do (push (aref vec1 i) result)
             (incf i))
    (loop while (< j (length vec2))
          do (push (aref vec2 j) result)
             (incf j))
    (coerce (nreverse result) 'vector)))



(defun jacc (pair)
  "Calculate Jaccard index for two sets of neighbors"
  (destructuring-bind (vertices neighbors-list) pair
    (let* ((set1 vertices)
           (set2 neighbors-list)
           (intersection-size (length (vector-intersection set1 set2)))
           (union-size (length (vector-union set1 set2))))
      (if (zerop union-size)
          0.0d0
          (/ intersection-size union-size 1.0d0)))))


(defun get-vertex-count (graph)
  (first (cl-ana.tensor:tensor-dimensions graph)))

;; Helper function to get neighbors of a vertex from cl-ana tensor
(defun get-neighbors (graph vertex)
  (let ((n (get-vertex-count graph))
        (neighbors nil))
    (dotimes (i n)
      (when (not (zerop (cl-ana.tensor:tref graph vertex i)))
        (push i neighbors)))
    (nreverse neighbors)))

;; Queue implementation for BFS
(defun make-queue () nil)

(defun queue-empty-p (stack)
  (equalp stack nil ))

(defmacro enqueue (item queue)
  `(if (queue-empty-p ,queue) (setf ,queue (list ,item)) (nconc ,queue (list ,item)))
  )

(defmacro dequeue (queue)
  `(if (queue-empty-p ,queue) nil (setf ,queue (cdr ,queue)))
  )

;; Stack implementation
(defun make-stack ()
  nil)


(defmacro push-stack (item queue)
  `(if (queue-empty-p ,queue) (setf ,queue (list ,item)) (nconc ,queue (list ,item)))
  )

(defmacro pop-stack (queue)
  `(if (queue-empty-p ,queue) nil (setf ,queue (reverse (cdr (reverse ,queue)))))
  )

(defun push-stack (item stack)
  (cons item stack)
  )

(defun stack-empty-p (stack)
  (equalp stack nil ))


;; Main Brandes algorithm implementation
;;(defun brandes (G)
;;  (let ((n (length G))
;;        (unfound -1)
;;        )
;;
;;
;;
;;    )
;;  )


(defparameter *rv-uniform* (distributions:r-uniform 0 1))
(print (distributions:draw *rv-uniform*))
(defun cumsum (list)
  (let ((newlist '(0)))
    (loop for i in (alexandria:iota (length list)) do
      (push (+ (first newlist) (nth i list)) newlist)
          )
    (cdr (reverse newlist))
    )
  )

(print (cumsum '(1 0 0 0 0)))

(defun vector-sample-n (tensor weights n)
  (let ((u01 (distributions:r-uniform 0 1))
        (dims (cl-ana.tensor:tensor-dimensions tensor))
        (vec (cl-ana.tensor:tensor-flatten tensor))
        (probs (cumsum (cl-ana.tensor:tensor-flatten weights)))
        )
    (loop for i in (alexandria:iota n)
          collect (position-if (lambda (x) (> x (distributions:draw u01))) probs :from-end nil)
          )
    )
  )

(defun matrix-to-adjacency-list (matrix)
  "Convert an adjacency matrix (list of lists) to an adjacency list representation.
   Each element in the result is a cons of (vertex . neighbors) where neighbors is
   a list of vertices that the vertex is connected to."
  (loop for row in matrix
        for vertex from 0
        collect (cons vertex
                      (loop for element in row
                            for neighbor from 0
                            when (= element 1)
                              collect neighbor))))

(defun adjacency-list-to-matrix (adj-list)
  "Convert an adjacency list representation ((vertex . neighbors) ...)
   to an adjacency matrix represented as a list of lists"
  (let* ((vertices (length adj-list))
         ;; Create a new matrix for each row to avoid sharing structure
         (matrix (loop for i below vertices
                       collect (make-list vertices :initial-element 0))))
    ;; For each vertex and its neighbors
    (dolist (vertex-entry adj-list)
      (let* ((vertex (car vertex-entry))
             (neighbors (cdr vertex-entry))
             ;; Create a new row to avoid modifying shared structure
             (new-row (copy-list (nth vertex matrix))))
        ;; For each neighbor, set the corresponding entry to 1
        (dolist (neighbor neighbors)
          (setf (nth neighbor new-row) 1))
        ;; Replace the entire row
        (setf (nth vertex matrix) new-row)))
    matrix))

;; Example usage:
(defun brandes (input &optional (directed nil))
  "Compute betweenness centrality in an unweighted graph.
   vertices is a list of vertices, adj is an adjacency list where each element is (vertex . neighbors)
   directed is optional boolean parameter, defaults to nil (undirected)"
  (let (
        (centrality (make-hash-table))
        (vertices (alexandria:iota (length input)))
        (adj (matrix-to-adjacency-list input))

        )
    ;; Initialize centrality scores to 0
    (dolist (v vertices)
      (setf (gethash v centrality) 0))

    (dolist (source vertices)
      (let ((stack nil)
            (pred (make-hash-table))
            (paths (make-hash-table))
            (dist (make-hash-table))
            (queue nil))

        ;; Initialize data structures
        (dolist (w vertices)
          (setf (gethash w pred) nil)
          (setf (gethash w paths) 0)
          (setf (gethash w dist) -1))

        ;; Set initial values for source
        (setf (gethash source paths) 1)
        (setf (gethash source dist) 0)
        (setf queue (list source))

        ;; BFS phase
        (loop while queue do
          (let ((v (pop queue)))
            (push v stack)
            (dolist (w (cdr (assoc v adj)))
              (when (< (gethash w dist) 0)
                (setf queue (append queue (list w)))
                (setf (gethash w dist) (1+ (gethash v dist))))
              (when (= (gethash w dist) (1+ (gethash v dist)))
                (setf (gethash w paths) (+ (gethash w paths) (gethash v paths)))
                (push v (gethash w pred))))))

        ;; Dependency accumulation
        (let ((depend (make-hash-table)))
          (dolist (v vertices)
            (setf (gethash v depend) 0))

          (loop while stack do
            (let ((w (pop stack)))
              (dolist (v (gethash w pred))
                (let ((term (* (/ (gethash v paths)
                                  (gethash w paths))
                               (1+ (gethash w depend)))))
                  (setf (gethash v depend)
                        (+ (gethash v depend) term))))
              (unless (eq w source)
                (setf (gethash w centrality)
                      (+ (gethash w centrality)
                         (gethash w depend))))))))

      ;; If graph is undirected, divide all centrality scores by 2


      ) (unless directed
      (dolist (v vertices)
        (setf (gethash v centrality)
              (/ (gethash v centrality) 2))))
    (mapcar #'(lambda (v) (gethash v centrality))
            vertices)
    ))
