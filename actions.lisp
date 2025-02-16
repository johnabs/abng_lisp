(defun scal (x y) (aops:each (alexandria:curry #'* x) y))

(defun action_identity (adjacency-matrix)
  (normalize-array adjacency-matrix) )


(defun action1 (adjacency-matrix)
  "Compute the degree-distribution of each node in a graph from its adjacency matrix.

   Parameters:
   adjacency-matrix - a square matrix representing the adjacency matrix of the graph

   Returns:
   a vector where the index corresponds to the node and the value is the degree of that node divided by the total degree of all nodes"
  (let*
      (
       (nodes (array-dimension adjacency-matrix 0))
       (deg (degree adjacency-matrix))
       (d2 (aops:each #'nan-div (lla:mm adjacency-matrix deg) deg))
       (norm (aops:each #'nan-div d2 (make-array nodes :initial-element (lla:asum d2))))
       (m (make-array (list nodes nodes) :initial-element 0.0)))
    ;; Repeat norm vector to create matrix
    (dotimes (i nodes)
      (dotimes (j nodes)
        (setf (aref m i j) (aref norm j))))

    ;; Zero out diagonal
    (dotimes (i nodes)
      (setf (aref m i i) 0.0))

    ;; Normalize rows then rows by number of nodes
    (dotimes (i nodes)
      (let ((row-sum (loop for j below nodes sum (aref m i j))))
        (when (> row-sum 0)
          (dotimes (j nodes)
            (setf (aref m i j) (/ (/ (aref m i j) row-sum) nodes))))))

    ;; Replace NaN with 0 (handled implicitly by the division check above)
    m
    )
  )




(defun action2 (adjacency-matrix)
  "Compute the degree-distribution of each node in a graph from its adjacency matrix.
   Parameters:
   adjacency-matrix - a square matrix representing the adjacency matrix of the graph
   Returns:
   a vector where the index corresponds to the node and the value is the degree of that node divided by the total degree of all nodes"
  (let*
      (
       (nodes (array-dimension adjacency-matrix 0))
       (deg (degree adjacency-matrix))
       (norm (aops:each #'nan-div deg (make-array nodes :initial-element (lla:asum deg))))
       (m (make-array (list nodes nodes) :initial-element 0.0)))
    ;; Repeat norm vector to create matrix
    (dotimes (i nodes)
      (dotimes (j nodes)
        (setf (aref m i j) (aref norm j))))

    ;; Zero out diagonal
    (dotimes (i nodes)
      (setf (aref m i i) 0.0))

    ;; Normalize rows then rows by number of nodes
    (dotimes (i nodes)
      (let ((row-sum (loop for j below nodes sum (aref m i j))))
        (when (> row-sum 0)
          (dotimes (j nodes)
            (setf (aref m i j) (/ (/ (aref m i j) row-sum) nodes))))))

    ;; Replace NaN with 0 (handled implicitly by the division check above)
    m
    )
  )


(defun action3 (adjacency-matrix)
  "Compute the degree-distribution of each node in a graph from its adjacency matrix.
   Parameters:
   adjacency-matrix - a square matrix representing the adjacency matrix of the graph
   Returns:
   a vector where the index corresponds to the node and the value is the degree of that node divided by the total degree of all nodes"
  (let*
      (
       (nodes (array-dimension adjacency-matrix 0))
       (pr (page-rank adjacency-matrix))
       (norm (aops:each #'nan-div pr  (make-array nodes :initial-element (lla:asum pr)) ))
       (m (make-array (list nodes nodes) :initial-element 0.0)))
    ;; Repeat norm vector to create matrix
    (dotimes (i nodes)
      (dotimes (j nodes)
        (setf (aref m i j) (aref norm j))))

    ;; Zero out diagonal
    (dotimes (i nodes)
      (setf (aref m i i) 0.0))

    ;; Normalize rows then rows by number of nodes
    (dotimes (i nodes)
      (let ((row-sum (loop for j below nodes sum (aref m i j))))
        (when (> row-sum 0)
          (dotimes (j nodes)
            (setf (aref m i j) (/ (/ (aref m i j) row-sum) nodes))))))

    ;; Replace NaN with 0 (handled implicitly by the division check above)
    m
    )
  )

(defun action4 (g)
  "Compute normalized betweenness-based matrix from adjacency matrix g"
  (let* ((n (array-dimension g 0))
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
                  (nan-div (aref m i j) matrix-sum))))))

    m))

(defun action5 (adj)
  (let* ((n (array-dimension adj 0))
         (a2 (lla:mm adj adj))
         (a3 (make-array (list n n) :initial-element 0.0))
         (sec-nei nil))

    ;; Set 1's where a2 > adj
    (dotimes (i n)
      (dotimes (j n)
        (when (and (> (aref a2 i j) 0)
                   (= (aref adj i j) 0))
          (setf (aref a3 i j) 1.0))))

    ;; Clear diagonal
    (dotimes (i n)
      (setf (aref a3 i i) 0.0))

    ;; Normalize by row sums
    (let ((row-sums (make-array n :initial-element 0.0)))
      (dotimes (i n)
        (dotimes (j n)
          (incf (aref row-sums i) (aref a3 i j))))

      (setf sec-nei (make-array (list n n) :initial-element 0.0))
      (dotimes (i n)
        (dotimes (j n)
          (if (> (aref row-sums i) 0)
              (setf (aref sec-nei i j)
                    (/ (aref a3 i j) (aref row-sums i)))
              (setf (aref sec-nei i j) 0.0)))))

    ;; Normalize by matrix size
    (dotimes (i n)
      (dotimes (j n)
        (setf (aref sec-nei i j) (/ (aref sec-nei i j) n))))

    ;; Final normalization
    (let ((total-sum (loop for i below n
                           sum (loop for j below n
                                     sum (aref sec-nei i j)))))
      (when (> total-sum 0)
        (dotimes (i n)
          (dotimes (j n)
            (setf (aref sec-nei i j)
                  (/ (aref sec-nei i j) total-sum))))))

    sec-nei))

(defun action6 (adj)
  "Process cl-ana tensor adjacency matrix using neighbor intersections"
  (let* ((n (array-dimension adj 0))
         (result-tensor (make-array (list n n)
                                    :initial-element 0.0d0))
         (neighbor-lists (make-array n)))

    ;; Precompute neighbor lists
    (dotimes (i n)
      (setf (aref neighbor-lists i) (neighbors adj i)))

    ;; Calculate inverse log values
    (dotimes (i n)
      (dotimes (j n)
        (setf (aref result-tensor i j)
              (inv-log adj i (aref neighbor-lists j)))))


    ;; Clear diagonal
    (dotimes (i n)
      (setf (aref result-tensor i i) 0.0d0))


    ;; Normalize matrix
    ;; (dotimes (i n)
    ;;   (let ((row-sum (loop for j below n
    ;;                        sum (cl-ana.tensor:tref result-tensor i j))))
    ;;     (unless (zerop row-sum)
    ;;       (dotimes (j n)
    ;;         (setf (cl-ana.tensor:tref result-tensor i j)
    ;;               (/ (cl-ana.tensor:tref result-tensor i j) row-sum))))))

    ;; Apply cumulative sum if requested
    (normalize-array result-tensor)))

(defun action7 (adj)
  "Process adjacency list using Jaccard index"
  (let* ((n (array-dimension adj 0))
         (result-matrix (make-array (list n n) :element-type 'double-float
                                               :initial-element 0.0d0))
         ;; (neighbor-lists (loop for i below n
         ;;                       collect (neighbors adj-list i)))
         )

    ;; Calculate Jaccard indices
    (dotimes (i n)
      (dotimes (j n)
        (setf (aref result-matrix i j)
              (jacc (list (neighbors adj i)
                          (neighbors adj j))))))

    ;; Clear diagonal
    (dotimes (i n)
      (setf (aref result-matrix i i) 0.0d0))

    (normalize-array result-matrix)))
