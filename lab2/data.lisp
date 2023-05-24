;----------------------------------
; копіювання масиву
(defun get-array-copy (arr size) 
    (setf copy (make-array N))
    (loop for i from 0 to (- size 1) do (setf (aref copy i) (aref arr i)))
    copy
)

;----------------------------------
; обчислення сумарної вірогідності кожного інтервалу з комбінації
(defun get-int-sum-probability (arr intervals M N)
    (setf probability 0.0)
    (loop for i from 0 to (- M 1) do 
        (setf probability 
            (fisher-z-p
                (aref arr (aref intervals i))
                (if (= i (- M 1)) (aref arr (- N 1)) (aref arr (- (aref intervals (+ i 1)) 1)))
            )
        )
    )
    probability
)

; здійснює перебір інтервалів
(defun get-possible-intervals (arr N M size-index int-sizes min-int-size)
    (setf cur-max-probability 0.0)
    (setf ind-intervals (make-array M))
    (setf local-int-sizes (get-array-copy int-sizes M))
    (loop
        (when (< (+ size-index 1) (- M 1)) 
            (setf tmp-ind-intervals (get-possible-intervals arr N M (+ size-index 1) local-int-sizes min-int-size))
            (setf cur-probability (get-int-sum-probability arr tmp-ind-intervals M N))
            (when (> cur-probability cur-max-probability) 
                (when (= cur-probability 1.0) (return-from get-possible-intervals tmp-ind-intervals))
                (setf cur-max-probability cur-probability)
                (setf ind-intervals (get-array-copy tmp-ind-intervals M))
            )
        )
        (when (= (+ size-index 1) (- M 1))
            (setf sum-size 0)
            (setf tmp-ind-intervals (make-array M))
            (loop for i from 0 to (- M 1) do
                (setf (aref tmp-ind-intervals i) sum-size)
                (setf sum-size (+ sum-size (aref local-int-sizes i)))
            )
            (setf cur-probability (get-int-sum-probability arr tmp-ind-intervals M N))
            (when (> cur-probability cur-max-probability) 
                (when (= cur-probability 1.0) (return-from get-possible-intervals tmp-ind-intervals))
                (setf cur-max-probability cur-probability)
                (setf ind-intervals (get-array-copy tmp-ind-intervals M))
            )
        )
        (setf (aref local-int-sizes size-index) (- (aref local-int-sizes size-index) 1))
        (setf (aref local-int-sizes (+ size-index 1)) (+ (aref local-int-sizes size-index) 1))
        (when (< (aref local-int-sizes size-index) min-int-size) 
            (return-from get-possible-intervals ind-intervals)
        )
    )
)

; повертає інтервали після перебору можливих
(defun divide-into-intervals (arr N M min-int-size) 
    (setf int-sizes (make-array M))
    (loop for i from 0 to (- M 1) do 
        (if (= i 0) 
            (setf (aref int-sizes i) (- N (* min-int-size (- M 1)))) 
            (setf (aref int-sizes i) min-int-size)
        )
    )
    (get-possible-intervals arr N M 0 int-sizes min-int-size)
)

; виокремлення інтервалів з індексів інтервалів
(defun get-intervals (arr intervals-ind N M) 
    (setf intervals (make-array M))
    (loop for i from 0 to (- M 1) do
        (setf tmp (make-array '(2)))
        (set (aref tmp 0)(aref arr (aref intervals i)))
        (setf (aref tmp 1) (if (= i (- M 1)) (aref arr (- N 1)) (aref arr (- (aref intervals (+ i 1)) 1))))
        (setf (aref intervals i) tmp)
    )
    intervals
)

;----------------------------------
; переведення числового ряду на лінгвістичний рядок
(defun translate-to-ling-chain (arr alphabet intervals N M)
    (setf ling-chain (make-array N))
    (loop for i from 0 to (- N 1) do
        (loop for j from 0 to (- M  1) do
            (when (and (<= (aref intervals j 0) (aref arr i)) (<= (aref arr i) (aref intervals j 1))))
        )
    )
    ling-chain
)

;----------------------------------
; виведення масиву
(defun print-array (arr size) (loop for i from 0 to (- size 1) do (write (aref arr i)) (princ " ")))

; виведення лінгвістичного рядка, що використовую виведення масиву
(defun print-ling-chain (arr size) (print "Linguistic chain") print-array(arr size))

;----------------------------------
; головна процедура
(defun main ()
    ; необхідні дані для виконання програми
    (defconstant M 4)
    (defconstant N 14)
    (defconstant MIN_INTERVAL_SIZE 14)
    (setf alphabet (make-array M :initial-contents '("A" "B" "C" "D")))
    (setf numbers (make-array N :initial-contents '(13 1 2 15 -10 20 3 12 4 5 6 26 7 8)))  
    
    (print-ling-chain ( 
        (translate-to-ling-chain arr alphabet 
            (get-intervals numbers 
                (divide-into-intervals (sort (get-array-copy numbers N) #'< ) N M)
            N M) 
        N M)
    N))
)

; вимірювання часу виконання програми
(time (main))

