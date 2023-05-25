; main.lisp

;----------------------------------
; допоміжні функції для обчислення функції розподілу

; допоміжна функція, що округляє float до n знаків після коми - необхідно проти помилки floating point undeflow
(defun cut-float-digits (a digits)
    (let (
        float-d
    )
        (setq float-d (expt 10 digits))
        (setq a (* a float-d))
        (setq a (fround a))
        (setq a (/ a float-d))
        a
    )
)

; допоміжна функція для обчислення непоіної бета функції
(defun inc-beta-fraction (a b x)
    (let (
        (FLOAT-DIGITS 38)

        (MAXIT 200)
        (EPS 3.0e-7)
        (FPMIN 1.0e-30)
        m m1
        h aa del
        qab qap qam
        c d
    )
        (setq qab (+ a b))
        (setq qap (+ a 1.0))
        (setq qam (- a 1.0))
        (setq c 1.0)
        (setq d (- 1.0 (/ (* qab x) qap)))
        (when (< (abs d) FPMIN) (setq d FPMIN))
        (setq d (/ 1.0 d))
        (setq h d)

        (loop for i from 1 to MAXIT do 
            (setq m i)
            (setq m2 (* 2.0 m))

            (setq aa (/ (* m (* (- b m) x))
                (* (+ qam m2) (+ a m2))))
            (setq aa (cut-float-digits aa FLOAT-DIGITS))

            (setq d (+ 1.0 (* aa d)))
            (setq d (cut-float-digits d FLOAT-DIGITS))
            (when (< (abs d) FPMIN) (setf d FPMIN))

            (setq c (+ 1.0 (/ aa c)))
            (setq c (cut-float-digits c FLOAT-DIGITS))
            (when (< (abs c) FPMIN) (setf c FPMIN))

            (setq d (/ 1.0 d))
            (setq d (cut-float-digits d FLOAT-DIGITS))

            (setq h (* h (* d c)))
            (setq h (cut-float-digits h FLOAT-DIGITS))

            (setq aa (/ 
                        (* (* (* (- 0 1) (+ a m)) (+ qab m)) x) 
                        (* (+ a m2) (+ qap m2))
                    )
            )
            (setq aa (cut-float-digits aa FLOAT-DIGITS))

            (setq d (+ 1.0 (* aa d)))
            (setq d (cut-float-digits d FLOAT-DIGITS))
            (when (< (abs d) FPMIN) (setf d FPMIN))

            (setq c (+ 1.0 (/ aa c)))
            (setq c (cut-float-digits c FLOAT-DIGITS))
            (when (< (abs c) FPMIN) (setf c FPMIN))

            (setq d (/ 1.0 d))
            (setq d (cut-float-digits d FLOAT-DIGITS))

            (setq del (* d c))
            (setq del (cut-float-digits del FLOAT-DIGITS))

            (setq h (* h del))
            (setq h (cut-float-digits h FLOAT-DIGITS))

            (when (< (abs (- del 1.0)) EPS) (return))
        )
        (if (> m MAXIT) Nil h)
    )
)

; обчислення неповної бета функції для обчислення I для функції розподілу
(defun incomplete-beta (a b x)
    (let (
        beta
    )
        (when (or (< x 0.0) (> x 1.0)) Nil)
        (when (or (= x 0) (= x 1)) (setq beta 0.0))
        (setq beta 
            (exp (+ (+ (- (- (lgamma (+ a b)) (lgamma a)) (lgamma b)) (* a (log x))) (* b (log (- 1 x)))) )
        )

        (if (< x (/ (+ a 1) (+ a (+ b 2)))) 
            (/ (* beta (inc-beta-fraction a b x)) a)
            (- 1 (/ (* beta (inc-beta-fraction a b (- 1 x))) b))
        )
    )
)

; обчислення повної бета функції для обчислення I для функції розподілу
(defun complete-beta (a b) 
    (exp (- 
        (+ (lgamma a) (lgamma b))  (lgamma (+ a b))
    ))
)

; обчислення x*
(defun count-xi (z d1 d2) 
    ( /
        (* d2 (exp (* 2 z))) 
        (+ d1 (* d2 (exp (* 2 z))))
    )
)

; обчислення z = (x - mu) / sigma
(defun count-z (x mu sigma) (/ (- x mu) sigma))

;----------------------------------
; обчислення функції розподілу
; F = Ix*
(defun fisher-z-fn (x) 
    (let (
        (MU 6.0)
        (SIGMA 7.0)
        (D1 3.0)
        (D2 2.0)
    )
        (/ 
            (incomplete-beta (/ D1 2) (/ D2 2) (count-xi (count-z x MU SIGMA) D1 D2)) 
            (complete-beta (/ D1 2) (/ D2 2))
        )
    )
)

; P(a, b) = F(b) - F(a)
(defun fisher-z-p (a b) (- (fisher-z-fn b) (fisher-z-fn a)) )

;----------------------------------
; функції для роботи з числовим рядом

;----------------------------------
; копіювання масиву
(defun get-array-copy (arr size) 
    (let (
        copy
    )
        (setq copy (make-array size))
        (loop for i from 0 to (- size 1) do (setf (aref copy i) (aref arr i)))
        copy
    )
)

;----------------------------------
; вибір підходящої комбінації інтервалів
(defun get-probable-intervals-from-list (arr intervals N M)
    (let (
        max-probability probable-intervals
    )
        (setq max-probability 0.0)
        (setq probable-intervals (make-array M))
        (loop for interval-indexes in intervals do
            (let (
                (cur-probability 0.0) 
                tmp-arr
            )
                (setq tmp-arr (make-array M))
                (loop for i from 0 to (- M 1) do
                    (let (
                        a b tmp-interval-arr
                    )
                        (setq tmp-interval-arr (make-array '(2)))

                        (setq a (aref arr (aref interval-indexes i)))
                        (setq b (if (equal i (- M 1)) (aref arr (- N 1)) (aref arr (- (aref interval-indexes (+ i 1)) 1))))

                        (setf (aref tmp-interval-arr 0) a)
                        (setf (aref tmp-interval-arr 1) b)
                        (setf (aref tmp-arr i) tmp-interval-arr)

                        (setq cur-probability (+ cur-probability (fisher-z-p a b)))
                    )
                )
                (when (= cur-probability 1.0) (return-from interval-sum-probability tmp-arr))
                (when (> cur-probability max-probability)
                    (setq probable-intervals tmp-arr)
                    (setq max-probability cur-probability)
                )
            )
        )
        probable-intervals
    )
)

; допоміжна функція для перебору інтервалів
(defun get-cur-interval-sizes (interval-sizes M size-index min-interval-size)
    (let (
        interval-size-list
        cur-interval-sizes
    )
        (setq interval-size-list (list))
        (setq cur-interval-sizes (get-array-copy interval-sizes M))
        (loop 
            (let (
                tmp-interval-sizes
                tmp-sizes-list
            )
                (setq tmp-interval-sizes (get-array-copy cur-interval-sizes M))
                (setq tmp-sizes-list (list Nil))
                (setf (car tmp-sizes-list) tmp-interval-sizes)
                (setq interval-size-list (append interval-size-list tmp-sizes-list))

                (decf (aref cur-interval-sizes size-index))
                (incf (aref cur-interval-sizes (+ size-index 1)))
                (when (< (aref cur-interval-sizes size-index) min-interval-size) (return))
            )
        )
        interval-size-list
    )
)

; здійснює перебір інтервалів
(defun get-possible-interval-sizes-list (interval-sizes M min-interval-size)
    (let (
        (size-index 0)
        cur-interval-sizes
        interval-size-list

    )
        (setq cur-interval-sizes (get-array-copy interval-sizes M))
        (setq interval-size-list (get-cur-interval-sizes cur-interval-sizes M size-index MIN-INTERVAL-SIZE))
        (loop
            (let (
                tmp-size-list
            )
                (setq tmp-size-list (list))
                (incf size-index)
                (loop for arr in interval-size-list do 
                    (setf tmp-size-list (append tmp-size-list (get-cur-interval-sizes arr M size-index MIN-INTERVAL-SIZE)))
                )
                (setq interval-size-list tmp-size-list)
                (when (= size-index (- M 2)) (return))
            )
        )
        interval-size-list
    )
)

; повертає розмірності інтервалів після перебору можливих
(defun divide-into-interval-sizes (N M min-interval-size) 
    (let (
        interval-sizes
    )
        (setq interval-sizes (make-array M))
        (loop for i from 0 to (- M 1) do 
            (if (= i 0) 
                (setf (aref interval-sizes i) (- N (* min-interval-size (- M 1)))) 
                (setf (aref interval-sizes i) min-interval-size)
            )
        )
        (get-possible-interval-sizes-list interval-sizes M MIN-INTERVAL-SIZE)
    )
)

; здійснює перетворення розмірностей інтервалів на індекси інтервалів
(defun transform-interval-size-to-index-list (interval-sizes M)
    (let (
        index-list
    )
        (setq index-list (list))
        (loop for cur-interval-sizes in interval-sizes do
            ( let (
                (sum-interval-size 0)
                tmp-list
                tmp-index-intervals
            )
                (setq tmp-list (list Nil))
                (setq sum-interval-size 0)
                (setq tmp-index-intervals (make-array M))
                (loop for i from 0 to (- M 1) do
                    (setf (aref tmp-index-intervals i) sum-interval-size)
                    (setq sum-interval-size (+ sum-interval-size (aref cur-interval-sizes i)))
                )
                (setf (car tmp-list) tmp-index-intervals)
                (setq index-list (append index-list tmp-list))
            )
        )
        index-list
    )
)

; повертає перетворений список розмірностей на індекси інтервалів
(defun get-possible-interval-index-list (N M min-interval-size)
    (transform-interval-size-to-index-list (divide-into-interval-sizes N M MIN-INTERVAL-SIZE) M)
)

; повертає вірогідний інтервал відповідно до функції розподілу
(defun get-probable-intervals (arr N M min-interval-size)
    (get-probable-intervals-from-list arr (get-possible-interval-index-list N M min-interval-size) N M)
)

;----------------------------------
; переведення числа на літеру
(defun translate-number (arr alphabet intervals i M)
    (loop for j from 0 to (- M 1) do 
        (when (and 
                (<= (aref (aref intervals j) 0) (aref arr i)) 
                (<= (aref arr i) (aref (aref intervals j) 1))
            )
            (return-from translate-number (aref alphabet j))
        )
    )
)

; переведення числового ряду на лінгвістичний рядок
(defun translate-to-ling-chain (arr alphabet intervals N M)
    (let (
        ling-chain
    )
        (setq ling-chain (make-array N))
        (loop for i from 0 to (- N 1) do
            (setf (aref ling-chain i) (translate-number arr alphabet intervals i M))
        )
        ling-chain
    )
    
)

;----------------------------------
; виведення масиву
(defun print-array (arr size) 
    (loop for i from 0 to (- size 1) do 
        (format t (aref arr i)) 
        (princ " ") 
    )
)

; виведення лінгвістичного рядка, що використовую виведення масиву
(defun print-ling-chain (arr size) 
    (format t "Linguistic chain~%") 
    (print-array arr size)
)

;----------------------------------
; головна процедура
(defun main ()
    (let (
        ; необхідні дані для виконання програми
        (M 4)
        (N 14)
        (MIN-INTERVAL-SIZE 2)
        alphabet numbers
    ) 
        (setq alphabet (make-array M :initial-contents '("A" "B" "C" "D")))
        (setq numbers (make-array N :initial-contents '(13 1 2 15 -10 20 3 12 4 5 6 26 7 8))) 
        
        (setq intervals (get-probable-intervals (sort (get-array-copy numbers N) #'< ) N M MIN-INTERVAL-SIZE))
        (print-ling-chain (
            translate-to-ling-chain numbers alphabet intervals N M
        ) N)
    )
)

; вимірювання часу виконання програми
 (time (main))



