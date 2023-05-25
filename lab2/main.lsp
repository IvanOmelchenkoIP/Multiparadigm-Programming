;----------------------------------
; допоміжні функції для обчислення функції розподілу

; допоміжна функція, що округляє float до n знаків після коми - необхідно проти помилки floating point undeflow
(defun cut-float-digits (a digits)
    (setf float-d (expt 10 digits))
    (setf a (* a float-d))
    (setf a (fround a))
    (setf a (/ a float-d))
    a
)

; допоміжна функція для обчислення непоіної бета функції
(defun inc-beta-fraction (a b x)
    (setf MAXIT 200)
    (setf EPS 3.0e-7)
    (setf FPMIN 1.0e-30)

    (setf m 1) (setf m1 0)
    (setf h 0.0)
    (setf aa 0.0) (setf del 0.0)
    (setf qab (+ a b))
    (setf qap (+ a 1.0))
    (setf qam (- a 1.0))
    (setf c 1.0)
    (setf d (- 1.0 (/ (* qab x) qap)))
    (when (< (abs d) FPMIN) (setf d FPMIN) )
    (setf d (/ 1.0 d))
    (setf h d)
    (loop for i from 1 to MAXIT do 
        (setf m i)
        (setf m2 (* 2.0 m))

        (setf aa (/ (* m (* (- b m) x))
            (* (+ qam m2) (+ a m2))))
        (setf aa (cut-float-digits aa 38))

        (setf d (+ 1.0 (/ a c)))
        (setf d (cut-float-digits aa 38))

        (when (< (abs d) FPMIN) (setf d FPMIN))

        (setf c (+ 1.0 (/ aa c)))
        (setf c (cut-float-digits aa 38))

        (when (< (abs c) FPMIN) (setf c FPMIN))

        (setf d (/ 1.0 d))
        (setf d (cut-float-digits aa 38))

        (setf h (* h (* d c)))
        (setf h (cut-float-digits aa 38))

        (setf aa (/ (* (* (* -1.0 (+ a m)) (+ qab m)) x) 
            (* (+ a m2) (+ qap m2))))
        (setf aa (cut-float-digits aa 38))

        (setf d (+ 1.0 (/ a c)))
        (setf d (cut-float-digits aa 38))

        (when (< (abs d) FPMIN) (setf d FPMIN))

        (setf c (+ 1.0 (/ aa c)))
        (setf c (cut-float-digits aa 38))

        (when (< (abs c) FPMIN) (setf c FPMIN))

        (setf d (/ 1.0 d))
        (setf d (cut-float-digits aa 38))

        (setf del (* d c))
        (setf del (cut-float-digits aa 38))

        (setf h (* h del))
        (setf h (cut-float-digits aa 38))

        (when (< (abs (- del 1.0)) EPS) return)
    )
    (if (> m MAXIT) Nil h)
)

; обчислення неповної бета функції для обчислення I для функції розподілу
(defun incomplete-beta (a b x)
    (when (or (< x 0.0) (> x 1.0)) Nil)
    (setf beta 0.0)
    (when (or (= x 0) (= x 1)) (setf beta 0.0))
    (setf beta 
        (exp (+ (+ (- (- (lgamma (+ a b)) (lgamma a)) (lgamma b)) (* a (log x))) (* b (log (- 1.0 x)))))
    )
    (if (< x (/ (+ a 1) (+ a (+ b 2)))) 
        (/ (* beta (inc-beta-fraction a b x)) a)
        (- 1 (/ (* beta (inc-beta-fraction a b (- 1 x)) b)))
    )
)

; обчислення повної бета функції для обчислення I для функції розподілу
(defun complete-beta (a b) 
    (exp (- (+ (lgamma a) (lgamma b) ) (lgamma (+ a b))))
)

; обчислення x*
(defun count-xi (z d1 d2) ( /
    (* d2 (exp (* 2 z))) 
    (+ d1 (* d2 (exp (* 2 z))))
    ))


; обчислення z = (x - mu) / sigma
(defun count-z (x mu sigma) (/ (- x mu) sigma))

; F = Ix*
(defun fisher-z-fn (x) 
    (defconstant MU 6.0)
    (defconstant SIGMA 7.0)
    (defconstant D1 3.0)
    (defconstant D2 2.0)

    (/ 
        (incomplete-beta (/ D1 2) (/ D2 2) (count-xi (count-z x MU SIGMA) D1 D2)) 
        (complete-beta (/ D1 2) (/ D2 2))
    )
)

; P(a, b) = F(b) - F(a)
(defun fisher-z-p (a b) (- (fisher-z-fn b) (fisher-z-fn a)) )

;(print  (fisher-z-fn -2) )


;----------------------------------
; копіювання масиву
(defun get-array-copy (arr size) 
    (setf copy (make-array size))
    (loop for i from 0 to (- size 1) do (setf (aref copy i) (aref arr i)))
    copy
)

;----------------------------------
; вибір підходящої комбінації інтервалів
(defun get-probable-intervals-from-list (arr intervals N M)
    (setf max-probability 0.0)
    (setf probable-intervals (make-array M))
    ;(print "arr")
    ;(print intervals)
    (loop for interval-indexes in intervals do
        ;(print interval-indexes)
        (setf a 0)
        (setf b 0)
        (setf cur-probability 0.0)
        (setf tmp-arr (make-array M))
        (loop for i from 0 to (- M 1) do
            (setf tmp-interval-arr (make-array '(2)))

            (setf a (aref arr (aref interval-indexes i)))
            (setf b (if (equal i (- M 1)) (aref arr (- N 1)) (aref arr (- (aref interval-indexes (+ i 1)) 1))))

            (setf (aref tmp-interval-arr 0) a)
            (setf (aref tmp-interval-arr 1) b)
            (setf (aref tmp-arr i) tmp-interval-arr)

            (setf cur-probability (+ cur-probability (fisher-z-p a b)))
        )

        (when (= cur-probability 1.0) (return-from interval-sum-probability tmp-arr))
        (when (> cur-probability max-probability)
            (setf probable-intervals tmp-arr)
            (setf max-probability cur-probability)
        )
    )
    probable-intervals
)

; допоміжна функція для перебору інтервалів
(defun get-cur-interval-sizes (interval-sizes M size-index min-interval-size)
    (setf result (make-array '(2)))
    (setf possible-intervals (list))
    (setf interval-size-list (list))
    (setf cur-interval-sizes (get-array-copy interval-sizes M))
    (loop 
        (setf tmp-interval-sizes (get-array-copy cur-interval-sizes M))
        (setf tmp-sizes-list (list Nil))
        (setf (car tmp-sizes-list) tmp-interval-sizes)
        (setf interval-size-list (append interval-size-list tmp-sizes-list))

        (decf (aref cur-interval-sizes size-index))
        (incf (aref cur-interval-sizes (+ size-index 1)))
        (when (< (aref cur-interval-sizes size-index) min-interval-size) (return))
    )
    interval-size-list
)

; здійснює перебір інтервалів
(defun get-possible-interval-sizes-list (interval-sizes M min-interval-size)
    (setf size-index 0)
    (setf cur-interval-sizes (get-array-copy interval-sizes M))
    (setf interval-size-list (get-cur-interval-sizes cur-interval-sizes M size-index MIN-INTERVAL-SIZE))
    (loop
        (incf size-index)
        (setf tmp-size-list (list))

        (loop for arr in interval-size-list do 
            (setf tmp-size-list (append tmp-size-list (get-cur-interval-sizes arr M size-index MIN-INTERVAL-SIZE)))
        )
        (setf interval-size-list tmp-size-list)
        (when (= size-index (- M 2)) (return))
    )
    interval-size-list
)

; повертає розмірності інтервалів після перебору можливих
(defun divide-into-interval-sizes (N M min-interval-size) 
    (setf interval-sizes (make-array M))
    (loop for i from 0 to (- M 1) do 
        (if (= i 0) 
            (setf (aref interval-sizes i) (- N (* min-interval-size (- M 1)))) 
            (setf (aref interval-sizes i) min-interval-size)
        )
    )
    (get-possible-interval-sizes-list interval-sizes M MIN-INTERVAL-SIZE)
)

; здійснює перетворення розмірностей інтервалів на індекси інтервалів
(defun transform-interval-size-to-index-list (interval-sizes M)
    (setf index-list (list))
    (loop for cur-interval-sizes in interval-sizes do
        (setf tmp-list (list Nil))
        (setf sum-interval-size 0)
        (setf tmp-index-intervals (make-array M))
        (loop for i from 0 to (- M 1) do
            (setf (aref tmp-index-intervals i) sum-interval-size)
            (setf sum-interval-size (+ sum-interval-size (aref cur-interval-sizes i)))
        )
        (setf (car tmp-list) tmp-index-intervals)
        (setf index-list (append index-list tmp-list))
    )
    index-list
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
; переведення числового ряду на лінгвістичний рядок
(defun translate-to-ling-chain (arr alphabet intervals N M)
    (setf ling-chain (make-array N))
    (loop for i from 0 to (- N 1) do
        (loop for j from 0 to (- M  1) do
            (when (and 
                    (<= (aref (aref intervals j) 0) (aref arr i)) 
                    (<= (aref arr i) (aref (aref intervals j) 1))
                )
                (setf (aref ling-chain i) (aref alphabet j))
                (return)
            )
        )
    )
    ling-chain
)

;----------------------------------
; виведення масиву
(defun print-array (arr size) (loop for i from 0 to (- size 1) do (write (aref arr i)) (princ " ")))

; виведення лінгвістичного рядка, що використовую виведення масиву
(defun print-ling-chain (arr size) (print "Linguistic chain") (print-array arr size))

;----------------------------------
; головна процедура
(defun main ()
    ; необхідні дані для виконання програми
    (setf M 4)
    (setf N 14)
    (setf MIN-INTERVAL-SIZE 2)
    (setf alphabet (make-array M :initial-contents '("A" "B" "C" "D")))
    (setf numbers (make-array N :initial-contents '(13 1 2 15 -10 20 3 12 4 5 6 26 7 8)))  
    
    (setf intervals (get-probable-intervals (sort (get-array-copy numbers N) #'< ) N M MIN-INTERVAL-SIZE))
    (print-ling-chain (
        translate-to-ling-chain numbers alphabet intervals N M
    ) N)
)

; вимірювання часу виконання програми
(time (main))


