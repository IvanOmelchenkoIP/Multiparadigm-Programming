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
    (defconstant MAXIT 30)
    (defconstant EPS 3.0e-7)
    (defconstant FPMIN 1.0e-30)

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

        (print h)
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
        (exp (+ (+ (- (- (lgamma (+ a b)) (lgamma a)) (lgamma b)) (* a (log x))) (* b (log (- 1.0 b)))))
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

;----------------------------------
; копіювання масиву
(defun get-array-copy (arr size) 
    (setf copy (make-array N))
    (loop for i from 0 to (- size 1) do (setf (aref copy i) (aref arr i)))
    copy
)

;----------------------------------
; обчислення сумарної вірогідності кожного інтервалу з комбінації
(defun get-int-sum-probability (arr intervals N M)
    (setf probability 0.0)
    (loop for i from 0 to (- M 1) do 
        (print (aref intervals (+ i 1)))
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
            (setf cur-probability (get-int-sum-probability arr tmp-ind-intervals N M))
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
            (setf cur-probability (get-int-sum-probability arr tmp-ind-intervals N M))
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
    (print N)
    (print M)
    (print min-int-size)
    (loop for i from 0 to (- M 1) do 
        (if (= i 0) 
            (setf (aref int-sizes i) (- N (* min-int-size (- M 1)))) 
            (setf (aref int-sizes i) min-int-size)
        )
    )
    (print int-sizes)
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
    (defconstant MIN-INTERVAL-SIZE 2)
    (setf alphabet (make-array M :initial-contents '("A" "B" "C" "D")))
    (setf numbers (make-array N :initial-contents '(13 1 2 15 -10 20 3 12 4 5 6 26 7 8)))  
    
    (setf sorted (sort (get-array-copy numbers N) #'< ))
    (setf interval-inds (divide-into-intervals sorted N M MIN-INTERVAL-SIZE))
    ;(print-ling-chain ( 
    ;    (translate-to-ling-chain arr alphabet 
    ;        (get-intervals numbers 
    ;            (divide-into-intervals  N M)
    ;        N M) 
    ;    N M)
    ;    )
    ;N)
)

; вимірювання часу виконання програми
(time (main))

