; допоміжна функція, що округляє float до n знаків після коми - необхідно проти помилки floating point undeflow
(defun cut-float-digits (a digits)
    (setf float-d (expt 10 digits))
    (setf a (* a float-d))
    (setf a (fround a))
    (setf a (/ a float-d))
    a
 )

(defun fisher-z-p (a b) (- (fisher-z-fn b) (fisher-z-fn a)) )

(defun fisher-z-fn (x) 
    (defconstant MU 6.0)
    (defconstant SIGMA 7.0)
    (defconstant D1 3.0)
    (defconstant D2 2.0)

    (/ (incomplete-beta (/ D1 2) (/ D2 2) (count-xi (count-z x MU SIGMA) D1 D2)) (complete-beta (/ D1 2) (/ D2 2)))
)

(defun count-z (x mu sigma) (/ (- x mu) sigma))

(defun count-xi (z d1 d2) ( /
    (* d2 (exp (* 2 z))) 
    (+ d1 (* d2 (exp (* 2 z))))
    ))

(defun complete-beta (a b) 
    (exp (- (+ (lgamma a) (lgamma b) ) (lgamma (+ a b))))
)

(defun incomplete-beta (a b x)
    (when (or (< x 0.0) (> x 1.0)) Nil)
    (setf beta 0.0)
    (when (or (= x 0) (= x 1)) (setf beta 0.0))
    (setf beta 
        (exp (+ (+ (- (- (lgamma (+ a b)) (lgamma a)) (lgamma b)) (* a (log x))) (* b (log (- 1.0 b))))))
    (if (< x (/ (+ a 1) (+ a (+ b 2)))) 
        (/ (* beta (inc-beta-fraction a b x)) a)
        (- 1 (/ (* beta (inc-beta-fraction a b (- 1 x)) b)))
    )
)

(defun inc-beta-fraction (a b x)
    (defconstant MAXIT 200)
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

(print "beta")
(print (inc-beta-fraction 1.5 1.0 0.5))
