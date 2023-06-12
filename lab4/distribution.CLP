; distribution.CLP

(defmodule DISTRIBUTION
    (export deffunction z-fisher-dist-p)
)

; допоміжні функції для обчислення gamma і lgamma
(deffunction gamma (?x) (* (sqrt (* (/ (pi) ?x) 2)) (** (/ ?x (exp 1)) ?x)))

(deffunction lgamma (?x) (if (= ?x 1) then 0 else (log (gamma ?x))))

; допоміжна функція для обчислення непоіної бета функції
(deffunction inc-beta-fraction (?a ?b ?x ?MAXIT ?EPS ?FPMIN)
    (bind ?qab (+ ?a ?b))
    (bind ?qap (+ ?a 1.0))
    (bind ?qam (- ?a 1.0))
    (bind ?c 1.0)
    (bind ?d (- 1.0 (/ (* ?qab ?x) ?qap)))
    (if (< (abs ?d) ?FPMIN) then (bind ?d ?FPMIN))
    (bind ?d (/ 1.0 ?d))
    (bind ?h ?d)
    (loop-for-count (?i 1 ?MAXIT)
        (bind ?m ?i)
        (bind ?m2 (* 2.0 ?m))
        (bind ?aa (/ (* ?m (* (- ?b ?m) ?x)) (* (+ ?qam ?m2) (+ ?a ?m2))))
        (bind ?d (+ 1.0 (* ?aa ?d)))
        (if (< (abs ?d) ?FPMIN) then (bind ?d ?FPMIN))
        (bind ?c (+ 1.0 (/ ?aa ?c)))
        (if (< (abs ?c) ?FPMIN) then (bind ?c ?FPMIN))
        (bind ?d (/ 1.0 ?d))
        (bind ?h (* ?h (* ?d ?c)))
        (bind ?aa (/ 
                    (* (* (* (- 0 1) (+ ?a ?m)) (+ ?qab ?m)) ?x) 
                    (* (+ ?a ?m2) (+ ?qap ?m2))
                )
        )
        (bind ?d (+ 1.0 (* ?aa ?d)))
        (if (< (abs ?d) ?FPMIN) then (bind ?d ?FPMIN))
        (bind ?c (+ 1.0 (/ ?aa ?c)))
        (if (< (abs ?c) ?FPMIN) then (bind ?c ?FPMIN))
        (bind ?d (/ 1.0 ?d))
        (bind ?del (* ?d ?c))
        (bind ?h (* ?h ?del))
        (if (< (abs (- ?del 1.0)) ?EPS) then (break))
    )
    (if (> ?m ?MAXIT) then nil else ?h)
)

; обчислення неповної бета-функції
(deffunction incomplete-beta (?a ?b ?x ?MAXIT ?EPS ?FPMIN)
    (if (or (< ?x 0.0) (> ?x 1.0)) then nil)
    (bind ?beta 0.0)
    (if (or (= ?x 0) (= ?x 1)) then (bind ?beta 0.0))
    (bind ?beta 
        (exp (+ (+ (- (- (lgamma (+ ?a ?b)) (lgamma ?a)) (lgamma ?b)) (* ?a (log ?x))) (* ?b (log (- 1 ?x)))) )
    )
    (if (< ?x (/ (+ ?a 1) (+ ?a (+ ?b 2)))) 
        then (/ (* ?beta (inc-beta-fraction ?a ?b ?x ?MAXIT ?EPS ?FPMIN)) ?a)
        else (- 1 (/ (* ?beta (inc-beta-fraction ?a ?b (- 1 ?x) ?MAXIT ?EPS ?FPMIN)) ?b))
    )
)

; обчислення повної бета-функції
(deffunction complete-beta (?a ?b) (/ (* (gamma ?a) (gamma ?b)) (gamma (+ ?a ?b))))

; допоможні функції для обчислення 
(deffunction count-z (?x ?mu ?sigma) (/ (- ?x ?mu) ?sigma))
(deffunction count-xi (?z ?d1 ?d2) 
    ( /
        (* ?d2 (exp (* 2 ?z))) 
        (+ ?d1 (* ?d2 (exp (* 2 ?z))))
    )
)

; обчислення функції z-розподілу Фішера
(deffunction z-fisher-dist-fn (?x) 
    (bind ?d1 3.0)
    (bind ?d2 2.0)
    (bind ?mu 6.0)
    (bind ?sigma 7.0)
    (bind ?MAXIT 200)
    (bind ?EPS 3.0e-7)
    (bind ?FPMIN 1.0e-30)
    (bind ?z (count-z ?x ?mu ?sigma))
    (bind ?xi (count-xi ?z ?d1 ?d2))
    (bind ?half-d1 (/ ?d1 2))
    (bind ?half-d2 (/ ?d2 2))
    (bind ?inc-beta (incomplete-beta ?half-d1 ?half-d2 ?xi ?MAXIT ?EPS ?FPMIN))
    (bind ?beta (complete-beta ?half-d1 ?half-d2))
    (/ ?inc-beta ?beta)
)

; обчислення ймовірності потрапляння в інтервал
(deffunction DISTRIBUTION::z-fisher-dist-p (?a ?b) (- (z-fisher-dist-fn ?b) (z-fisher-dist-fn ?a)))