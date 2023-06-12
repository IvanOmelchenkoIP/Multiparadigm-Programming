; ling.CLP

(defmodule LING 
    (export deffunction translate-to-ling-chain)
)

; формування лінгвістичного ряду
(deffunction translate-to-ling-chain (?N ?M ?alphabet ?probable-intervals ?numbers)
    (bind ?ling-chain (create$ TMP))
    (loop-for-count (?i 1 ?N)
        (bind ?val (nth ?i ?numbers))
        (loop-for-count (?j 1 ?M)
            (bind ?first-border (nth (- (* ?j 2) 1) ?probable-intervals))
            (bind ?second-border (nth (* ?j 2) ?probable-intervals))
            (if (and (<= ?first-border ?val) (<= ?val ?second-border)) then
                (bind ?ling-chain (insert$ ?ling-chain (+ ?i 1) (nth ?j ?alphabet)))
                (break)
            )
        )
    )
    (delete$ ?ling-chain 1 1)
)
