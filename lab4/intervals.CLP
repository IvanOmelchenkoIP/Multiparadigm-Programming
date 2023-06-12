; intervals.CLP

(defmodule INTERVALS
    (export deffunction divide-into-intervals)
    (import UTILS deffunction copy-multifield)
    (import DISTRIBUTION deffunction z-fisher-dist-p)
)

; визначення початкових розмірів інтервалів
(deffunction init-interval-sizes (?N ?M ?min-interval-size)
    (bind ?interval-sizes (create$ (- ?N (* ?min-interval-size (- ?M 1)))))
    (loop-for-count (?i 2 ?M) (bind ?interval-sizes (insert$ ?interval-sizes ?i ?min-interval-size)))
    ?interval-sizes
)

; допоміжна функція для визначення всіх можливих індексів інтервалів
(deffunction possible-interval-indexes (?given-interval-sizes ?M ?size-index ?min-interval-size)
    (bind ?cur-interval-sizes (copy-multifield ?given-interval-sizes))
    (bind ?next-size-index (+ ?size-index 1))
    (bind ?ind 2)
    (bind ?interval-indexes (create$ TMP))
    (while (>= (nth ?size-index ?cur-interval-sizes) ?min-interval-size)
        (if (< ?next-size-index (- ?M 1)) then
            (bind ?partial-interval-indexes 
                (possible-interval-indexes ?cur-interval-sizes ?M ?next-size-index ?min-interval-size)
            )
            (bind ?interval-indexes 
                (insert$ ?interval-indexes ?ind ?partial-interval-indexes)
            )
            (bind ?tmp-size (length$ ?partial-interval-indexes))
            (bind ?ind (+ ?ind ?tmp-size))
        )
        (if (= ?next-size-index (- ?M 1)) then
            (bind ?sum-size 1)
            (loop-for-count (?i 1 ?M)
                (bind ?interval-indexes (insert$ ?interval-indexes ?ind ?sum-size))
                (bind ?sum-size (+ ?sum-size (nth ?i ?cur-interval-sizes)))
                (bind ?ind (+ ?ind 1))
            )
        )
        (bind ?cur-interval-sizes
            (replace$ ?cur-interval-sizes ?size-index ?size-index (- (nth ?size-index ?cur-interval-sizes) 1))
        )
        (bind ?cur-interval-sizes
            (replace$ ?cur-interval-sizes ?next-size-index ?next-size-index (+ (nth ?next-size-index ?cur-interval-sizes) 1))
        )
    )
    (bind ?interval-indexes (delete$ ?interval-indexes 1 1))
    ?interval-indexes
)

; визначення найбільш подходящих інтервалів
(deffunction get-most-probable-intervals (?N ?M ?intervals ?numbers) 
    (bind ?intervals-ind 1)
    (bind ?size (length$ ?intervals))
    (bind ?probable-intervals (create$ TMP))
    (bind ?max-probability 0.0)
    (while (< ?intervals-ind ?size)
        (bind ?cur-probability 0.0)
        (bind ?cur-ind 2)
        (bind ?cur-intervals (create$ TMP))
        (loop-for-count (?i 1 ?M)
            (bind ?a tmp)
            (bind ?b tmp)
            (bind ?a (nth (nth ?intervals-ind ?intervals) ?numbers))
            (if (= ?i ?M) 
                then (bind ?b (nth ?N ?numbers)) 
                else (bind ?b (nth (- (nth (+ ?intervals-ind 1) ?intervals) 1) ?numbers))
            )
            (bind ?intervals-ind (+ ?intervals-ind 1))

            (bind ?cur-probability (+ ?cur-probability (z-fisher-dist-p ?a ?b)))

            (bind ?cur-intervals (insert$ ?cur-intervals ?cur-ind ?a))
            (bind ?cur-ind (+ ?cur-ind 1))
            (bind ?cur-intervals (insert$ ?cur-intervals ?cur-ind ?b))
            (bind ?cur-ind (+ ?cur-ind 1))
        )
        (if (> ?cur-probability ?max-probability) then
            (bind ?max-probability ?cur-probability)
            (bind ?probable-intervals (delete$ ?cur-intervals 1 1))
        )
    )
    ?probable-intervals
)

; знаходження необхідного інтервалу
(deffunction INTERVALS::divide-into-intervals (?N ?M ?min-interval-size ?numbers)
    (bind ?interval-sizes (init-interval-sizes ?N ?M ?min-interval-size))
    (bind ?possible-intervals (possible-interval-indexes ?interval-sizes ?M 1 ?min-interval-size))
    (get-most-probable-intervals ?N ?M ?possible-intervals ?numbers)
)