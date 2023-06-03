; utils.CLP

(defmodule UTILS
    (export deffunction sort-multifield)
    (export deffunction copy-multifield)
    (export deffunction print-multifield)
)

; допоміжна функція для сортування значень з багатьох полів
(deffunction UTILS::sort-multifield (?multifield)
    (bind ?size (length$ ?multifield))
    (loop-for-count (?i 1 (- ?size 1))
        (loop-for-count (?j (+ ?i 1) ?size)
            (bind ?a (+ (nth ?i ?multifield) 0))
            (bind ?b (+ (nth ?j ?multifield) 0))
            (if (> ?a ?b) then
                (bind ?tmp ?a)
                (bind ?multifield (replace$ ?multifield ?i ?i ?b))
                (bind ?multifield (replace$ ?multifield ?j ?j ?tmp))
            )
        )
    )
    ?multifield
)

; допоміжна функція для копіювання значень з багатьох полів
(deffunction UTILS::copy-multifield (?source)
    (bind ?dest (create$ (nth 1 ?source)))
    (bind ?size (length$ ?source))
    (loop-for-count (?i 2 ?size) (bind ?dest (insert$ ?dest ?i (nth ?i ?source))))
    ?dest
)

; функція для відображення значень з багатьох полів 
(deffunction UTILS::print-multifield (?multifield) 
    (bind ?size (length$ ?multifield))
    (loop-for-count (?i 1 ?size)
        (printout t "" (nth ?i ?multifield) " ")
    )
)
