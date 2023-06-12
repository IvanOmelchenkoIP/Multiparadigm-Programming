; main.CLP

(defmodule MAIN
    (import DISTRIBUTION deffunction z-fisher-dist-p)
    (import UTILS deffunction sort-multifield)
    (import UTILS deffunction copy-multifield)
    (import UTILS deffunction print-multifield)
    (import LING deffunction translate-to-ling-chain)
    (import INTERVALS deffunction divide-into-intervals)
)


; головна функція
(deffunction main (?numbers ?alphabet)
    (bind ?N (length$ ?numbers))
    (bind ?M (length$ ?alphabet))
    (bind ?MIN-INTERVAL-SIZE 2)

    (bind ?sorted-numbers (sort-multifield (copy-multifield ?numbers)))
    (bind ?intervals (divide-into-intervals ?N ?M ?MIN-INTERVAL-SIZE ?sorted-numbers))
    (translate-to-ling-chain ?N ?M ?alphabet ?intervals ?numbers)
)

(deftemplate init-data
    (multislot numbers)
    (multislot alphabet)
)
(deffacts init-main-data
    (init-data
        (numbers (create$ 13 1 2 15 -10 20 3 12 4 5 6 26 7 8))
        (alphabet (create$ A B C D))
    )
)

(defrule main-rule
    ?fact <- (init-data)
    => 
    (bind ?start (time))
    (bind ?ling-chain 
        (main 
            (fact-slot-value ?fact numbers)
            (fact-slot-value ?fact alphabet)
        )
    )
    (bind ?end (time))
    (printout t "Linguistic chain" crlf)
    (print-multifield ?ling-chain)
    (printout t "" crlf)
    (printout t "Time: " (- ?end ?start) "s" crlf)
)