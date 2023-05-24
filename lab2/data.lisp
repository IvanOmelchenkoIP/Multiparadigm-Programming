
;(defun loop-break (max) (
 ;   (loop for i from 1 to 3 do ((print i))) ;when (< m max)) 
;))

;(loop-break 4)
(loop for x from 1 to 20
   if(evenp x)
   do (print x)
)

(defun power (x y) (if (> y 0) (power-count x y) (/ 1.0 (power-count x (* -1 y)))))
(defun power-count (x y) (if (= y 0) 1 (* x (power x (- y 1)))))
(print (power 10 -1))