
(defun cut-float-digits (a digits)
    (setf float-d (expt 10 digits))
    (setf a (* a (* float-d)))
    (setf a1 a)
    (setf a (fround a))
    (setf a (/ a float-d))
    a
 )

(print (cut-float-digits 0.235 35))
