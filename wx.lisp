;; load this for wxMaxima which doesn't understand how to
;; display the other code properly.
(setf (get '$bra 'dissym) '((#\<) #\|))
(setf (get '$bra 'dimension) 'dimension-match)
(setf (get '$ket 'dissym) '((#\|) #\>))
(setf (get '$ket 'dimension) 'dimension-match)
