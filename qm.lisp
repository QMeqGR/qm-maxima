(defun $complex_number_p (e)
  (complex-number-p e #'$numberp))

;; arrange for bra(a) to be displayed as <a|

(setf (get '$bra 'dissym) '((#\<) #\|))
(setf (get '$bra 'dimension) 'dimension-match)

;; arrange for ket(a) to be displayed as |a>

(setf (get '$ket 'dissym) '((#\|) #\>))
(setf (get '$ket 'dimension) 'dimension-match)

;; arrange for dagger(a) to be displayed as a†

(setf (get '$dagger 'dissym) '(#\†))
(setf (get '$dagger 'dimension) 'dimension-postfix)
