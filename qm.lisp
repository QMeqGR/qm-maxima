
(defun $complex_number_p (e)
  (complex-number-p e #'$numberp))

;; This code is courtesy of Robert Dodier
;; arrange for bra(a) to be displayed as <a|

(defun dimension-bra (e so-far) (dimension-match `((bra) ,@(rest
(first (rest e)))) so-far))
(setf (get '$bra 'dimension) 'dimension-bra)
(setf (get 'bra 'dissym) '((#\<) #\|))

;; (setf (get '$bra 'dissym) '((#\<) #\|))
;; (setf (get '$bra 'dimension) 'dimension-match)

;; arrange for ket(a) to be displayed as |a>

(defun dimension-ket (e so-far) (dimension-match `((ket) ,@(rest
(first (rest e)))) so-far))
(setf (get '$ket 'dimension) 'dimension-ket)
(setf (get 'ket 'dissym) '((#\|) #\>))

;; (setf (get '$ket 'dissym) '((#\|) #\>))
;; (setf (get '$ket 'dimension) 'dimension-match)

;; arrange for dagger(a) to be displayed as a†

(setf (get '$dagger 'dissym) '(#\†))
(setf (get '$dagger 'dimension) 'dimension-postfix)
