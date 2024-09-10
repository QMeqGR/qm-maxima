
(defun $complex_number_p (e)
  (complex-number-p e #'$numberp))

;; This code is courtesy of Robert Dodier

;; arrange for dagger(a) to be displayed as a†
(setf (get '$dagger 'dissym) '(#\†))
(setf (get '$dagger 'dimension) 'dimension-postfix)

;; arrange for bra(a) to be displayed as <a|
(defun dimension-bra (e so-far) (dimension-match `((bra) ,@(rest
(first (rest e)))) so-far))
(setf (get '$bra 'dimension) 'dimension-bra)
(setf (get 'bra 'dissym) '((#\<) #\|))

;; arrange for ket(a) to be displayed as |a>
(defun dimension-ket (e so-far) (dimension-match `((ket) ,@(rest
(first (rest e)))) so-far))
(setf (get '$ket 'dimension) 'dimension-ket)
(setf (get 'ket 'dissym) '((#\|) #\>))

;; detect wxmaxima and use simpler bra/ket representations
;; (the above ones don't work with xml for some reason)
 ;; (if (> (first (last $wxplot_size)) 0)
 ;;     (progn
 ;;  (setf (get '$bra 'dissym) '((#\<) #\|))
 ;;  (setf (get '$bra 'dimension) 'dimension-match)
 ;;  (setf (get '$ket 'dissym) '((#\|) #\>))
 ;;  (setf (get '$ket 'dimension) 'dimension-match)
 ;; ))
