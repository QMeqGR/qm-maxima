(displa-def motimes dim-motimes)
(defprop munaryotimes (#\⊗ #\space) dissym)
(defprop cdot (#\· #\space) dissym)

;; This is a bit of a hack. I used the mplus code from displa.lisp
;; and made some small modifications. ehm
(defun dim-motimes (form result) 
  ;; If only 0 or 1 arguments, then print "⊗"() or ⊗A
  (cond ((and (null (cddr form))
              (not (member (cadar form) '(trunc exact) :test #'eq)))
         (if (null (cdr form))
             (dimension-function form result)
             (dimension-prefix (cons '(cdot) (cdr form)) result)))
        (t (setq result (dimension (cadr form) result lop 'mplus 0 0))
           (checkbreak result width)
           (do ((l (cddr form) (cdr l))
                (w width) (h height) (d depth) (cnt 0)
                (trunc (member 'trunc (cdar form) :test #'eq)) (dissym))
               ((null l) (cond (trunc
                                (setq width (+ 8 w) height h depth d)
                                (push-string " ⊗ . . ." result)))
                result)
             (if (mmminusp (car l))
                 (setq dissym '(#\space #\- #\space) form (cadar l))
		 (if (eq cnt 0)
		     (progn 
		       (setq dissym '(#\· ) form (car l))
		       (setf cnt (+ cnt 1) ) )
		     (setq dissym '(#\space #\⊗ #\space) form (car l)) )
		 )
             (cond ((and (not trunc) (null (cdr l)))
                    (setq result (dimension form (append dissym result)
                                            'mplus rop (+ 3 w) right)
                          width (+ 3 w width)
                          height (max h height)
                          depth (max d depth))
                    (return result))
                   (t (setq result
                            (dimension form (append dissym result)
                                       'mplus 'mplus (+ 3 w) 0)
                            w (+ 3 w width)
                            h (max h height)
                            d (max d depth))
                      (checkbreak result w)))))))

(setf (get '$tpket 'dimension) 'dim-motimes)
(setf (get '$tpbra 'dimension) 'dim-motimes)

