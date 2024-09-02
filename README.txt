Package qm
**********

1 Package qm
************

1.1 Introduction to package qm
==============================

The ‘qm’ package was written by Eric Majzoub, University of Missouri.
Email: majzoube-at-umsystem.edu

   The package is loaded with: ‘load(qm);’

   The ‘qm’ package provides functions and standard definitions to solve
quantum mechanics problems in a finite dimensional Hilbert space.  For
example, one can calculate the outcome of Stern-Gerlach experiments
using built-in definitions of the Sx, Sy, and Sz operators for arbitrary
spin, e.g.  ‘s={1/2, 1, 3/2, ...}’.  For spin-1/2 the standard basis
states in the ‘x’, ‘y’, and ‘z’-basis are available as ‘{xp,xm}’,
‘{yp,ym}’, and ‘{zp,zm}’, respectively.  One can create general ket
vectors with arbitrary but finite dimension and perform standard
computations such as expectation value, variance, etc.  The angular
momentum <|j,m>> representation of kets is also available.  Tensor
product states for multiparticle systems can be created to perform
calculations on those systems.

   Let us consider a simple example involving spin-1/2 particles.  A bra
vector in the ‘z’-basis may be written as

           ‘<psi| = a <z+| + b <z-|’.

   The bra will be represented in Maxima by the row vector ‘[a b]’,
where the basis vectors are

           ‘<z+| = [1 0]’

and

           ‘<z-| = [0 1]’.

1.1.1 Types of kets and bras
----------------------------

There are two types of kets and bras available in this package, the
first type is given by a _matrix representation_, as in the above
example.  ‘mket’s are column vectors and ‘mbra’s are row vectors, and
their components are entered as Maxima _lists_ in the ‘mbra’ and ‘mket’
functions.  The second type of bra or ket is _abstract_; there is no
matrix representation.  Abstract bras and kets are entered using the
‘bra’ and ‘ket’ functions using Maxima lists for the elements.  These
general kets are displayed in Dirac notation.  Abstract bras and kets
are used for both the ‘(j,m)’ representation of states and also for
tensor products.  For example, a tensor product of two ket vectors ‘|a>’
and ‘|b>’ is input as ‘ket([a,b])’ and displayed as

           ‘|[a,b]>’      (general ket)

Note that abstract kets and bras are _assumed to be orthonormal_.  These
general bras and kets may be used to build arbitrarily large tensor
product states.  Tensor product states in the matrix representation are
also available through the ‘tpket’ and ‘tpbra’ commands.

   The following examples illustrate some of the basic capabilities of
the ‘qm’ package.  Here both abstract, and concrete (matrix
representation) kets are shown.

     (%i1) ket([a,b])+ket([c,d]);
     (%o1)                         |[c, d]> + |[a, b]>
     (%i2) mket([a,b])+mket([c,d]);
                                        [ c + a ]
     (%o2)                              [       ]
                                        [ d + b ]
     (%i3) bell:(1/sqrt(2))*(ket([u,d])-ket([d,u]));
                                   |[u, d]> - |[d, u]>
     (%o3)                         -------------------
                                         sqrt(2)

   Note that ‘ket([a,b])’ is treated as tensor product of states ‘a’ and
‘b’ as shown below.

     (%i1) braket(bra([a1,b1]),ket([a2,b2]));
     (%o1)                kron_delta(a1, a2) kron_delta(b1, b2)

   Next, tensor products of the spin-1/2 basis states ‘{zp,zm}’ are
shown in the matrix representation.

     (%i1) tpket([zp,zm]);
                                          [ 1 ]  [ 0 ]
     (%o1)                       [tpket, [[   ], [   ]]]
                                          [ 0 ]  [ 1 ]

   Constants that multiply kets and bras must be declared complex by the
user in order for the dagger function to properly conjugate such
constants.  The example below illustrates this behavior.

     (%i1) declare([a,b],complex);
     (%o1)                                done
     (%i2) psi:a*ket([1])+b*ket([2]);
     (%o2)                          |[2]> b + |[1]> a
     (%i3) psidag:dagger(psi);
     (%o3)               <[2]| conjugate(b) + <[1]| conjugate(a)
     (%i4) psidag . psi;
     (%o4)                   b conjugate(b) + a conjugate(a)

   The following shows how to declare a ket with both real and complex
components in the matrix representation.

     (%i1) declare([c1,c2],complex,r,real);
     (%o1)                                done
     (%i2) k:mket([c1,c2,r]);
                                         [ c1 ]
                                         [    ]
     (%o2)                               [ c2 ]
                                         [    ]
                                         [ r  ]
     (%i3) b:dagger(k);
     (%o3)                 [ conjugate(c1)  conjugate(c2)  r ]
     (%i4) b . k;
                         2
     (%o4)              r  + c2 conjugate(c2) + c1 conjugate(c1)

1.2 Functions and Variables for qm
==================================

 -- Variable: hbar
     Planck's constant divided by ‘2*%pi’.  ‘hbar’ is not given a
     floating point value, but is declared to be a real number greater
     than zero.

 -- Function: ket ([k_{1},k_{2},...])
     ‘ket’ creates a general state ket, or tensor product, with symbols
     ‘k_{i}’ representing the states.  The state kets ‘k_{i}’ are
     assumed to be orthonormal.

     (%i1) k:ket([u,d]);
     (%o1)                              |[u, d]>
     (%i2) b:bra([u,d]);
     (%o2)                              <[u, d]|
     (%i3) b . k;
     (%o3)                                  1

 -- Function: ketp (abstract ket)
     ‘ketp’ is a predicate function for abstract kets.  It returns
     ‘true’ for abstract ‘ket’s and ‘false’ for anything else.

 -- Function: bra ([b_{1},b_{2},...])
     ‘bra’ creates a general state bra, or tensor product, with symbols
     ‘b_{i}’ representing the states.  The state bras ‘b_{i}’ are
     assumed to be orthonormal.

     (%i1) k:ket([u,d]);
     (%o1)                              |[u, d]>
     (%i2) b:bra([u,d]);
     (%o2)                              <[u, d]|
     (%i3) b . k;
     (%o3)                                  1

 -- Function: brap (abstract bra)
     ‘brap’ is a predicate function for abstract bras.  It returns
     ‘true’ for abstract ‘bra’s and ‘false’ for anything else.

 -- Function: mket ([c_{1},c_{2},...])
     ‘mket’ creates a _column_ vector of arbitrary finite dimension.
     The entries ‘c_{i}’ can be any Maxima expression.  The user must
     ‘declare’ any relevant constants to be complex.  For a matrix
     representation the elements must be entered as a list in ‘[...]’
     square brackets.

     (%i1) declare([c1,c2],complex);
     (%o1)                                done
     (%i2) mket([c1,c2]);
                                         [ c1 ]
     (%o2)                               [    ]
                                         [ c2 ]
     (%i3) facts();
     (%o3) [kind(hbar, real), hbar > 0, kind(c1, complex), kind(c2, complex)]

 -- Function: mketp (_vector_)
     ‘mketp’ is a predicate function that checks if its input is an
     mket, in which case it returns ‘true’, else it returns ‘false’.
     ‘mketp’ only returns ‘true’ for the matrix representation of a ket.

     (%i1) k:ket([a,b]);
     (%o1)                              |[a, b]>
     (%i2) mketp(k);
     (%o2)                                false
     (%i3) k:mket([a,b]);
                                          [ a ]
     (%o3)                                [   ]
                                          [ b ]
     (%i4) mketp(k);
     (%o4)                                true

 -- Function: mbra ([c_{1},c_{2},...])
     ‘mbra’ creates a _row_ vector of arbitrary finite dimension.  The
     entries ‘c_{i}’ can be any Maxima expression.  The user must
     ‘declare’ any relevant constants to be complex.  For a matrix
     representation the elements must be entered as a list in ‘[...]’
     square brackets.

     (%i1) kill(c1,c2);
     (%o1)                                done
     (%i2) mbra([c1,c2]);
     (%o2)                             [ c1  c2 ]
     (%i3) facts();
     (%o3)                    [kind(hbar, real), hbar > 0]

 -- Function: mbrap (_vector_)
     ‘mbrap’ is a predicate function that checks if its input is an
     mbra, in which case it returns ‘true’, else it returns ‘false’.
     ‘mbrap’ only returns ‘true’ for the matrix representation of a bra.

     (%i1) b:mbra([a,b]);
     (%o1)                              [ a  b ]
     (%i2) mbrap(b);
     (%o2)                                true

   Two additional functions are provided to create kets and bras in the
matrix representation.  These functions conveniently attempt to
automatically ‘declare’ constants as complex.  For example, if a list
entry is ‘a*sin(x)+b*cos(x)’ then only ‘a’ and ‘b’ will be ‘declare’-d
complex and not ‘x’.

 -- Function: autoket ([a_{1},a_{2},...])
     ‘autoket’ takes a list [‘a_{1},a_{2},...’] and returns a ket with
     the coefficients ‘a_{i}’ ‘declare’-d complex.  Simple expressions
     such as ‘a*sin(x)+b*cos(x)’ are allowed and will ‘declare’ only the
     coefficients as complex.

     (%i1) autoket([a,b]);
                                          [ a ]
     (%o1)                                [   ]
                                          [ b ]
     (%i2) facts();
     (%o2)  [kind(hbar, real), hbar > 0, kind(a, complex), kind(b, complex)]
     (%i1) autoket([a*sin(x),b*sin(x)]);
                                      [ a sin(x) ]
     (%o1)                            [          ]
                                      [ b sin(x) ]
     (%i2) facts();
     (%o2)  [kind(hbar, real), hbar > 0, kind(a, complex), kind(b, complex)]

 -- Function: autobra ([a_{1},a_{2},...])
     ‘autobra’ takes a list [‘a_{1},a_{2},...’] and returns a bra with
     the coefficients ‘a_{i}’ ‘declare’-d complex.  Simple expressions
     such as ‘a*sin(x)+b*cos(x)’ are allowed and will ‘declare’ only the
     coefficients as complex.

     (%i1) autobra([a,b]);
     (%o1)                              [ a  b ]
     (%i2) facts();
     (%o2)  [kind(hbar, real), hbar > 0, kind(a, complex), kind(b, complex)]
     (%i1) autobra([a*sin(x),b]);
     (%o1)                           [ a sin(x)  b ]
     (%i2) facts();
     (%o2)  [kind(hbar, real), hbar > 0, kind(a, complex), kind(b, complex)]

 -- Function: dagger (_vector_)
     ‘dagger’ is the quantum mechanical _dagger_ function and returns
     the ‘conjugate’ ‘transpose’ of its input.  Arbitrary constants must
     be ‘declare’-d complex for dagger to produce the conjugate.

     (%i1) dagger(mbra([%i,2]));
                                        [ - %i ]
     (%o1)                              [      ]
                                        [  2   ]

 -- Function: braket (psi,phi)
     Given two kets ‘psi’ and ‘phi’, ‘braket’ returns the quantum
     mechanical bracket ‘<psi|phi>’.  The vector ‘psi’ may be input as
     either a ‘ket’ or ‘bra’.  If it is a ‘ket’ it will be turned into a
     ‘bra’ with the ‘dagger’ function before the inner product is taken.
     The vector ‘phi’ must always be a ‘ket’.

     (%i1) declare([a,b,c],complex);
     (%o1)                                done
     (%i2) braket(mket([a,b,c]),mket([a,b,c]));
     (%o2)          c conjugate(c) + b conjugate(b) + a conjugate(a)
     (%i3) braket(ket([a1,b1,c1]),ket([a2,b2,c2]));
     (%o3)      kron_delta(a1, a2) kron_delta(b1, b2) kron_delta(c1, c2)

 -- Function: norm (psi)
     Given a ‘ket’ or ‘bra’ ‘psi’, ‘norm’ returns the square root of the
     quantum mechanical bracket ‘<psi|psi>’.  The vector ‘psi’ must
     always be a ‘ket’, otherwise the function will return ‘false’.

     (%i1) declare([a,b,c],complex);
     (%o1)                                done
     (%i2) norm(mket([a,b,c]));
     (%o2)       sqrt(c conjugate(c) + b conjugate(b) + a conjugate(a))

 -- Function: magsqr (c)
     ‘magsqr’ returns ‘conjugate(c)*c’, the magnitude squared of a
     complex number.

     (%i1) declare([a,b,c,d],complex);
     (%o1)                                done
     (%i2) A:braket(mket([a,b]),mket([c,d]));
     (%o2)                   conjugate(b) d + conjugate(a) c
     (%i3) P:magsqr(A);
     (%o3) (conjugate(b) d + conjugate(a) c) (b conjugate(d) + a conjugate(c))

1.2.1 Handling general kets and bras
------------------------------------

General kets and bras are, as discussed, created without using a list
when giving the arguments.  The following examples show how general kets
and bras can be manipulated.

     (%i1) ket([a])+ket([b]);
     (%o1)                            |[b]> + |[a]>
     (%i2) braket(bra([a]),ket([b]));
     (%o2)                          kron_delta(a, b)
     (%i3) braket(bra([a])+bra([c]),ket([b]));
     (%o3)                 kron_delta(b, c) + kron_delta(a, b)

1.2.2 Spin-1/2 state kets and associated operators
--------------------------------------------------

Spin-1/2 particles are characterized by a simple 2-dimensional Hilbert
space of states.  It is spanned by two vectors.  In the <z>-basis these
vectors are ‘{zp,zm}’, and the basis kets in the <z>-basis are ‘{xp,xm}’
and ‘{yp,ym}’ respectively.

 -- Function: zp
     Return the <|z+>> ket in the <z>-basis.

 -- Function: zm
     Return the <|z->> ket in the <z>-basis.

 -- Function: xp
     Return the <|x+>> ket in the <z>-basis.

 -- Function: xm
     Return the <|x->> ket in the <z>-basis.

 -- Function: yp
     Return the <|y+>> ket in the <z>-basis.

 -- Function: ym
     Return the <|y->> ket in the <z>-basis.

     (%i1) zp;
                                          [ 1 ]
     (%o1)                                [   ]
                                          [ 0 ]
     (%i2) zm;
                                          [ 0 ]
     (%o2)                                [   ]
                                          [ 1 ]
     (%i1) yp;
                                       [    1    ]
                                       [ ------- ]
                                       [ sqrt(2) ]
     (%o1)                             [         ]
                                       [   %i    ]
                                       [ ------- ]
                                       [ sqrt(2) ]
     (%i2) ym;
                                      [     1     ]
                                      [  -------  ]
                                      [  sqrt(2)  ]
     (%o2)                            [           ]
                                      [     %i    ]
                                      [ - ------- ]
                                      [   sqrt(2) ]
     (%i1) braket(xp,zp);
                                            1
     (%o1)                               -------
                                         sqrt(2)

   Switching bases is done in the following example where a <z>-basis
ket is constructed and the <x>-basis ket is computed.

     (%i1) declare([a,b],complex);
     (%o1)                                done
     (%i2) psi:mket([a,b]);
                                          [ a ]
     (%o2)                                [   ]
                                          [ b ]
     (%i3) psi_x:'xp*braket(xp,psi)+'xm*braket(xm,psi);
                         b         a              a         b
     (%o3)           (------- + -------) xp + (------- - -------) xm
                      sqrt(2)   sqrt(2)        sqrt(2)   sqrt(2)

1.2.3 Pauli matrices and Sz, Sx, Sy operators
---------------------------------------------

 -- Function: sigmax
     Returns the Pauli <x> matrix.

 -- Function: sigmay
     Returns the Pauli <y> matrix.

 -- Function: sigmaz
     Returns the Pauli <z> matrix.

 -- Function: Sx
     Returns the spin-1/2 <Sx> matrix.

 -- Function: Sy
     Returns the spin-1/2 <Sy> matrix.

 -- Function: Sz
     Returns the spin-1/2 <Sz> matrix.

     (%i1) sigmay;
                                      [ 0   - %i ]
     (%o1)                            [          ]
                                      [ %i   0   ]
     (%i2) Sy;
                                 [            %i hbar ]
                                 [    0     - ------- ]
                                 [               2    ]
     (%o2)                       [                    ]
                                 [ %i hbar            ]
                                 [ -------      0     ]
                                 [    2               ]

 -- Function: commutator (X,Y)
     Given two operators ‘X’ and ‘Y’, return the commutator ‘X . Y - Y .
     X’.

     (%i1) commutator(Sx,Sy);
                                [        2             ]
                                [ %i hbar              ]
                                [ --------      0      ]
                                [    2                 ]
     (%o1)                      [                      ]
                                [                    2 ]
                                [             %i hbar  ]
                                [    0      - -------- ]
                                [                2     ]

1.2.4 SX, SY, SZ operators for any spin
---------------------------------------

 -- Function: SX (s)
     ‘SX(s)’ for spin ‘s’ returns the matrix representation of the spin
     operator ‘Sx’.  Shortcuts for spin-1/2 are ‘Sx,Sy,Sz’, and for
     spin-1 are ‘Sx1,Sy1,Sz1’.

 -- Function: SY (s)
     ‘SY(s)’ for spin ‘s’ returns the matrix representation of the spin
     operator ‘Sy’.  Shortcuts for spin-1/2 are ‘Sx,Sy,Sz’, and for
     spin-1 are ‘Sx1,Sy1,Sz1’.

 -- Function: SZ (s)
     ‘SZ(s)’ for spin ‘s’ returns the matrix representation of the spin
     operator ‘Sz’.  Shortcuts for spin-1/2 are ‘Sx,Sy,Sz’, and for
     spin-1 are ‘Sx1,Sy1,Sz1’.

   Example:

     (%i1) SY(1/2);
                                 [            %i hbar ]
                                 [    0     - ------- ]
                                 [               2    ]
     (%o1)                       [                    ]
                                 [ %i hbar            ]
                                 [ -------      0     ]
                                 [    2               ]
     (%i2) SX(1);
                              [           hbar            ]
                              [    0     -------     0    ]
                              [          sqrt(2)          ]
                              [                           ]
                              [  hbar              hbar   ]
     (%o2)                    [ -------     0     ------- ]
                              [ sqrt(2)           sqrt(2) ]
                              [                           ]
                              [           hbar            ]
                              [    0     -------     0    ]
                              [          sqrt(2)          ]

1.2.5 Expectation value and variance
------------------------------------

 -- Function: expect (O,psi)
     Computes the quantum mechanical expectation value of the operator
     ‘O’ in state ‘psi’, ‘<psi|O|psi>’.

     (%i1) ev(expect(Sy,xp+ym),ratsimp);
     (%o1)                               - hbar

 -- Function: qm_variance (O,psi)
     Computes the quantum mechanical variance of the operator ‘O’ in
     state ‘psi’, ‘sqrt(<psi|O^{2}|psi> - <psi|O|psi>^{2})’.

     (%i1) ev(qm_variance(Sy,xp+ym),ratsimp);
                                         %i hbar
     (%o1)                               -------
                                            2

1.2.6 Angular momentum representation of kets and bras
------------------------------------------------------

1.2.6.1 Matrix representation of (j,m)-kets and bras
....................................................

The matrix representation of kets and bras in the ‘qm’ package are
represented in the ‘z’-basis.  To create a matrix representation of of a
ket or bra in the (j,m)-basis one uses the ‘spin_mket’ and ‘spin_mbra’
functions.

 -- Function: spin_mket (s,m_{s},[1,2])
     ‘spin_mket’ returns a ket in the ‘z’-basis for spin ‘s’ and
     z-projection ‘m_{s}’, for axis 1=X or 2=Y.

 -- Function: spin_mbra (s,m_{s},[1,2])
     ‘spin_mbra’ returns a bra in the ‘z’-basis for spin ‘s’ and
     z-projection ‘m_{s}’, for axis 1=X or 2=Y.

     (%i1) spin_mket(3/2,1/2,2);
                                     [  sqrt(3)   ]
                                     [  -------   ]
                                     [    3/2     ]
                                     [   2        ]
                                     [            ]
                                     [     %i     ]
                                     [    ----    ]
                                     [     3/2    ]
                                     [    2       ]
     (%o1)                           [            ]
                                     [     1      ]
                                     [    ----    ]
                                     [     3/2    ]
                                     [    2       ]
                                     [            ]
                                     [ sqrt(3) %i ]
                                     [ ---------- ]
                                     [     3/2    ]
                                     [    2       ]
     (%i2) spin_mbra(1,1,1);
                                    [ 1     1     1 ]
     (%o2)                          [ -  -------  - ]
                                    [ 2  sqrt(2)  2 ]

1.2.6.2 Abstract (j,m)-kets and bras
....................................

To create kets and bras in the <|j,m>> representation you use the
abstract ‘ket’ and ‘bra’ functions with ‘j,m’ as arguments, as in
‘ket([j,m])’ and ‘bra([j,m])’.

     (%i1) bra([3/2,1/2]);
                                          3  1
     (%o1)                              <[-, -]|
                                          2  2
     (%i2) ket([3/2,1/2]);
                                          3  1
     (%o2)                              |[-, -]>
                                          2  2

 -- Function: jmketp (jmket)
     ‘jmketp’ checks to see that the ket has an ‘m’-value that is in the
     set ‘{-j,-j+1,...,+j}’.

     (%i1) jmketp(jmket([j,m]));
     (%o1)                                false
     (%i2) jmketp(jmket([3/2,1/2]));
     (%o2)                                true

 -- Function: jmbrap (jmbra)
     ‘jmbrap’ checks to see that the bra has an ‘m’-value that is in the
     set ‘{-j,-j+1,...,+j}’.

 -- Function: jmcheck (j,m)
     ‘jmcheck’ checks to see that <m> is one of {-j, ..., +j}.

     (%i1) jmcheck(3/2,1/2);
     (%o1)                                true

 -- Function: JP (_jmket_)
     ‘JP’ is the ‘J_{+}’ operator.  It takes a ‘jmket’ ‘jmket(j,m)’ and
     returns ‘sqrt(j*(j+1)-m*(m+1))*hbar*jmket(j,m+1)’.

 -- Function: JM (_jmket_)
     ‘JM’ is the ‘J_{-}’ operator.  It takes a ‘jmket’ ‘jmket(j,m)’ and
     returns ‘sqrt(j*(j+1)-m*(m-1))*hbar*jmket(j,m-1)’.

 -- Function: Jsqr (_jmket_)
     ‘Jsqr’ is the ‘J^{2}’ operator.  It takes a ‘jmket’ ‘jmket(j,m)’
     and returns ‘(j*(j+1)*hbar^{2}*jmket(j,m)’.

 -- Function: Jz (_jmket_)
     ‘Jz’ is the ‘J_{z}’ operator.  It takes a ‘jmket’ ‘jmket(j,m)’ and
     returns ‘m*hbar*jmket(j,m)’.

   These functions are illustrated below.

     (%i1) k:ket([j,m]);
     (%o1)                              |[j, m]>
     (%i2) JP(k);
     (%o2)            hbar |[j, m + 1]> sqrt(j (j + 1) - m (m + 1))
     (%i3) JM(k);
     (%o3)            hbar |[j, m - 1]> sqrt(j (j + 1) - (m - 1) m)
     (%i4) Jsqr(k);
                                    2
     (%o4)                      hbar  j (j + 1) |[j, m]>
     (%i5) Jz(k);
     (%o5)                           hbar |[j, m]> m

1.2.7 Angular momentum and ladder operators
-------------------------------------------

 -- Function: SP (s)
     ‘SP’ is the raising ladder operator <S_{+}> for spin ‘s’.

 -- Function: SM (s)
     ‘SM’ is the raising ladder operator <S_{-}> for spin ‘s’.

   Examples of the ladder operators:

     (%i1) SP(1);
                            [ 0  sqrt(2) hbar       0       ]
                            [                               ]
     (%o1)                  [ 0       0        sqrt(2) hbar ]
                            [                               ]
                            [ 0       0             0       ]
     (%i2) SM(1);
                            [      0             0        0 ]
                            [                               ]
     (%o2)                  [ sqrt(2) hbar       0        0 ]
                            [                               ]
                            [      0        sqrt(2) hbar  0 ]

1.3 Rotation operators
======================

 -- Function: RX (s,t)
     ‘RX(s)’ for spin ‘s’ returns the matrix representation of the
     rotation operator ‘Rx’ for rotation through angle ‘t’.

 -- Function: RY (s,t)
     ‘RY(s)’ for spin ‘s’ returns the matrix representation of the
     rotation operator ‘Ry’ for rotation through angle ‘t’.

 -- Function: RZ (s,t)
     ‘RZ(s)’ for spin ‘s’ returns the matrix representation of the
     rotation operator ‘Rz’ for rotation through angle ‘t’.

     (%i1) RZ(1/2,t);
     Proviso: assuming 64*t # 0
                                  [     %i t         ]
                                  [   - ----         ]
                                  [      2           ]
                                  [ %e          0    ]
     (%o1)                        [                  ]
                                  [             %i t ]
                                  [             ---- ]
                                  [              2   ]
                                  [    0      %e     ]

1.4 Time-evolution operator
===========================

 -- Function: UU (H,t)
     ‘UU(H,t)’ is the time evolution operator for Hamiltonian ‘H’.  It
     is defined as the matrix exponential ‘matrixexp(-%i*H*t/hbar)’.

     (%i1) UU(w*Sy,t);
     Proviso: assuming 64*t*w # 0
                                [     t w         t w  ]
                                [ cos(---)  - sin(---) ]
                                [      2           2   ]
     (%o1)                      [                      ]
                                [     t w        t w   ]
                                [ sin(---)   cos(---)  ]
                                [      2          2    ]

1.5 Tensor products
===================

Tensor products are represented as lists in the ‘qm’ package.  The ket
tensor product ‘|z+,z+>’ could be represented as ‘ket([u,d])’, for
example, and the bra tensor product ‘<a,b|’ is represented as
‘bra([a,b])’ for states ‘a’ and ‘b’.  For a tensor product where the
identity is one of the elements of the product, substitute the string
‘Id’ in the ket or bra at the desired location.  See the examples below
for the use of the identity in tensor products.

 -- Function: tpket ([k_{1}, k_{2}, ...])
     ‘tpket’ produces a tensor product of kets ‘k_{i}’.  All of the
     elements must pass the ‘ketp’ predicate test to be accepted.  If a
     list is not used for the input kets, the tpket will be an abstract
     tensor product ket.

 -- Function: tpbra ([b_{1}, b_{2}, ...])
     ‘tpbra’ produces a tensor product of bras ‘b_{i}’.  All of the
     elements must pass the ‘brap’ predicate test to be accepted.  If a
     list is not used for the input bras, the tpbra will be an abstract
     tensor product bra.

 -- Function: tpketp (tpket)
     ‘tpketp’ checks to see that the ket has the 'tpket' marker.  Only
     the matrix representation will pass this test.

 -- Function: tpbrap (tpbra)
     ‘tpbrap’ checks to see that the bra has the 'tpbra' marker.  Only
     the matrix representation will pass this test.

 -- Function: tpbraket (B,K)
     ‘tpbraket’ takes the inner product of the tensor products ‘B’ and
     ‘K’.  The tensor products must be of the same length (number of
     kets must equal the number of bras).

   Examples below show how to create concrete (matrix representation)
tensor products and take the bracket of tensor products.

     (%i1) kill(a,b,c,d);
     (%o1)                                done
     (%i2) declare([a,b,c,d],complex);
     (%o2)                                done
     (%i3) tpbra([mbra([a,b]),mbra([c,d])]);
     (%o3)                    [tpbra, [[ a  b ], [ c  d ]]]
     (%i4) tpbra([dagger(zp),mbra([c,d])]);
     (%o4)                    [tpbra, [[ 1  0 ], [ c  d ]]]
     (%i1) K:tpket([zp,zm]);
                                          [ 1 ]  [ 0 ]
     (%o1)                       [tpket, [[   ], [   ]]]
                                          [ 0 ]  [ 1 ]
     (%i2) zpb:dagger(zp);
     (%o2)                              [ 1  0 ]
     (%i3) zmb:dagger(zm);
     (%o3)                              [ 0  1 ]
     (%i4) B:tpbra([zpb,zmb]);
     (%o4)                    [tpbra, [[ 1  0 ], [ 0  1 ]]]
     (%i5) tpbraket(K,B);
     (%o5)                                false
     (%i6) tpbraket(B,K);
     (%o6)                                  1

   Examples below show how to create abstract tensor products that
contain the identity element ‘Id’ and how to take the bracket of these
tensor products.

     (%i1) K:ket([a1,b1]);
     (%o1)                             |[a1, b1]>
     (%i2) B:bra([a2,b2]);
     (%o2)                             <[a2, b2]|
     (%i3) braket(B,K);
     (%o3)                kron_delta(a1, a2) kron_delta(b1, b2)
     (%i1) bra([a1,Id,c1]) . ket([a2,b2,c2]);
     (%o1)         |[-, b2, -]> kron_delta(a1, a2) kron_delta(c1, c2)
     (%i2) bra([a1,b1,c1]) . ket([Id,b2,c2]);
     (%o2)         <[a1, -, -]| kron_delta(b1, b2) kron_delta(c1, c2)

Appendix A Function and Variable index
**************************************

* Menu:

* autobra:                               Functions and Variables for qm.
                                                              (line 243)
* autoket:                               Functions and Variables for qm.
                                                              (line 224)
* bra:                                   Functions and Variables for qm.
                                                              (line 146)
* braket:                                Functions and Variables for qm.
                                                              (line 268)
* brap:                                  Functions and Variables for qm.
                                                              (line 158)
* commutator:                            Functions and Variables for qm.
                                                              (line 420)
* dagger:                                Functions and Variables for qm.
                                                              (line 258)
* expect:                                Functions and Variables for qm.
                                                              (line 479)
* JM:                                    Functions and Variables for qm.
                                                              (line 578)
* jmbrap:                                Functions and Variables for qm.
                                                              (line 564)
* jmcheck:                               Functions and Variables for qm.
                                                              (line 568)
* jmketp:                                Functions and Variables for qm.
                                                              (line 555)
* JP:                                    Functions and Variables for qm.
                                                              (line 574)
* Jsqr:                                  Functions and Variables for qm.
                                                              (line 582)
* Jz:                                    Functions and Variables for qm.
                                                              (line 586)
* ket:                                   Functions and Variables for qm.
                                                              (line 130)
* ketp:                                  Functions and Variables for qm.
                                                              (line 142)
* magsqr:                                Functions and Variables for qm.
                                                              (line 292)
* mbra:                                  Functions and Variables for qm.
                                                              (line 194)
* mbrap:                                 Functions and Variables for qm.
                                                              (line 208)
* mket:                                  Functions and Variables for qm.
                                                              (line 162)
* mketp:                                 Functions and Variables for qm.
                                                              (line 178)
* norm:                                  Functions and Variables for qm.
                                                              (line 282)
* qm_variance:                           Functions and Variables for qm.
                                                              (line 486)
* RX:                                    Functions and Variables for qm.
                                                              (line 631)
* RY:                                    Functions and Variables for qm.
                                                              (line 635)
* RZ:                                    Functions and Variables for qm.
                                                              (line 639)
* sigmax:                                Functions and Variables for qm.
                                                              (line 389)
* sigmay:                                Functions and Variables for qm.
                                                              (line 392)
* sigmaz:                                Functions and Variables for qm.
                                                              (line 395)
* SM:                                    Functions and Variables for qm.
                                                              (line 610)
* SP:                                    Functions and Variables for qm.
                                                              (line 607)
* spin_mbra:                             Functions and Variables for qm.
                                                              (line 510)
* spin_mket:                             Functions and Variables for qm.
                                                              (line 506)
* Sx:                                    Functions and Variables for qm.
                                                              (line 398)
* SX:                                    Functions and Variables for qm.
                                                              (line 438)
* Sy:                                    Functions and Variables for qm.
                                                              (line 401)
* SY:                                    Functions and Variables for qm.
                                                              (line 443)
* Sz:                                    Functions and Variables for qm.
                                                              (line 404)
* SZ:                                    Functions and Variables for qm.
                                                              (line 448)
* tpbra:                                 Functions and Variables for qm.
                                                              (line 689)
* tpbraket:                              Functions and Variables for qm.
                                                              (line 703)
* tpbrap:                                Functions and Variables for qm.
                                                              (line 699)
* tpket:                                 Functions and Variables for qm.
                                                              (line 683)
* tpketp:                                Functions and Variables for qm.
                                                              (line 695)
* UU:                                    Functions and Variables for qm.
                                                              (line 658)
* xm:                                    Functions and Variables for qm.
                                                              (line 334)
* xp:                                    Functions and Variables for qm.
                                                              (line 331)
* ym:                                    Functions and Variables for qm.
                                                              (line 340)
* yp:                                    Functions and Variables for qm.
                                                              (line 337)
* zm:                                    Functions and Variables for qm.
                                                              (line 328)
* zp:                                    Functions and Variables for qm.
                                                              (line 325)

* Menu:

* hbar:                                  Functions and Variables for qm.
                                                              (line 125)

