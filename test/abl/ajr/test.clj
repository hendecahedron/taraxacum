(ns abl.ajr.test
  (:require
    [clojure.test :refer :all]
    [abl.ajr.core :refer :all]))

(require :reload 'abl.ajr.core)

(comment


  (in-ga 4 0 0 (∨ [(/ 1 2) e12] [2 e23]))

  (b⊻ 1 2)


  (in-ga 4 0 0 (∼ [1 e1 2 e2]))

  (in-ga 4 0 0 (⍟ [1 e1 2 e2]))

  (in-ga 4 0 0 (∧  [0.1 e0 0.2 e1 0.3 e2 -1.5 e3]
                   [0.7 e0 0.4 e12 0.1 e2 -1 e13]))



  (in-ga 4 0 0 (•∧ [0.1 e0 0.2 e1 0.3 e2 -1.5 e3]
                   [0.7 e0 0.4 e1 0.1 e2 -1 e3]))

  (in-ga 4 0 0 (* [0.1 e0 0.2 e1 0.3 e2 -1.5 e3]
                  [0.7 e0 0.4 e1 0.1 e2 9e-7 e3]))

  (in-ga 2 0 0 (* e0 e1))
  (in-ga 2 0 0 (* e1 e0))

  (in-ga 4 0 0 (* [0.1 e0 0.2 e1 0.3 e2 -1.5 e3] [I]))

  (print-cause-trace *e)

  (let [{{:syms [e_ e0 e1 e2 e3 e12 e01 e02 e03 e13 e23 e012 e023 e123 e0123]} :basis
       {:syms [* • ∼ v - _ s]} :ops
       {:syms [I I-]} :specials
         :as ga} (ga 'e 3 0 0)
     ]
  (qr ga
    (map (partial multivector ga)
      [[0.1 e0 0.4 e1 0.4 e2]
       [0.4 e0 0.2 e1 0.4 e2]
       [0.4 e0 0.3 e1 0.2 e2]])))

  (let [{{:syms [e_ e0 e1 e2]} :basis
       {:syms [* • ∧]} :ops} (ga 'e 3 0 0)]
    [(* (G e1 0.4) (G e2 0.7))
     (* (G e2 0.7) (G e1 0.4))])

  (G )

  (in-ga 4 0 0
    (basis-range 2 3
     [[0 e0 0 e1 0 e2 -1 e3]
      [0 e0 1 e1 0 e2 0 e3]
      [0 e0 0 e1 1 e2 0 e3]
      [-1 e0 0 e1 0 e2 0 e3]]))

  ; could also select bases like this
  (in-ga 4 0 0
    (simplify ga
      (map *
       [0 e0 1 e1 1 e2 1 e3]
       [1 e0 -1 e1 1 e2 -1 e3]
       [1 e0 1 e1 1 e2 1 e3])))

  ; or as bivectors
  (in-ga 4 0 0
    [[0 e0 0 e1 0 e2 -1 e03]
     [1 e03]])

  (in-ga 3 0 1
    (+ [0.3 e01] [0.3 e02]))

  (in-ga 3 0 1 (:basis-by-bitmap ga))

  ; in each iteration:
  ;
  ; select bases range k n
  ; take first mv v
  ; find the hyperplane h which reflects v to xe0:
  ;   find bisector of v and e0, that defines h
  ; reflect all mvs in h


  (in-ga 3 0 0 (* [1 e01] [0.5 e01 0.5 e12]))

  (in-ga e 1 3 0 1 (* e12 e12))


  (bases-of 'e 1 4)

  (macroexpand '(in-ga e 3 0 1 (* (e0/+0.1e-3) [4 e1/-0.6 e12/-0.0134])))

  (in-ga e 3 0 1
    (reduce * [ [5 e2 4 e1] [e0/+0.1e-3] (4 e1/-0.6 e12/-0.0134)]))

  (in-ga e 3 0 1 (let [x 3] (* (x e0) (2 e1))))



  (every? (fn [[c b]] (number? c) (symbol? b)) (partition 2 '[a e1]))
 )



(comment



:e3/+0.7

e3/+0.7


; how to manage high dimensional algebras ? 131010
(count (bases-of (+ 17 0 0)))

; Algebra(0,1,()=>(3+2e1)*(1+4e1));

 (in-ga 0 1 0 (* [3 e_ 2 e0] [1 e_ 4 e0]))

; Algebra(0,1,() => (3 + 2e1) * (1 + 4e1));

  ; 2,3i

  ; 3e_

  (in-ga 0 1 0
    (reductions * (repeat 4 (0 e0/+1))))

  (in-ga e 2 0 1 (* [e_/+3 e0/-2e-19] [e_/+1 e0/+4]))

  (in-ga e 7 1 0 (* [3 e0] [5 e0] [5 e0] [5 e0]))

  (in-ga e 3 0 1 (* [e12] (∼ [e12] [e2])))

  (in-ga e 4 0 0 (* [e123] [e012]))


  (in-ga e 3 0 0 (• (∼ [e12]) [e01]))

  (in-ga e 3 0 0 (V [e01] [e12]))

  (in-ga e 4 0 0 (* I I-))

  (in-ga e 1 3 0 0 (* e12 e123 e23 e123 e123))


  (in-ga e 4 0 0 (V [e2 e3] (* (<- [e1 e2]) [I])))

 (in-ga e 2 0 1 (* [3 e12, 2e-19 e0] [1.3e-5, e2 4 e0]))

 (in-ga v 2 0 1 (* [3 v1 2 v0] [1.3e-5 v1 4 v0]))

(in-ga e 3 0 1 (:help ga))


(let [{{:syms [e_ e0 e1 e2 e3 e12 e01 e02 e03 e13 e23 e012 e023 e123 e0123]} :basis
       {:syms [* • ∼ v - _ s]} :ops
       {:syms [I I-]} :specials
         :as ga} (ga 'e 3 0 1)
     ]
  {
   :reflection (reduce * [[e1 e2] (_ [(G e1 0.2) (G e2 0.8)]) (- [e1 e2])])
   :reflection2 (reduce * [(* [e1] [e0]) [(G e1 0.2) (G e2 0.8)] (* [e1] [e0])])
   :dual-e0 (∼ [e1 e2])
   :i (* [e01] [e12])
   :*v (reduce * [[(G e1 0.2) (G e2 0.8)] [(G e1 0.9) (G e2 0.1)]])
   })

  (in-ga 3 0 1
    {
     :reflection (reduce * [[1 e1 1 e2] (_ [0.3 e1 0.8 e2]) (- [1 e1 1 e2])])
     :reflection2 (reduce * [(* [1 e1] [1 e0]) [0.2 e1 0.8 e2] (* [1 e1] [1 e0])])
     :dual-e0 (∼ [1 e1 1 e2])
     :i (* [1 e01] [1 e12])
     :*v (reduce * [[0.2 e1 0.8 e2] [0.9 e1 0.1 e2]])
     })

(in-ga {:prefix v :p 3 :q 0 :r 1}
  {
   :reflection (reduce * [[1 v1 1 v2] (_ [0.3 v1 0.8 v2]) (- [1 v1 1 v2])])
   :reflection2 (reduce * [(* [1 v1] [1 v0]) [0.2 v1 0.8 v2] (* [1 v1] [1 v0])])
   :dual-e0 (∼ [1 v1 1 v2])
   :d1 (* [0.3 v1 0.7 v2] [0.8 v1 0.2 v2])
   :i (* [1 v01] [1 v12])
   :x (* [v0 v1 v2] [v0 v1 v2])
   :p (sort compare-blades (*' [0.1 v0 0.7 v1 0.2 v2] [0.1 v0 0.7 v1 0.2 v2]))
   })

  (:examples (ga {:prefix 'v :p 3 :q 0 :r 0}))

  (:help (ga {:prefix 'v :p 3 :q 0 :r 0}))

  (in-ga v 3 0 0 (:help ga))

  (in-ga e 3 0 1 (* [e0 e1 e2] [e0]))

  (in-ga e 3 0 0
    (• [0.3 e0 0.4 e1 0.5 e2] [0.3 e0 0.4 e1 0.5 e2]))

  (Math/sqrt 19.3)

  (reduce + (map * [0.3 0.4 0.5] [0.3 0.4 0.5]))

  (in-ga e 3 0 1 (* [2 e0 3 e1 1 e2] [1 e0 0.7 e1 0.8 e2]))

  (in-ga e 3 0 1 (* [e3] [e12] (∼ e1)))

  (in-ga {:prefix v :p 3 :q 0 :r 0} (V [0.5 v01 0.5 v12] [1 v12]))

  (in-ga {:prefix v :p 3 :q 0 :r 0} (h [0.5 v01 0.5 v12]))

(:basis-by-grade (ga 3 0 1))

  ; Hestenes - blades which are duals of vectors
  ; are called covectors which define hyperplanes

  (in-ga e 4 0 0 (∼ [2 e0]))

  ; (∼ [2 e0]), [2 e123], is linear functional
  ; representing a 3d hyperplane in R4
  (in-ga e 4 0 0 (• [0.5 e0 0.5 e1] (∼ [2 e0])))
  ; [1 e23] is the interior product which tells us
  ; the amount of [0.5 e0 0.5 e1] that intersects [2 e123]
  ; e23 being the component of e123 not in the query vector

  (in-ga e 3 0 0 (• [0.5 e0 0.5 e1] (∼ [3 e0 2 e1])))

  (in-ga e 3 0 0 (V [0.5 e0 0.5 e1] (- [3 e0 2 e1])))

  (in-ga 9 0 0 (∼ [0 e0 1 e1 2 e2 3 e3 4 e4 5 e5 6 e6 7 e7 8 e8]))

  (in-ga 5 0 0 (∼ [0 e0 1 e1 2 e2 3 e3 4 e4]))


  (/ 1 (Math/sqrt 2))
"

⋆

another way to do interior product

(v a (hodge b))

project point P onto line l  (l.P)l

project line l onto point P  (l.P)P

reflect x in hyperplane normal to h

(* [-1 h] [1 x] [1 h])

An orthogonal metric with respect to the
orthonormal basis vectors ei is one where

 eᵢ * eⱼ = mᵢ𝛿ᵢⱼ


 x  y  z  w

 !  0  0  0
 0  !  0  0
 0  0  !  0
 0  0  0  0


"

()


(let [{{:syms [v_ v3 v02 v01 v12]} :basis
       [v_ v0 e1 v2] :basis-in-order
       {:syms [* • 𝑒 ∼ ⍣]} :ops :as ga} (ga 'v 2 0 1)
       m (fn [a p] (𝑒 (* [(G v_ a)] p)))
       p (fn [x y] [(G v02 x) (G v01 y) v12])
       rpi (fn [x] (/ Math/PI x))
      ]
  (take 8
   (iterate (partial ⍣ (m (rpi 8) (p 0 3)))
     (p 1 1))))


  (in-ga v 2 0 1
    (let [m (fn [a p] (𝑒 (* [a v_] p)))
          p (fn [x y] [x v02 y v01 1 v12])
          rpi (fn [x] (/ Math/PI x))
          ]
      (take 8
        (iterate (partial ⍣ (m (rpi 8) (p 0 3)))
          (p 1 1)))))


  ; a 2-d vector
  (in-ga 3 0 0 [[0.7 e0] [-0.3 e2]])

  ; a bivector in 3d
  (in-ga 3 0 0 [[0.7 e02] [-0.3 e01]])



"

I = e12

Ie1 = −e2

Ie2 = e1

I(−e1) = e2

I(−e2) = −e1

"


; A multivector is a linear combination of diﬀerent k-blades
[(G e0 0.1) (G e1 0.2) (G e3 0.7) (G e23 0.6)]

  ; there would 2ⁿ basis blades where n=sum(p,q,r)



)