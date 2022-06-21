(ns abl.ajr.test
  (:require
    [clojure.stacktrace :refer [print-cause-trace]]
    [clojure.test :refer :all]
    [abl.ajr.core :refer :all]
    [clojure.string :as string]))


(comment

  (print-cause-trace *e)


  (in-ga 3 0 1
    (∧
      [1 e0 2 e1 3 e2 1 e3]
      [4 e0 3 e1 2 e2 1 e3]
      [-1 e0 -7 e1 1 e2 1 e3]))

  ; composition two reflections
  (in-ga 2 0 1
    (*
      [-5 e2 1 e0]
      [-7 e2 1 e1]))

  ; to form a point
  ; => [1.0e01 -7.0e02 5.0e12]

  ; (note the due to the order of basis in
  ; the naming, mine are negative those of Enki's
  ; because he uses 0 for my 2 - my results agree)

  (in-ga 2 0 1
    (*
      (- [1 e01 4 e0 -4 e1])
      [-5 e2 1 e0]
      [-7 e2 1 e1]
      [1 e01 4 e0 -4 e1]))




 )



(comment



:e3/+0.7

'e3/+0.7


; how to manage high dimensional algebras ? 131010
(count (bases-of (+ 17 0 0)))

; Algebra(0,1,()=>(3+2e1)*(1+4e1));

 (in-ga 0 1 0 (* [3 e_ 2 e0] [1 e_ 4 e0]))



(let [{{:syms [v_ v3 v02 v01 v12]} :basis
       [v_ v0 e1 v2] :basis-by-grade
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
        (iterate (partial ⍣ (m (rpi 8.0) (p 0.0 3.0)))
          (p 1.0 1.0)))))



  (in-ga 3 0 0
    (qr ga
      [
        [1 e0 2 e1 3 e2]
        [4 e0 5 e1 6 e2]
        [7 e0 8 e1 9 e2]]))


  (in-ga 3 0 0
    (* 2.0 3.0))

  ; table 5.2 Gunn - correct results
  (in-ga 2 0 1
    (let [E2 e01 E0 e12 E1 e02]
      (* [E2] [e2])))

  (in-ga 3 0 1
    (let [E2 e01 E0 e12 E1 e02]
      (* [E2] [e2])))

  ; ∥ ∥ ∞ ideal norm

  '(("Intersection point of two lines" "a ∧ b")
     ("Angle of two intersecting lines" "cos−1(a · b), sin−1(∥a ∧ b∥)")
     ("Distance of two || lines" "∥a ∧ b∥∞")
     ("Joining line of two points" "P ∨ Q")
     ("⊥ direction to join of two points" "P × Q")
     ("Distance between two points" "∥P ∨ Q∥, ∥P × Q∥∞")
     ("Oriented distance point to line" "∥a ∧ P∥")
     ("Angle of ideal point to line" "sin−1 (∥a ∧ P∥∞)")
     ("Line through point ⊥ to line" "P · a")
     ("Nearest point on line to point" "(P · a)a")
     ("Line through point || to line" "(P · a)P")
     ("Area of triangle ABC" "½∥A∨B∨C∥")
     ("Reflection in line (X = point or line)" "aXa")
     ("Rotation around point of angle 2α" "RXR̃ (R := eαP)")
     ("Translation by 2d in direction V⊥" "TXT̃ (T := 1 + dV)"))


  ; (* (• x a) (⁻ a)) projects x onto a, x and a can be any objects
  ; (* (∧ x a) (⁻ a)) object through x perpendicular to a




 (in-ga 2 0 1
   (for [a basis-by-grade b basis-by-grade]
     [a b (* a b)]))

 (in-ga 3 0 0
   (* [e12] [e012]))

(defn intersect-lines [{{:syms [∼ • ∧]} :ops :as ga} a b]
  (• (∼ a) b))


 (in-ga 3 0 0
   (• (∼ (∧ [2 e0 0 e1 1 e2] [2 e0 4 e1 1 e2]))
      (∧ [0 e0 2 e1 1 e2] [4 e0 2 e1 1 e2])))



  ; 1blade 2blade orthogonal 1blade
  [1.0 e0 1.0 e01 1.0 e1]

  ; 1blade 3blade orthogonal 2blade
  [1.0 e0 1.0 e012 1.0 e12]

  ; 3blade 1blade othogonal 2blade
  [1.0 e012 1.0 e0 1.0 e12]
)