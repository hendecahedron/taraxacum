(ns ^{:doc "Geometric Algebra (first sketch)"
      :author "Matthew Chadwick"}
  hendecahedron.taraxacum.core
  (:require
    [clojure.string :as string :refer [starts-with?]]
    [clojure.math :as maths :refer [pow sqrt signum]]
    [clojure.math.combinatorics :as x]
    [clojure.walk :as w]))

(def b& bit-and)
(def b| bit-or)
(def b< bit-shift-left)
(def b> bit-shift-right)
(def b>> unsigned-bit-shift-right)
(def b¬ bit-not)
(def b⊻ bit-xor)

(defrecord Basis [bitmap scale grade name])

; Java's method reimplemented for CLJS
(defn bit-count [i]
  (let [i (- i (b& (b>> i 1) 0x5555555555555555))
        i (+ (b& i 0x3333333333333333) (b& (b>> i 2) 0x3333333333333333))
        i (b& (+ i (b>> i 4)) 0x0f0f0f0f0f0f0f0f)
        i (+ i (b>> i 8))
        i (+ i (b>> i 16))
        i (+ i (b>> i 32))
        ]
    (b& i 0x7f)))

; grade is the number of factors, or elements
; which have been wedged together - also the
; dimensionality of the subspace represented
; by the k-basis/vector
(def grade bit-count)

; make a basis. 'G' is like a circle with an arrow
; which is like the circular bivector visualization
; in GA4CS, and 'G' is also for Geometric Number
(defn G
  ([from scale]
    (G {} (or (:bitmap from) from) scale (:name from)))
  ([ga bitmap scale]
    (Basis. bitmap scale (grade bitmap) (get-in ga [:basis-by-bitmap bitmap])))
  ([ga bitmap scale name]
    (Basis. bitmap scale (grade bitmap) name)))

(defn multivector [{basis :basis :as ga} elements]
  (mapv
    (fn [[co b]]
      (G (get basis b b) co))
    (partition 2 elements)))

(defn <-- [{:keys [scale grade] :as basis}]
  (let [n (* 0.5 grade (dec grade))
        r (if (== 0 n) 1.0 (last (take n (cycle [-1.0 1.0]))))]
    (G basis (*' scale r))))

(defn <- [multivector]
  (mapv <-- multivector))

(defn involute [{:keys [scale grade] :as basis}]
  (G basis (* scale (pow -1 grade))))

(defn <_
  ([ga mv] (<_ mv))
  ([multivector]
    (mapv involute multivector)))

(defn negate [{:keys [scale grade] :as basis}]
  (G basis (* scale -1N)))

; also called conjugation
(defn negate-mv
  ([ga mv] (negate-mv mv))
  ([mv] (mapv negate mv)))

(defn inverse [{{:syms [• *]} :ops {S 'S} :specials :as ga} mv]
  (let [r (<- mv)
        [{s :scale}] (• r r)]
     (if s
       (* r [(G S (/ 1N s))])
       (throw
         (ex-info (str "non-invertable multivector ["
           (string/join " " (map (fn [{:keys [scale name]}] (str scale name)) mv)) "]") {:non-invertable mv})))))

(defn length
  ([{{:syms [•]} :ops :as ga} mv]
    (length ga mv 16))
  ([{{:syms [•]} :ops :as ga} mv n]
    (if (seq mv)
       (let [[{l :scale}] (• mv mv)] (if l (sqrt (abs l)) 0))
       0)))

(defn normalize [{{:syms [•]} :ops :as ga} mv]
  (if (seq mv)
    (let [l (length ga mv) d (if (== 0 l) 1 (/ 1 l))]
      (mapv (fn [e] (update e :scale * d)) mv))
     mv))

(defn scale [ga mv s]
  (mapv (fn [b] (update b :scale * s)) mv))

; todo check that the bitmaps made by xoring here are < count bases
; also this is only good for small numbers of dimensions
; soon need to work out how to manage large spaces
(defn bases-of
   ([n] (bases-of "e" 0 n))
   ([prefix n] (bases-of prefix 0 n))
   ([prefix start n]
     (reduce
       (fn [r components]
         (let [n (symbol (str prefix (string/join "" (map (fn [[i]] (+ start i)) components))))]
          (assoc r n (G {} (reduce b⊻ (map (comp :bitmap last) components)) 1.0 n))))
       {(symbol (str prefix "_")) (G {} 0 1.0 (symbol (str prefix "_")))}
       (let [b (map (fn [i] (let [n (symbol (str prefix i))] [i n (G {}  (b< 1 i) 1.0 n)])) (range n))]
          (mapcat (fn [k] (x/combinations b k))
            (range 1 (inc n)))))))

(defn bit-flips [a b]
  (->> (b> a 1)
    (iterate (fn [x] (b> x 1)))
    (take-while (fn [x] (> x 0)))
    (map (fn [x] (bit-count (b& x b))))
    (reduce +)))

; page 514 GA4CS
(defn canonical-order [a b]
  (if (== 0 (b& (bit-flips a b) 1)) +1.0 -1.0))

(defn <=grade [max]
  (filter (fn [{g :grade}] (<= g max))))

(defn consolidate-basis [ga]
  (comp
    (partition-by :bitmap)
    (map (fn [[fb :as basis]] (G fb (reduce + (map :scale basis)))))))

(def int-xf (filter (fn [[{ag :grade} {bg :grade} {pg :grade}]] (== pg (- bg ag)))))
(def ext-xf (filter (fn [[{ag :grade} {bg :grade} {pg :grade}]] (== pg (+ ag bg)))))

; this is temporary I promise myself
(defn tiny? [x]
  (< (abs x) 1e-6))

(defn approx-same [a b]
  (every? tiny? (map - (map :scale a) (map :scale b))))

(def remove-scalars-xf (remove (comp zero? :grade)))

(def remove0s (remove (comp tiny? :scale)))

(defn consolidate&remove0s [ga]
  (comp
    (consolidate-basis ga)
    remove0s))

(defn consolidate [mv]
  (vals
    (reduce
      (fn [r {:keys [bitmap scale] :as b}]
        (update r bitmap (fn [x] (if x (update x :scale + scale) b))))
      {} mv)))

(defn simplify- [xf mv]
  (into [] xf mv))

(defn simplify
  ([ga mv]
   (into [] remove0s (consolidate mv))))

(defn simplify0 [ga mv]
  (consolidate mv))

(defn <>
  "return the multivector by grade"
  [mv]
  (into {} (map (juxt :grade identity) mv)))

(defn <>r
  "return the grade r part of mv"
  [r mv]
  ((<> mv) r))

(defn **
  {:doc ""}
  [{{:syms [*]} :ops [e_] :basis-by-grade :as ga} mva mvb]
  (for [a mva b mvb]
    [a b (* a b)]))

(defn ⌋
  {:doc "Left contraction" :ref "§2.2.4 eq 2.6 IPOGA"}
  [{{:syms [*]} :ops [e_] :basis-by-grade :as ga} mva mvb]
  (simplify-
    (comp
      (filter (fn [[{ag :grade} {bg :grade} {g :grade}]] (== g (- bg ag))))
      (map peek)
      remove0s)
    (** ga mva mvb)))

(defn ⌊
  {:doc "Right contraction" :ref "§2.2.4 eq 2.7 IPOGA"}
  [{{:syms [*]} :ops [e_] :basis-by-grade :as ga} mva mvb]
  (simplify-
    (comp
      (filter (fn [[{ag :grade} {bg :grade} {g :grade}]] (== g (- ag bg))))
      (map peek)
      remove0s)
    (** ga mva mvb)))

; left contraction or inner product §3.5.3
(defn ⌋• [{{:syms [*]} :ops [e_] :basis-by-grade :as ga} mva mvb]
  (simplify ga
    (map
       (fn [[{ag :grade :as a} {bg :grade :as b} {g :grade :as p}]]
         (if (== g (abs (- bg ag))) p (G e_ 0)))
       (for [{ag :grade :as a} mva {bg :grade :as b} mvb :when (and (> ag 0) (> bg 0))]
         [a b (* a b)]))))

(defn ⌋' [{{:syms [*]} :ops [e_] :basis-by-grade :as ga} mva mvb]
  (mapv peek
    (filter
       (fn [[{ag :grade :as a} {bg :grade :as b} {g :grade :as p}]]
         (== g (- bg ag)))
       (for [{ag :grade :as a} mva {bg :grade :as b} mvb :when (not= ag bg)]
         [a b (* a b)]))))

(defn basis-range
  "select bases having basis in range f t "
  [mv f t]
  (vec (filter
         (fn [{b :bitmap}]
           ((into (hash-set) (map (fn [i] (int (pow 2 i))) (range f t))) b)) mv)))

(defn imv [mvs]
  (mapv (fn [i mv] (mapv (fn [j e] (G e (if (== i j) 1 0))) (range) mv)) (range) mvs))

(defn point [{{:syms [∼ ∨]} :ops bbg :basis-by-grade {z0 'z0} :specials :as ga} components]
  (∼ (into (mapv (fn [p b] (G b p)) components (rest bbg)) [z0])))

(defn sandwich [{{:syms [* - ⁻]} :ops :as ga} h x]
  (* (- h) x (⁻ h)))

(defn qr
  "a GA implementation of QR decomposition by Householder reflection"
  ([{{:syms [+ - ⁻ * *- *0 • ∧ V ∼ • ⍣ ⃠]} :ops
    [e_ e0 :as bg] :basis-by-grade :as ga} [fmv :as mvs]]
   (loop [n (count fmv) d 0 r mvs q identity qs []]
     (if (< d (dec n))
       (let [
              vd  (mapv (fn [b i] (if (< i d) (assoc b :scale 0) b)) (r d) (range)) ; dth basis vector, zeroed out up to d
              ed  [(update (vd d) :scale (fn [x] (let [sn (signum x)] (* -1.0 (if (== 0 sn) 1.0 sn)))))]
              bi' (+ (⃠ vd) ed)                         ; bisector of unit v and ei
              bi  (if (seq bi') bi' ed)                 ; if v is ei then bisector will be empty
              hy  (∼ bi)                                ; reflection hyperplane
              qd  (fn [x] (* (- hy) x (⁻ hy)))
              qs' (into (vec (repeat d identity)) (repeat (clojure.core/- n d) qd))
            ]
         (recur n (inc d)
           (mapv (fn [f x] (f x)) qs' r)
           (comp q qd)
           qs'))
       {:q (mapv (fn [v] (basis-range (q v) 0 n)) (imv mvs))
        :qfn (fn [mvs] (mapv (fn [v] (basis-range (q v) 0 n)) mvs))
        :r (mapv (fn [v] (basis-range v 0 n)) r)}))))

(defn eigendecompose
  ([{mm :metric-mvs mmga :mmga :as ga}]
    (eigendecompose mmga mm))
  ([ga mvs]
    (eigendecompose ga (qr ga mvs) mvs))
  ([ga {:keys [q r]} mvs]
   {
    :eigenvectors q
    :eigenvalues (vec (map-indexed (fn [i mv] (mv i)) r))
    }))

; this will use the metric of the given GA so that
; will need to be identity
(defn ->basis
  "Change of basis. Based on the GA4CS reference impl"
  ([{basis :eigenvectors mmga :mmga :as ga} basis]
    (->basis mmga basis basis))
  ([{[e_] :basis-by-grade bbb :basis-in-order
    {* '* ∧ '∧ *- '*- •∧ '•∧} :ops :as ga}
    metric-mvs {:keys [bitmap scale] :as basis}]
   (loop [r [(G e_ scale)] i 0 b bitmap]
     (if (== b 0)
       r
       (if (even? b)
         (recur r (inc i) (b>> b 1))
         (recur
           (reduce
             (fn [rr [mv j]]
               (let [s (:scale (mv i))]
                 (if (== s 0)
                   rr
                   (reduce
                     (fn [rrr x]
                       (into rrr (∧ [x] [(G (bbb (b< 1 j)) s)])))
                     rr r))))
             [] (map vector metric-mvs (range)))
           (inc i) (b>> b 1)))))))


(defn op-error
  ([op {help :help :as ga} a b]
   (throw (ex-info (str op " (" (help op) ") can't take " (or (type a) "nil") a " & " (or (type b) "nil") b) {:op op :args [a b]})))
  ([op {help :help :as ga} a]
   (throw (ex-info (str op " (" (help op) ") can't take " (or (type a) "nil") a) {:op op :arg a}))))

(defn compare-G
  ([op]
    (fn dispatch
      ([{ops :ops :as m} a]
        ((ops (compare-G op m a) (partial op-error op)) (assoc m :op op) a))
      ([{ops :ops :as m} a b]
        ((ops (compare-G op m a b) (partial op-error op)) (assoc m :op op) a b))
      ([{ops :ops :as ga} a b & more]
       (cond
         (ops [op :multivectors]) ((ops [op :multivectors]) ga (cons a (cons b more)))
         :default (reduce (partial (ops (compare-G op ga a b) (partial op-error op)) (assoc ga :op op)) (cons a (cons b more)))))))
  ([op ga a b]
    (let [meet (if (and (:bitmap a) (:bitmap b)) (b& (:bitmap a) (:bitmap b)) nil)
          dependency (if (and meet (zero? meet)) :independent :dependent)
          ag (if (vector? a) :grades (min 1 (or (:grade a) 0)))
          bg (if (vector? b) :grades (min 1 (or (:grade b) 0)))
          ta (cond (number? a) :number (:bitmap a) :basis (vector? a) :multivector)
          tb (cond (number? b) :number (:bitmap b) :basis (vector? b) :multivector)
          ]
       [op dependency ta tb ag bg]))
   ([op ga a]
    (cond
      (seq? a) [op :multivectors]
      (vector? a) [op :multivector]
      (:grade a) [op :basis]
      :else [op (type a)])))

(defn ga-ops []
  {
   [:no-such-op :independent :basis :basis 0 0]
   (fn g* [ga a b]
     (update a :scale * (:scale b)))

   ['* :dependent :number :number 0 0]
   (fn g* [ga a b]
     (* a b))

   ['* :independent :basis :basis 0 0]
   (fn g* [ga a b]
     (update a :scale * (:scale b)))

   ['* :independent :basis :basis 0 1]
   (fn g* [ga a b]
     (update b :scale * (:scale a)))

   ['* :independent :basis :basis 1 0]
   (fn g* [ga a b]
     (update a :scale * (:scale b)))

   ['* :independent :basis :basis 1 1]
   (fn g* [ga {bma :bitmap va :scale :as a} {bmb :bitmap vb :scale :as b}]
     (G ga (b⊻ bma bmb)
       (* (canonical-order bma bmb) va vb)))

   ['* :dependent :basis :basis 1 1]
   (fn g* [{metric :metric :as ga} {bma :bitmap va :scale :as a} {bmb :bitmap vb :scale :as b}]
     (G ga (b⊻ bma bmb)
       (loop [i 0 m (b& (:bitmap a) (:bitmap b))
              s (* (canonical-order bma bmb) va vb)]
         (if (== 0 m)
           s
           (recur (inc i) (b> m 1)
             (if (== 1 (b& m 1)) (* s (metric i)) s))))))

   ^{:doc "raw Geometric product"}
   ['*'' :dependent :multivector :multivector :grades :grades]
   (fn g*'' [{{* '*} :ops :as ga} mva mvb]
     (for [a mva b mvb] [a b (* a b)]))

   ; see 22.3.1 for future optimizations
   ^{:doc "Geometric product" :e.g. '(* [1 v1 2 v2] [3 v0 4 v2])}
   ['* :dependent :multivector :multivector :grades :grades]
   (fn g* [{{* '*} :ops :as ga} mva mvb]
     (simplify ga (for [a mva b mvb] (* ga a b))))

   ^{:doc "unsimplified Geometric product"}
   ['*- :dependent :multivector :multivector :grades :grades]
   (fn g*- [{{* '*} :ops :as ga} a b]
     (vec (for [a a b b] (* a b))))

   ^{:doc "simplified Geometric product keep zeros"}
   ['*0 :dependent :multivector :multivector :grades :grades]
   (fn g*0 [{{* '*} :ops :as ga} a b]
     (simplify0 ga (for [a a b b] (* a b))))

   ^{:doc "Interior and exterior products"}
   ['•∧ :dependent :multivector :multivector :grades :grades]
   (fn ip [{{:syms [*'']} :ops :as ga} mva mvb]
     (let [gp (*'' mva mvb)]
       {:• (into [] remove0s (consolidate (simplify- (comp int-xf (map peek)) gp)))
        :∧ (into [] remove0s (consolidate (simplify- (comp ext-xf (map peek)) gp)))}))

   ^{:doc "Interior product ·"}
   ['⌋• :dependent :multivector :multivector :grades :grades]
   (fn left-contraction [{{:syms [•∧]} :ops :as ga} a b]
     (⌋• ga a b))

   ^{:doc "Interior product"}
   ['• :dependent :multivector :multivector :grades :grades]
   (fn ip [{{:syms [•∧]} :ops :as ga} a b]
     (:• (•∧ a b)))

   ^{:doc "Exterior product or meet, largest common subspace, intersection"}
   ['∧ :dependent :multivector :multivector :grades :grades]
   (fn ∧ [{{:syms [•∧]} :ops :as ga} a b]
     (:∧ (•∧ a b)))

   ^{:doc "Regressive product or join, smallest common superspace, union"
     :note "Gunn arXiv:1501.06511v8 §3.1"}
   ['∨' :dependent :multivector :multivector :grades :grades]
   (fn ∨' [{{:syms [* • ∼]} :ops :as ga} a b]
     (∼ (* (∼ a) (∼ b))))

   ^{:doc "Regressive product or join, smallest common superspace, union"
     :note "Gunn arXiv:1501.06511v8 §3.1"}
   ['∨ :dependent :multivector :multivector :grades :grades]
   (fn ∨ [{{:syms [∧ ∼]} :ops {:keys [I I-]} :specials :as ga} a b]
     (simplify ga (∼ (∧ (∼ b) (∼ a)))))

   ; remove this impl
   ^{:doc "remove"}
   ['∨ :multivectors]
   (fn ∨ [{{:syms [∧ ∼]} :ops {:keys [I I-]} :specials :as ga} mvs]
     (let [r (simplify ga (∼ (reduce ∧ (map ∼ mvs))))]
       (if (odd? (count mvs))
         r
         (mapv (fn [b] (update b :scale * -1)) r))))

   ^{:doc ""}
   ['h∨ :dependent :multivector :multivector :grades :grades]
   (fn h∨ [{{:syms [* • ∼]} :ops :as ga} a b]
     ; Hestenes (13) defines ∨ as (• (∼ a) b) which doesn't give the same result
     (• (∼ a) b))

   ^{:doc "Sum, bisect 2 planes, bisect 2 normalized lines"
     :ascii '+ :short 'sum :verbose 'sum :gs '+}
   ['+ :dependent :multivector :multivector :grades :grades]
   (fn g+ [ga a b]
     (simplify ga (into a b)))

   ^{:doc "Sum, bisect 2 planes, bisect 2 normalized lines"
     :ascii '+ :short 'sum :verbose 'sum :gs '+}
   ['+ :multivector]
   (fn g+ [ga a] a)

   ^{:doc "Exponential"
     :ascii 'e :short 'exp :verbose 'exponential}
   ['𝑒 :multivector]
   (fn exp [{{:syms [• *]} :ops
             basis :basis [G_] :basis-in-order :as ga} a]
     (let [[{max :scale}] (• a (<- a))
           scale (loop [s (if (> max 1) (b< 1 1) 1) m max]
                   (if (> m 1) (recur (b< s 1) (/ m 2)) s))
           scaled (* a [(G G_ (/ 1 scale))])
           r (simplify ga (reduce
                            (fn exp1 [r i]
                              (into r (* [(peek r)] (* scaled [(G G_ (/ 1 i))]))))
                            [G_] (range 1 16)))
           ]
       (loop [r r s scale]
         (if (> s 1)
           (recur (* r r) (b>> s 1))
           (simplify ga r)))))

   ^{:doc "Sandwich product"}
   ['|*| :dependent :multivector :multivector :grades :grades]
   sandwich

   ^{:doc "Inverse"}
   ['⁻ :multivector]
   inverse

   ^{:doc "Inverse"}
   ['⁻ :basis]
   (fn [ga b] (<-- b))

   ^{:doc "Negate"}
   ['- :multivector]
   negate-mv

   ^{:doc "Involution"}
   ['_ :multivector]
   <_

   ^{:doc "Dual"}
   ['∼ :multivector]
   (fn dual [{{• '•} :ops duals :duals ds :duals- :as ga} mv]
     ; (⌋ ga a [I-])
     (mapv
       (fn dual-component [{bm :bitmap s :scale :as a}]
         (when (nil? (duals (G a 1.0)))
           (println "  no dual of " a " in this algebra"))
         (assoc (duals (G a 1.0)) :scale (* (ds (G a 1.0)) s))) mv))

   ^{:doc "Dual"}
   ['∼ Basis]
   (fn dual [{{⌋ '⌋ • '•} :ops duals :duals ds :duals- :as ga} {bm :bitmap s :scale :as a}]
     (assoc (duals (G a 1.0)) :scale (* (ds (G a 1.0)) s)))

   ^{:doc "Hodge dual ★"}
   ['★ :multivector]
   (fn hodge [{{* '*} :ops {I 'I} :specials :as ga} mv]
     (* (<- mv) [I]))

   ^{:doc "Hodge dual ★"}
   ['★ Basis]
   (fn hodge [{{* '*} :ops {I 'I} :specials :as ga} x]
     (* (<-- x) I))

   ^{:doc "Normalize"}
   ['⃠ :multivector]
   normalize
   })

; p.80 the inverse of the pseudoscalar is simply the reverse
(defn with-specials
  "Specials for want of a better name"
  [{b :basis-by-grade bio :basis-in-order md :metric :as ga}]
  (let [
        I (peek b)
        I- (<-- I)
        S (first b)
        ]
    (assoc ga :specials
      (reduce
        (fn [r [j [i m]]]
          (if (== 0 m)
            (assoc r (symbol (str "z" j)) (b (inc i)))
            r))
        {'I I 'I- I- 'S S} (map-indexed vector (filter (comp zero? peek) (map-indexed vector md)))))))

; the sign of the dual such that (= I (* x (∼ x)))
(defn with-duals [{bbg :basis-by-grade duals :duals {:syms [*]} :ops :as ga}]
  (assoc ga :duals-
    (into {} (map (fn [[a b]] [a (:scale (* ga a b))]) duals))))

(defn compare-basis [{ag :grade ab :bitmap :as a} {bg :grade bb :bitmap :as b}]
  (if (== ag bg)
    (compare ab bb)
    (compare ag bg)))

(defn with-help [{ops :ops :as a-ga}]
  (-> a-ga
    (assoc :help
      (into {} (map (juxt first (comp :doc meta)) (filter meta (keys ops)))))
    (assoc :examples
      (into {} (map (juxt first (comp :e.g. meta)) (filter meta (keys ops)))))))

(defn with-eigendecomposition [{mm :metric-mvs :as ga}]
  (if mm
    (let [{:keys [eigenvalues eigenvectors]} (eigendecompose ga)]
      (assoc ga :eigenvalues eigenvalues :eigenvectors eigenvectors))
    ga))

(defn ga-
  "
   Create a new GA from the given params:

   md - metric diagonal

   Metrics:

   Euclidean - all diagonal elements are 1
   Diagonal - metric factors only along diagonal
   Othonormal -

   use :pqr to permute pqr e.g. [:r :q :p] so e0 is 0^2
  "
  ([{:keys [prefix base p q r pm qm rm md pqr]
     :or {prefix "e" base 0 pm 1 qm -1 rm 0 pqr [:p :q :r]}}]
   (let [md (or md (vec (apply concat (map {:p (repeat p pm) :q (repeat q qm) :r (repeat r rm)} pqr))))
         p (or p (count (filter (partial == 1) md)))
         q (or q (count (filter (partial == -1) md)))
         r (or r (count (filter (partial == 0) md)))
         d (+ p q r)
         bases (bases-of prefix base d)
         bio (reduce (fn [r [n b]] (assoc r (:bitmap b) b)) (vec (repeat (count bases) 0)) bases)
         bbg (vec (sort compare-basis (vals bases)))
         zv (vec (repeat (count bases) 0))
         m {
            :p p :q q :r r
            :metric md
            :basis bases
            :dimensions d
            :size (int (pow 2 d))
            :zv zv
            :basis-by-bitmap (reduce (fn [r [n b]] (assoc r (:bitmap b) n)) zv bases)
            :duals (zipmap bbg (rseq bbg))
            :basis-by-grade bbg
            :basis-in-order bio
            :ops (ga-ops)
            }
         ; note ops must be in order of dependence because of the partial later
         ops '[+ * ⍣ - _ *'' *- *0 •∧ • •' ⁻ ∧ ∼ ∨ ∨' h∨ ★ ⃠ 𝑒 op]
         ]
     (ga- m ops)))
  ([m ops]
   (ga- m
     (reduce
       (fn [r op] (assoc-in r [:ops op] (compare-G op))) m ops) ops))
  ([m g ops]
   (reduce
     (fn [r op]
       (update-in r [:ops op]
         (fn [g] (partial g r))))
     (-> g with-duals with-specials with-help) ops)))

(defn ga
  ([p q r]
   (ga "e" p q r))
  ([prefix p q r]
   (ga {:prefix prefix :p p :q q :r r}))
  ([prefix base p q r]
   (ga {:prefix prefix :base base :p p :q q :r r}))
  ([{:keys [prefix base p q r pm qm rm md mm mmga]
     :or {prefix "e" base 0 pm 1 qm -1 rm 0} :as params}]
   (let [
         ; if creating a GA from a metric-multivectors, precompute their eigendecomposition
         ; now for later changes of basis. Such a change of basis needs a GA itself,
         ; so either use the one supplied or create one
         {:keys [eigenvalues eigenvectors]} (if mm (eigendecompose mmga mm) nil)
         g (ga- (if mm (assoc params :md (mapv :scale eigenvalues)) params))
         g' (assoc g :eigenvectors eigenvectors :eigenvalues eigenvalues
              :metric-mvs mm :mmga (or mmga (ga- {:prefix 'b :p (count mm) :q 0 :r 0})))
         ]
     g')))

(defn multivector? [f]
  (and
    (or
      (and (list? f) (every? (fn [[c b]] (and (number? c) (and (symbol? b) (nil? (namespace b))))) (partition 2 f)))
      (and (vector? f) (every? (fn [[c b]] (and (or (number? c) (list? c)) (symbol? b))) (partition 2 f))))
    (== 0 (mod (count f) 2))))

(defn basislike [x]
  (if (symbol? x)
    (let [basis (namespace x) n (name x)]
      (try {:basis (symbol basis) :scale (parse-double n)}
        (catch Exception e nil)))
    nil))

(defn multivector1? [f]
  (and (list? f)
    (every? (fn [x] (basislike x)) f)))

; TODO: function names beginning with 'e' don't work
(defmacro in-ga
  ([prefix p q r body]
   `(in-ga {:prefix ~prefix :base 0 :p ~p :q ~q :r ~r} ~body))
  ([prefix base p q r body]
   `(in-ga {:prefix ~prefix :base ~base :p ~p :q ~q :r ~r} ~body))
  ([p q r body]
   `(in-ga {:prefix e :base 0 :p ~p :q ~q :r ~r} ~body))
  ([{:keys [prefix base p q r mm pqr] :or {base 0 prefix 'e pqr [:p :q :r]}} body]
   (let [
         prefix (name prefix)
         ops '[+ * 𝑒 ⍣ - _ *'' *- *0 •∧ • •' ⁻ ∧ ∼ ∨ ∨' h∨ ★ ⃠ op]
         specials '[I I- S]
         opz (into #{} ops)
         o (complement opz)
         s (filter symbol? (tree-seq seqable? seq body))
         e (vec (cons (symbol (str prefix "_")) (filter (fn [i] (or (starts-with? (name i) prefix)
                                                                  (and (namespace i) (starts-with? (namespace i) prefix)))) s)))
         e (mapv (fn [x] (if-let [{b :basis} (basislike x)] b x)) e)
           ]
       `(let [{{:syms ~e :as ~'basis} :basis
               {:syms ~ops} :ops
               ~'basis-by-grade :basis-by-grade
               ~'basis-by-bitmap :basis-by-bitmap
               ~'basis-in-order :basis-in-order
               ~'basis :basis
               ~'duals :duals
               ~'duals- :duals-
               {:syms ~specials} :specials :as ~'-ga} (ga {:prefix ~prefix :base ~base :p ~p :q ~q :r ~r :mm ~mm :pqr ~pqr})]
          ~(w/postwalk
             (fn [f]
               (cond
                 (and (multivector? f) (o (first f)))
                   (mapv (fn [[c b]] `(G ~b ~c)) (partition 2 f))
                 (and (multivector1? f) (o (first f)))
                   (mapv (fn [x] (let [{b :basis c :scale}
                                       (if (number? x)
                                         `{:basis ~(symbol (str prefix "_")) :scale ~x}
                                         (basislike x))]
                                   `(G ~b ~c))) f)
                 :else f))
               body)))))


(defmethod print-method Basis [{:keys [name scale]} writer]
  (doto writer
    (.write (str scale))
    (.write " ")
    (.write (str name))))