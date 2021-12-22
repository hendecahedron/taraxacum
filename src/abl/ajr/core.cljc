(ns abl.ajr.core
  "

  "
  (:require
    [clojure.string :as string :refer [starts-with?]]
    [clojure.math.combinatorics :as x]
    [clojure.java.math :refer :all :exclude [min max]]
    [clojure.walk :as w])
  (:import
     #?(:clj  [clojure.lang PersistentVector]
        :cljs [cljs.core    PersistentVector])))

(def b& bit-and)
(def b| bit-or)
(def b< bit-shift-left)
(def b> bit-shift-right)
(def b>> unsigned-bit-shift-right)
(def b¬ bit-not)
(def b⊻ bit-xor)


; metric - §19.3
;
; k-vectors (weighted sums of k-blades)
; are represented by vectors of k-blades
; Or, could add a new type and use protocols
;
; ⋀ wedge is also called outer product or exterior product
; -- although some say that outer and exterior are different
; because they produce different types of objects;
; "inner" and "outer" products are associated with
; vectors and matrices while
; interior and exterior products are for GA
; because they produce blades not matrices
;
; the geometric product is interior + exterior products
; exterior product increases grade, interior product decreases grade

(defrecord Blade [bitmap scale grade])

; grade is the number of factors, or elements
; which have been wedged together - also the
; dimensionality of the subspace represented
; by the k-blade/vector
(defn grade [bits]
  (Long/bitCount bits))

; make a blade. 'G' is like a circle with an arrow
; which is like the circular bivector visualization
; in GA4CS, and 'G' is also for Geometric Number
(defn G
  ([basis scale]
    (G {} basis scale))
  ([ga basis scale]
    (let [b (or (:bitmap basis) basis)]
      (assoc (Blade. b (double scale) (grade b))
        :basis
        (or (get-in ga [:basis-by-bitmap b])
          (:basis basis))))))

(defn edalb [{:keys [scale grade] :as blade}]
  (G blade (* scale (pow -1.0 (* 0.5 grade (dec grade))))))

(defn <- [multivector]
  (mapv edalb multivector))

(defn involute [{:keys [scale grade] :as blade}]
  (G blade (* scale (pow -1.0 grade))))

(defn <_
  ([ga mv] (<_ mv))
  ([multivector]
    (mapv involute multivector)))

(defn negate [{:keys [scale grade] :as blade}]
  (G blade (* scale -1.0)))

(defn negate-mv
  ([ga mv] (negate-mv mv))
  ([mv] (mapv negate mv)))

(defn inverse [{{:syms [• *]} :ops {S 'S} :specials :as ga} mv]
  (let [r (<- mv)
        [{s :scale}] (• r r)]
     (if s
       (* r [(G S (/ 1.0 s))])
       (throw (ex-info (str "non-invertable multivector "
                         (string/join " " (map (fn [{:keys [scale basis]}] (str scale basis)) mv))) {:non-invertable mv})))))

(defn ⌋ [{[s] :basis-by-bitmap basis :basis :as ga} {ag :grade :as a} {bg :grade :as b} {rg :grade :as r}]
     (if (or (> ag bg) (not (== rg (- bg ag))))
       (G (basis s) 0) r))

(defn normalize [{{:syms [•]} :ops :as ga} mv]
  (if (seq mv)
    (let [[{l :scale}] (• mv mv)
           d (sqrt l)] (mapv (fn [e] (G e (/ (:scale e) d))) mv))
     mv))

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
          (assoc r n
            (assoc (G (reduce b⊻ (map (comp :bitmap last) components)) 1)
              :basis n))))
       {(symbol (str prefix "_")) (assoc (G 0 1) :basis (symbol (str prefix "_")))}
       (let [b (map (fn [i] [i (symbol (str prefix i)) (G (b< 1 i) 1)]) (range n))]
          (mapcat (fn [k] (x/combinations b k))
            (range 1 (inc n)))))))

(defn bit-flips [a b]
  (reduce +
    (map (fn [x] (Long/bitCount (b& x b)))
      (take-while (fn [x] (> x 0))
        (iterate (fn [x] (b> x 1)) (b> a 1))))))

; page 514 GA4CS
(defn canonical-order [a b]
  (if (== 0 (b& (bit-flips a b) 1)) +1.0 -1.0))

(defn simplify [ga blades]
  (into []
    (comp
      (partition-by :bitmap)
      (map (fn [[fb :as bb]] (G ga fb (reduce + (map :scale bb)))))
      (remove (comp zero? :scale)))
    (sort-by :bitmap blades)))

(defn simplify0 [ga blades]
  (into []
    (comp
      (partition-by :bitmap)
      (map (fn [[fb :as bb]] (G ga fb (reduce + (map :scale bb))))))
    (sort-by :bitmap blades)))


(defn basis-range
  "select blades having basis in range f t"
  [mv f t]
  (vec (filter
         (fn [{b :bitmap}]
           ((into (hash-set) (map (fn [i] (int (pow 2 i))) (range f t))) b)) mv)))

(defn imv [mvs]
  (mapv (fn [i mv] (mapv (fn [j e] (G e (if (== i j) 1 0))) (range) mv)) (range) mvs))


; upper trianguler matrix is one where the basis vectors
; are of increasing dimensionality, nth basis vector is n dimensional
;
; GA4CS§7.1
(defn qr
  ([{{:syms [+ - ⁻ * *' *0 • ∧ V ∼ • ⍣ ⧄]} :ops
    [e_ e0 :as bg] :basis-by-grade :as ga} [fmv :as mvs]]
   (loop [n (count fmv) d 0 r mvs q identity qs []]
     (if (< d (dec n))
       (let [
              v (mapv (fn [b i] (if (< i d) (assoc b :scale 0) b)) (r d) (range)) ; vector
              e [(update (v d) :scale (fn [x] (* -1.0 (signum x))))]
              bi (+ (⧄ v) e)                         ; bisector of unit v and ei
              h (∼ bi)                               ; reflection hyperplane
              qi (fn [x] (*0 (- h) x (⁻ h)))
              qs' (into (vec (repeat d identity)) (repeat (clojure.core/- n d) qi))
            ]
         (recur n (inc d)
           (mapv (fn [f x] (f x)) qs' r)
           (comp q qi) ; check if there's a way to compose hyperplane reflections directly, i.e. compose into one reflection
           qs'))
       {:q (mapv (fn [v] (basis-range (q v) 0 n)) (imv mvs))
        :qf (fn [mvs] (mapv (fn [v] (basis-range (q v) 0 n)) mvs))
        :r (mapv (fn [v] (basis-range v 0 n)) r)}))))

(defn op-error
  ([op {help :help :as ga} a b]
   (throw (ex-info (str "operation: " op " (" (help op) ") can't take " (type a) " & " (type b)) {:op op :args [a b]})))
  ([op {help :help :as ga} a]
   (throw (ex-info (str "operation: " op " (" (help op) ") can't take " (type a)) {:op op :arg a}))))

(defn compare-G
  ([op]
    (fn dispatch
      ([{ops :ops :as m} a]
        ((ops (compare-G op m a) (partial op-error op)) (assoc m :op op) a))
      ([{ops :ops :as m} a b]
        ((ops (compare-G op m a b) (partial op-error op)) (assoc m :op op) a b))
      ([{ops :ops :as m} a b & more]
        (reduce (partial (ops (compare-G op m a b) (partial op-error op)) (assoc m :op op)) (cons a (cons b more))))))
  ([op ga a b]
    (let [meet (if (and (:bitmap a) (:bitmap b)) (b& (:bitmap a) (:bitmap b)) nil)
          dependency (if (and meet (zero? meet)) :independent :dependent)
          ag (if (vector? a) :grades (min 1 (or (:grade a) 0)))
          bg (if (vector? b) :grades (min 1 (or (:grade b) 0)))]
       [op dependency (type a) (type b) ag bg]))
   ([op ga a]
    (cond
      (seq? a) [op :multivectors]
      (vector? a) [op :multivector]
      :else [op (type a)])))

(defn ga-ops []
  {
   [:no-such-op :independent Blade Blade 0 0]
   (fn g* [ga a b]
     (update a :scale * (:scale b)))

   ['* :dependent Double Double 0 0]
   (fn g* [ga a b]
     (* a b))

   ['* :independent Blade Blade 0 0]
   (fn g* [ga a b]
     (update a :scale * (:scale b)))

   ['* :independent Blade Blade 0 1]
   (fn g* [ga a b]
     (update b :scale * (:scale a)))

   ['* :independent Blade Blade 1 0]
   (fn g* [ga a b]
     (update a :scale * (:scale b)))

   ['* :independent Blade Blade 1 1]
   (fn g* [ga {bma :bitmap va :scale :as a} {bmb :bitmap vb :scale :as b}]
     (G ga (b⊻ bma bmb)
       (* (canonical-order bma bmb) va vb)))

   ['* :dependent Blade Blade 1 1]
   (fn g* [{metric :metric :as ga} {bma :bitmap va :scale :as a} {bmb :bitmap vb :scale :as b}]
     (G ga (b⊻ bma bmb)
       (loop [i 0 m (b& (:bitmap a) (:bitmap b))
              s (* (canonical-order bma bmb) va vb)]
         (if (== 0 m)
           s
           (recur (inc i) (b> m 1)
             (if (== 1 (b& m 1)) (* s (metric i)) s))))))

   ; see 22.3.1 for future optimizations
   ^{:doc "Geometric product" :e.g. '(* [1 v1 2 v2] [3 v0 4 v2])}
   ['* :dependent PersistentVector PersistentVector :grades :grades]
   (fn g* [{{* '*} :ops :as ga} a b]
     (simplify ga (for [a a b b] (* ga a b))))

   ^{:doc "unsimplified Geometric product"}
   ['*' :dependent PersistentVector PersistentVector :grades :grades]
   (fn g*' [{{* '*} :ops :as ga} a b]
     (vec (for [a a b b] (* a b))))

   ^{:doc "simplified Geometric product keep zeros"}
   ['*0 :dependent PersistentVector PersistentVector :grades :grades]
   (fn g*0 [{{* '*} :ops :as ga} a b]
     (simplify0 ga (for [a a b b] (* a b))))

   ^{:doc "Hodge dual ⍟"}
   ['⍟ :multivector]
   (fn hodge [{{* '*} :ops {I 'I} :specials :as ga} mv]
     (* (<- mv) [I]))

   ^{:doc "Interior product"
     :ascii 'o :short 'ip :verbose 'interior-product :gs '><}
   ['• :dependent PersistentVector PersistentVector :grades :grades]
   (fn ip [{{:syms [*]} :ops :as ga} a b]
     (simplify ga
       (for [a a b b]
         (⌋ ga a b (* a b)))))

   ^{:doc "Interior and exterior products"
     :ascii 'ox :short 'ox :verbose 'interior-and-exterior-product :gs '><<>}
   ['•∧ :dependent PersistentVector PersistentVector :grades :grades]
   (fn ip [{{:syms [*]} :ops :as ga} ae be]
     (let [ie (vec (for [a ae b be]
                     (let [gp (* a b)]
                       [(⌋ ga a b gp) gp])))
           usi (into #{} (map first ie))
           int (simplify ga usi)
           ext (simplify ga (remove usi (map peek ie)))]
       {:• int :∧ ext}))

   ^{:doc "Exterior product or meet, largest common subspace, intersection"
     :ascii 'x :short 'xp :verbose 'exterior-product}
   ['∧ :dependent PersistentVector PersistentVector :grades :grades]
   (fn ∧ [{{:syms [•∧]} :ops :as ga} a b]
     (:∧ (•∧ a b)))

   ^{:doc "Regressive product or join, smallest common superspace, union"
     :ascii 'v :short 'rp :verbose 'regressive-product}
   ['∨ :dependent PersistentVector PersistentVector :grades :grades]
   (fn ∨ [{{:syms [* • ∼]} :ops :as ga} a b]
     ; Hestenes (13) defines ∨ as (• (∼ a) b) which doesn't give the same result
     (∼ (* (∼ a) (∼ b))))

   ^{:doc "Sum, bisect 2 planes, bisect 2 normalized lines"
     :ascii '+ :short 'sum :verbose 'sum :gs '+}
   ['+ :dependent PersistentVector PersistentVector :grades :grades]
   (fn g+ [ga a b]
     (simplify ga (into a b)))

   ^{:doc "Exponential"
     :ascii 'e :short 'exp :verbose 'exponential}
   ['𝑒 :multivector]
   (fn exp [{{:syms [• *]} :ops
             basis :basis [G_] :basis-by-bitmap :as ga} a]
     (let [[{max :scale}] (• a (<- a))
           scale (loop [s (if (> max 1) (b< 1 1) 1) m max]
                   (if (> m 1) (recur (b< s 1) (/ m 2)) s))
           scaled (* a [(G (basis G_) (/ 1 scale))])
           r (simplify ga (reduce
                            (fn exp1 [r i]
                              (into r (* [(peek r)] (* scaled [(G (basis G_) (/ 1 i))]))))
                            [(basis G_)] (range 1 16)))
           ]
       (loop [r r s scale]
         (if (> s 1)
           (recur (* r r) (b>> s 1))
           (simplify ga r)))))

   ^{:doc "Sandwich product"}
   ['⍣ :dependent PersistentVector PersistentVector :grades :grades]
   (fn |s| [{{* '*} :ops :as ga} r mv]
     (reduce * [(<- r) mv r]))

   ^{:doc "Inverse"}
   ['⁻ :multivector]
   inverse

   ^{:doc "Negate"}
   ['- :multivector]
   negate-mv

   ^{:doc "Involution"}
   ['_ :multivector]
   <_

   ^{:doc "Dual"}
   ['∼ :multivector]
   (fn dual [{{• '•} :ops {I- 'I-} :specials :as ga} a]
     (• a [I-]))

   ^{:doc "Dual"}
   ['∼ Blade]
   (fn dual [{{• '•} :ops {I- 'I-} :specials :as ga} a]
     (• [a] [I-]))

   ^{:doc "Normalize"}
   ['⧄ :multivector]
   normalize
   })

; p.80 the inverse of the pseudoscalar is simply the reverse
(defn with-specials [{b :basis-by-grade :as ga}]
  (let [
         I (peek b)
         I- (edalb I)
         S (first b)
       ]
    (assoc ga :specials
      {'I I 'I- I- 'S S})))

(defn compare-blades [{ag :grade ab :bitmap :as a} {bg :grade bb :bitmap :as b}]
  (if (== ag bg)
    (compare ab bb)
    (compare ag bg)))

(defn with-help [{ops :ops :as a-ga}]
  (-> a-ga
    (assoc :help
       (into {} (map (juxt first (comp :doc meta)) (filter meta (keys ops)))))
    (assoc :examples
       (into {} (map (juxt first (comp :e.g. meta)) (filter meta (keys ops)))))))

(defn ga
  ([p q r]
    (ga "e" p q r))
  ([prefix p q r]
    (ga {:prefix prefix :p p :q q :r r}))
  ([prefix base p q r]
    (ga {:prefix prefix :base base :p p :q q :r r}))
  ([{:keys [prefix base p q r pm qm rm md] :or {base 0 pm 1 qm -1 rm 0}}]
    (let [bases-of (bases-of prefix base (+ p q r))
          m {
             :metric (or md (vec (concat (repeat p pm) (repeat q qm) (repeat r rm))))
             :basis bases-of
             :basis-by-bitmap (reduce (fn [r [n b]] (assoc r (:bitmap b) n)) (vec (repeat (count bases-of) 0)) bases-of)
             :basis-by-grade (vec (sort compare-blades (vals bases-of)))
             :basis-in-order (reduce (fn [r [n b]] (assoc r (:bitmap b) b)) (vec (repeat (count bases-of) 0)) bases-of)
             :ops (ga-ops)
             }
          ; note ops must be in order of dependence because of the partial later
           ops '[+ * • ∼ 𝑒 ⍣ - ⁻ _ *' *0 •∧ ∧ ∨ ⍟ ⧄ op]
         ]
      (reduce
        (fn [r op]
          (update-in r [:ops op]
            (fn [g] (partial g r))))
        (with-help
          (reduce
            (fn [r op] (assoc-in r [:ops op] (compare-G op)))
            (with-specials m) ops)) ops))))

(defn multivector [{basis :basis} elements]
  (mapv
    (fn [[co bb]]
      (G (get basis bb bb) co))
    (partition 2 elements)))

(defn multivector? [f]
  (and
    (or
       (and (list? f) (every? (fn [[c b]] (and (number? c) (and (symbol? b) (nil? (namespace b))))) (partition 2 f)))
       (and (vector? f) (every? (fn [[c b]] (and (or (number? c) (list? c)) (symbol? b))) (partition 2 f))))
    (== 0 (mod (count f) 2))))

(defn bladelike [x]
  (if (symbol? x)
    (let [basis (namespace x) n (name x)]
      (try {:basis (symbol basis) :scale (Double/parseDouble n)}
        (catch Exception e nil)))
    nil))

(defn multivector1? [f]
  (and (list? f)
    (every? (fn [x] (or (bladelike x) (number? x))) f)))

(defmacro in-ga
  ([prefix p q r body]
    `(in-ga {:prefix ~prefix :base 0 :p ~p :q ~q :r ~r} ~body))
  ([prefix base p q r body]
    `(in-ga {:prefix ~prefix :base ~base :p ~p :q ~q :r ~r} ~body))
  ([p q r body]
    `(in-ga {:prefix e :base 0 :p ~p :q ~q :r ~r} ~body))
  ([{:keys [prefix base p q r metric] :or {base 0 prefix 'e}} body]
    (let [
          prefix (name prefix)
          {:keys [ops specials] :as ga} (ga 'e base 0 0 0)
           opz (into #{} (filter symbol? (keys ops)))
           o (complement (into #{} (filter symbol? (keys ops))))
           s (filter symbol? (tree-seq seqable? seq body))
           e (vec (cons (symbol (str prefix "_")) (filter (fn [i] (or (starts-with? (name i) prefix)
                                                                    (and (namespace i) (starts-with? (namespace i) prefix)))) s)))
           e (mapv (fn [x] (if-let [{b :basis} (bladelike x)] b x)) e)
           ]
       `(let [{{:syms ~e :as ~'basis} :basis
               {:syms ~(vec (filter symbol? (keys ops)))} :ops
               {:syms ~(vec (keys specials))} :specials :as ~'ga} (ga {:prefix ~prefix :base ~base :p ~p :q ~q :r ~r :metric ~metric})]
          ~(w/postwalk
             (fn [f]
               (cond
                 (and (multivector? f) (o (first f)))
                   (mapv (fn [[c b]] `(G ~b ~c)) (partition 2 f))
                 (and (multivector1? f) (o (first f)))
                   (mapv (fn [x] (let [{b :basis c :scale}
                                       (if (number? x)
                                         `{:basis ~(symbol (str prefix "_")) :scale ~x}
                                         (bladelike x))]
                                   `(G ~b ~c))) f)
                 :else f))
               body)))))

(defmethod print-method Blade [{:keys [basis scale]} writer]
  (doto writer
    (.write (str scale))
    (.write " ")
    (.write (str basis))))