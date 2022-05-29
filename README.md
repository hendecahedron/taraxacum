# abl-ajr

This library is an experimental work in progress from my notes taken while learning Geometric Algebra. At this stage it has the core GA functions in Euclidean metric:
```
  * geometric product
  • interior product
  ∧ exterior product
  ∼ dual
  ∨ regressive product    
  𝑒 exponential
```   


## The Geometric Numbers

`e₀` `e₁` `e₂`...`eᵢ`

These are basis elements of a GA and most people tend to use`e`for "element" as the prefix although some prefer something else like`v`

Basis elements are of type`Blade`which represents oriented volumes of space

`Blade{:bitmap 1, :scale 1.0, :grade 1, :basis e0}` is a basis blade of grade 1,or 1-blade

`Blade{:bitmap 5, :scale 1.0, :grade 2, :basis e02}` is a bivector of grade 2, or 2-blade 

`e_` (or `v_` or `whatever_`) is the scalar basis element of grade 0

the grade is the number of independent bases comprising the blade, the bitmap is the xor of those bases' bitmaps and the scale is its weight and orientation

create a blade directly`(Blade. bitmap scale grade)`or with`G`(for Geometric number)

`(require 'abl.ajr.core)`

`(G blade-or-bitmap scale)`

Multivectors are linear combinations of blades and are the general objects used in GA computations. Multivectors are represented with (Clojure) vectors of blades or vectors of coefficients and blades

`[0.3 e0 0.7 e1 3e-7 e2]` is a multivector


Create an algebra `(let [ga3d (ga 'e 3 0 0))` for Euclidean 3d space, `(ga prefix p q r)` where `p q r` is the metric signature: `p` number dimensions squaring to +1, `q` number dimensions squaring to -1 and `r` number dimensions squaring to 0.

This algebra is a map containing the basis elements and operations for computing with multivectors

`(:basis ga3d)` is a map of basis elements

`(let [{{:syms [e_ e0 e1 e2 e12 e01 e02]} :basis} (ga 'e 3 0 0))`

`(:ops ga3d)` is a map of operations for that GA

`(let [{{:syms [* • ∧]} :ops} ga3d])` some operations

the macro `in-ga` evaluates GA expressions in the context of a given GA

```
(in-ga 4 0 0
  (* [0.1 e0 0.2 e1 0.3 e2 -1.5 e3]
     [0.7 e0 0.4 e1 0.1 e2 9e-7 e3]))
```
                 
returns the geometric product `*` of two multivectors in 4d space with a default prefix`e`for the basis elements. `(in-ga prefix p q r body)` evaluates the body in the context of a GA of the specified metric.

The geometric product is anticommutative

```
(in-ga 2 0 0 (* e0 e1))
=> #abl.ajr.Blade{:bitmap 3, :scale 1.0, :grade 2, :basis e01}
(in-ga 2 0 0 (* e1 e0))
=> #abl.ajr.Blade{:bitmap 3, :scale -1.0, :grade 2, :basis e01}
```   

this is written e₀e₁ meaning the geometric product of e₀ and e₁ and the result is written e₀₁ which is a bivector spanning e₀ and e₁ and since e₀₁ = -e₁₀ we'll write basis blades eᵢⱼ `(i < j)`



The geometric product is the sum of the interior and exterior products

`(+ (• a b) (∧ a b))`

```
(in-ga 3 0 1 
    (let [a [5 e12 3 e0] 
          b [2 e3 4 e0]] 
      (+ (• a b) (∧ a b))))
=> [12.0 e_ 20.0 e012 6.0 e03 10.0 e123]

```
or we could use `•∧` which returns the parts of the geometric product separately

```
(in-ga 3 0 1 (•∧ [5 e12 3 e0] [2 e3 4 e0]))
=> {:• [12.0 e_], :∧ [20.0 e012 6.0 e03 10.0 e123]}

```

the interior product is of lower grade

```
(in-ga 3 0 1 (• [5 e12 3 e0] [2 e3 4 e0]))
=> [12.0 e_]

(in-ga 3 0 0 (• [5 e0] [2 e01]))
=> [10.0 e1]
```

the exterior product is of higher grade

```
(in-ga 3 0 1 (∧ [5 e12 3 e0] [2 e3 4 e0]))
=> [20.0 e012 6.0 e03 10.0 e123]
```


(give more examples)


### Notes on terminology

(notes here)

#### GA Input Source / Keyboard Layout

This library uses the symbols `• ∼ 𝑒 ⍣ ∧ ∨` so a GA-specific input source is needed for macOS -- I used Ukele to create one and mapped the fn key to switch sources. For Linux, there'll be a way to do it.


#### References:


Geometric Algebra for Computer Science, Leo Dorst, Daniel Fontijne, Stephen Mann
 — main reference for this code, along with https://geometricalgebra.org, [ganja.js](https://github.com/enkimute/ganja.js) and the Java reference implementation

A Guided Tour to the Plane-Based Geometric Algebra PGA, Leo Dorst

The Meet and Join, Constraints Between Them, and Their Transformation under Outermorphism, A supplementary discussion by Greg Grunberg of Sections 5.1—5.7 of the textbook Geometric Algebra for Computer Science, by Dorst, Fontijne, & Mann

Universal Geometric Algebra, David Hestenes

Geometric Algebra Primer, Jaap Suter

and the bivector discord, particularly Enki


#### Acknowledgements

Thanks to Georgiana and Alex for listening to my ideas. Thanks to Dr Geoff James who many years ago gave an inspiring presentation on Clifford Algebra.


#### License

Copyright © 2021 Matthew Chadwick

Distributed under the Eclipse Public License either version 1.0 or (at
your option) any later version.