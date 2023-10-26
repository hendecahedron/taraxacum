# Taraxacum

This library is an experimental work in progress from my notes taken while learning Geometric Algebra. At this stage it has the core GA functions in Euclidean metric:
```
  * geometric product
  • interior product
  ∧ exterior product
  ∼ dual
  ∨ regressive product    
  𝑒 exponential
```   


### Use

`hendecahedron/taraxacum {:git/url "https://github.com/hendecahedron/taraxacum.git" :sha "0b6fab772eb5577069c2b1eabfa9fa9e62f60066"` in your `deps.edn`

Basis elements are of type`Basis`which represents oriented volumes of space

`(require '[hendecahedron.taraxacum.core :refer :all])`

Taraxacum defaults to an unconventional 0-based indexing for bases, so null vectors start at p, but one can use the conventional 1-based indexing too.

Multivectors are linear combinations of bases and are the general objects used in GA computations. Multivectors are represented with (Clojure) vectors of bases or vectors of coefficients and bases

The macro `in-ga` evaluates GA expressions in the context of a given GA, e.g.

```
(in-ga 4 0 1
  (∧ [0.1 e0 0.2 e1 0.3 e2 -1.5 e3 1 e4]
     [0.7 e0 0.4 e1 0.1 e2 9e-7 e3 1 e4]))
```
                 
returns the intersection of two planes. `(in-ga p q r body)` evaluates the body in the context of a GA of the specified metric.


This multivector representation only works inside the macro `in-ga`, otherwise use `(G basis scale)` to make a vector of bases.


One can also create an algebra with `ga` e.g. `(let [ga3d (ga 3 0 1)])` for 3d PGA, or `(ga prefix p q r)` where `p q r` is the metric signature: `p` number dimensions squaring to +1, `q` number dimensions squaring to -1 and `r` number dimensions squaring to 0.

This algebra is a map containing the basis elements and operations for computing with multivectors

`(:basis ga3d)` is a map of basis elements

`(let [{{:syms [e_ e0 e1 e2 e12 e01 e02]} :basis} (ga 3 0 1))`

`(:ops ga3d)` is a map of operations for that GA

`(let [{{:syms [* + - • ∨ ∧ ⁻ ∼]} :ops} ga3d])` the main GA operations





(give more examples)


#### Status

Experimental and definitely going to change, however I'm using the current version for my projects and it works fine.

#### Notes

I changed the name of this library from abl-ajr to Taraxacum which is also a corruption of the original Arabic, and dandelion is a corruption of the original French.

My first implementation of this used multimethods, then I rewrote it using an experimental dispatch system and I'll probably rewrite it again.

#### GA Input Source / Keyboard Layout

This library uses the symbols `• ∼ 𝑒 ⍣ ∧ ∨` so a GA-specific input source is needed for macOS -- I used Ukele to create one and mapped the a key to switch sources. For Linux, there'll be a way to do it.


#### References:


Geometric Algebra for Computer Science, Leo Dorst, Daniel Fontijne, Stephen Mann
 — main reference for this code, along with https://geometricalgebra.org, [ganja.js](https://github.com/enkimute/ganja.js) and the Java reference implementation

A Guided Tour to the Plane-Based Geometric Algebra PGA, Leo Dorst

The Meet and Join, Constraints Between Them, and Their Transformation under Outermorphism, A supplementary discussion by Greg Grunberg of Sections 5.1—5.7 of the textbook Geometric Algebra for Computer Science, by Dorst, Fontijne, & Mann

Universal Geometric Algebra, David Hestenes

Geometric Algebra Primer, Jaap Suter

and the bivector discord


#### Acknowledgements

Thanks to Georgiana and Alex for listening to my ideas. Thanks to Dr Geoff James who many years ago gave an inspiring presentation on Clifford Algebra. And thanks to all in the bivector discord


#### License

Copyright © 2021 Matthew Chadwick

Distributed under the Eclipse Public License either version 1.0 or (at
your option) any later version.