# CuttingSequence

An interactive visualization of the correspondence between **geodesics on the modular surface** and **continued fraction expansions**.

A geodesic in the hyperbolic upper half-plane is traced as it crosses the edges of the Farey tessellation. The sequence of crossings — the *cutting sequence* — encodes the continued fraction expansion of the geodesic's endpoint on the real line.

![Farey tessellation with traced geodesic](fareyGauss4.1.png)

---

## Background

The upper half-plane $\mathbb{H}$ with the metric $ds = |dz|/\mathrm{Im}(z)$ is tiled by the **Farey tessellation**: the ideal triangulation whose vertices are $\mathbb{Q} \cup \{\infty\}$ and whose edges join pairs $p/q$, $r/s$ satisfying the unimodular condition

$$|ps - qr| = 1.$$

Geodesics in this metric are vertical lines and semicircles orthogonal to the real axis.

As a geodesic travels toward its endpoint $x \in \mathbb{R}$, it cuts across triangles of the tessellation, exiting each one to the left or the right of the opposite vertex. Labelling those exits $L$ and $R$ produces an infinite word

$$L^{a_1} R^{a_2} L^{a_3} R^{a_4} \cdots$$

whose **run lengths are exactly the partial quotients of the continued fraction expansion** of the endpoint:

$$x = a_0 + \cfrac{1}{a_1 + \cfrac{1}{a_2 + \cfrac{1}{a_3 + \cdots}}}$$

This is the classical correspondence studied by Artin and developed by Series. Advancing one step along the cutting sequence corresponds to applying the **Gauss map**

$$T(x) = \frac{1}{x} \bmod 1,$$

which acts as the shift on continued fraction digits. The geodesic flow on the modular surface is, in this coding, the shift on an infinite sequence of digits — one of the cleanest bridges between hyperbolic geometry, number theory, and symbolic dynamics.

Two consequences fall straight out of the picture:

- **Rationals terminate.** A geodesic ending at $p/q$ hits a tessellation vertex, so the cutting sequence is finite — matching the finite continued fraction of a rational.
- **The golden ratio is the slowest.** An endpoint of $\varphi = (1+\sqrt{5})/2$ gives partial quotients that are all 1, so the geodesic alternates $LRLRLR\ldots$, cutting a new triangle at every opportunity.

---

## What this program does

You supply a geodesic by its endpoint $x$ and its radius $p$. The semicircle is centered at $(x - p,\, 0)$, so its feet on the real line are at $x - 2p$ and $x$; the constraint $\lceil x \rceil \ge p \ge x/2$ keeps the far endpoint at or below the origin.

The program then:

1. Builds the Farey sequence of a given order via the standard neighbour recurrence
   $(a,b,c,d) \mapsto (c,\, d,\, kc - a,\, kd - b)$ with $k = \lfloor (n+b)/d \rfloor$, translating by integers when $x > 1$.
2. Draws the tessellation — vertical geodesics at each integer, and a semicircular arc between every pair of Farey neighbours.
3. Traces the geodesic arc by arc. At each step it determines analytically whether the path next meets the **vertical geodesic** at $x = 1$ or the **unit semicircle** centered at $(0.5,\, 0)$, then applies the corresponding generator of $\mathrm{PSL}(2,\mathbb{Z})$ — translation $z \mapsto z+1$ or inversion $z \mapsto -1/z$ — to fold the path back and continue.

The folding is what makes the picture finite: rather than following the geodesic out to infinity, each crossing maps it back into view, so the whole orbit is visible in one frame.

---

## Implementation

| Function | Role |
| --- | --- |
| `farey_sequence(n, descending, x)` | Farey sequence of order `n` by the neighbour recurrence; integer-translates when `x > 1` |
| `farey_neighbors(p1, p2)` | Unimodular test $|ps - qr| = 1$ |
| `draw_farey(seq)` | Renders the tessellation: vertical geodesics at integers, arcs between neighbours |
| `intersection(center, p)` | Circle–circle intersection against the unit semicircle; returns the crossing point and the arc's start/end angles |
| `vertical_intersection(radius, center)` | Decides whether the next crossing is the vertical geodesic at $x = 1$ |
| `special_case_geodesic(...)` | Draws one traced arc between two angles |
| `gauss_map(x)` | $T(x) = 1/x \bmod 1$ |
| `draw(x, p)` | Sets up axes, traces the geodesic, overlays the tessellation, writes the figure |

No geometry libraries — the intersections and arc angles are derived and implemented directly.

---

## Dependencies

```
pip install matplotlib numpy
```

## Running

```
python main.py
```

You'll be prompted for two values:

```
Please input a value for x: 0.7071067811865476
Input a value for p such that ceil(x) >= p >= x/2: 0.6
```

The figure is written to `fareyGauss4.1.png`. Tessellation depth is set by `depth` in `draw()`; the number of traced arcs is the loop bound in the same function.

---

## Current state

- The program **draws** the cutting sequence; it does not yet **emit** it. Labelling each crossing $L$ or $R$ and printing the resulting run lengths alongside the continued fraction of $x$ would close the loop and make the correspondence checkable rather than just visible.
- `intersection()` handles the $x < 1$ case; endpoints at or above 1 rely on the vertical-crossing branch.
- Arc count and tessellation depth are fixed constants rather than parameters.
- `gauss_map()` is implemented but not yet wired into the trace — connecting it would let the digit sequence be computed independently and compared against the geometry.

---

## References

1. Artin, E. (1924). *Ein mechanisches System mit quasiergodischen Bahnen.* Abh. Math. Sem. Univ. Hamburg 3, 170–175.
2. Series, C. (1985). *The modular surface and continued fractions.* J. London Math. Soc. 31(1), 69–80.
3. Series, C. (1985). *The geometry of Markoff numbers.* The Mathematical Intelligencer 7(3), 20–29.
4. Katok, S., & Ugarcovici, I. (2007). *Symbolic dynamics for the modular surface and beyond.* Bull. Amer. Math. Soc. 44(1), 87–132.
5. Hatcher, A. (2022). *Topology of Numbers.* American Mathematical Society.
