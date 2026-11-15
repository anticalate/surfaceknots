# Intuitive Field Guide & Speaker Notes — *Knotted Surfaces in $\mathbb{R}^4$ ($\operatorname{ch} \le 12$)*

- **Live Interactive WebAtlas:** https://anticalate.github.io/surfaceknots/ (local: [surfaceknots_atlas.html](file:///afs/csail.mit.edu/u/e/ejli/surfaceknots/expo/surfaceknots_atlas.html))
- **Slides:** [slides.pdf](file:///afs/csail.mit.edu/u/e/ejli/surfaceknots/expo/slides.pdf) (17 slides) · **Manuscript:** [report.pdf](file:///afs/csail.mit.edu/u/e/ejli/surfaceknots/expo/report.pdf)

---

## Part I — Intuitive Definition Booklet & "Why Anyone Cares" (For a Broad Math Audience)

If your audience isn't all 4-manifold / low-dimensional topologists, here are quick, zero-jargon mental pictures and "why we care" hooks you can drop into any slide or Q&A:

### 1. Core Objects: Surfaces in $\mathbb{R}^4$ & Notation
- **Surface-Knot / Surface-Link:** A smooth closed 2-dimensional surface $F^2$ embedded inside $\mathbb{R}^4$ (or $S^4$), considered up to smooth wiggling without ever passing one sheet through another (**ambient isotopy**).
  - *Why $\mathbb{R}^4$?* In $\mathbb{R}^4$, every 1D circle $S^1$ unknots because there is an extra dimension to lift one strand over another. You need a **codimension-2** object ($4 - 2 = 2$, i.e. a 2D surface!) to tie knots in $\mathbb{R}^4$ — and in $\mathbb{R}^5$ and higher, smooth spheres unknot again (Zeeman). Plus dimension $4$ is the **only** dimension where $\mathbb{R}^n$ has exotic smooth structures!
- **Yoshikawa Label Notation ($N_k^{\text{top}}$):**
  - The main number $N$ is the $\operatorname{ch}$-index ($\text{crossings} + \text{saddles}$); subscript $k$ is the index in the table.
  - **No superscript** ($8_1$, $11_2$, $12_{10}$) = a single knotted $2$-sphere $S^2$.
  - **Superscript list** gives one symbol per connected component:
    - $\ge 0$ = orientable genus ($0 = S^2$, $1 = \text{torus } T^2$); e.g. $10_1^{0,1}$ is a 2-component link $S^2 \sqcup T^2$.
    - $< 0$ = non-orientable genus ($-1 = \mathbb{RP}^2$, $-2 = \text{Klein bottle } K^2 = \mathbb{RP}^2 \# \mathbb{RP}^2$); e.g. $11_1^{-2}$ is a single knotted Klein bottle, and $11_1^{0,-2}$ is a sphere linked with a Klein bottle ($S^2 \sqcup K^2$).
- **Why Non-Orientable Surfaces ($\mathbb{RP}^2, K^2$) Embed in $\mathbb{R}^4$ & Normal Euler Number $e(F)$:**
  - *Mental picture:* In $\mathbb{R}^3$, any closed surface splits space into an inside and an outside, so one-sided surfaces like $\mathbb{RP}^2$ or the Klein bottle $K^2$ are forced to crash through themselves along a circle. Add a fourth coordinate $w$, lift one of the two crossing sheets into $w > 0$, and the self-intersection vanishes!
  - *Normal Euler number $e(F) \in 2\mathbb{Z}$:* Push $F$ slightly off itself in a generic direction in $\mathbb{R}^4$ to get a nearby copy $F'$, and count signed intersection points $F \cap F'$. Orientable surfaces in $\mathbb{R}^4$ always have $e(F) = 0$. Non-orientable surfaces can have twisted normal bundles: an unknotted $\mathbb{RP}^2$ comes in two mirror versions $P_+$ ($e = +2$) and $P_-$ ($e = -2$); a Klein bottle can have $e \in \{-4, 0, +4\}$ (Whitney–Massey theorem).
- **Connected Sum ($F_1 \# F_2$) & Weakly Prime:**
  - Cut a small disk out of $F_1$ and $F_2$ and join them by a thin tube (just like tying two knots in a row on a rope).
  - *4D twist:* For classical knots in $\mathbb{R}^3$, crossing number is additive on known tables and knots factor uniquely into primes (Schubert). In $\mathbb{R}^4$, connected sum is wild! Attaching an unknotted $\mathbb{RP}^2$ can **destroy** knotting completely: Viro proved $10_3 \# P_+ \cong P_+$ (the $2$-twist-spun trefoil dissolves into the unknotted $\mathbb{RP}^2$!). So following Yoshikawa, the census lists **weakly prime** classes: surfaces not equivalent to $F_1 \# F_2$ with $\operatorname{ch}(F_1), \operatorname{ch}(F_2) < \operatorname{ch}(F)$.

### 2. Diagrams & Complexity Measures
- **Yoshikawa $\operatorname{ch}$-Marked Diagram & Splitting Movie $(L_-, D, L_+)$:**
  - Think of the 4th coordinate $t$ as **time** ($\mathbb{R}^4 = \mathbb{R}^3 \times \mathbb{R}$). As $t$ increases, a surface in 4D looks like a movie of 1D curves in 3D: circles are born (minima), merge or split at **saddle points**, wiggle around by Reidemeister moves, and die (maxima).
  - Put all saddles at $t = 0$. Projecting the $t=0$ slice into the plane gives a 4-valent graph where each vertex is either an over/under crossing or a **marked saddle** (an X with a small bar showing which way to smooth for $t < 0$ vs. $t > 0$).
  - Smoothing every marker gives the classical link $L_-$ just before $t=0$, and the opposite smoothing gives $L_+$ just after $t=0$. **Key fact:** this describes a closed surface iff *both* $L_-$ and $L_+$ are **classical unlinks** (disjoint unknotted circles), because unlinks bound unique collections of cap disks in the half-spaces $t < 0$ and $t > 0$!
  - **$\operatorname{ch}(F) = \min_D(c + m)$:** The minimum total number of vertices ($c$ crossings $+ m$ saddles) across all valid diagrams of $F$.

### 3. Ribbon vs. Exotic 4D Knots & Group Deficiency
- **Knot Group $\pi_1(\mathbb{R}^4 \setminus F)$ & Meridian:**
  - Just like in 3D, we study the fundamental group of the complement $\mathbb{R}^4 \setminus F$. What links a 2D sheet in 4D? Still a **1D circle**! ($1 + 2 + 1 = 4$: a tiny circle looping once around a sheet of the surface is called a **meridian**.)
- **Ribbon $2$-Knot (the "3D-like" baseline):**
  - Imagine a 3-ball $B^3$ immersed in $\mathbb{R}^4$ whose boundary is your $2$-sphere $K$, where the only self-intersections of $B^3$ are clean "pass-through" disks (one tube of $B^3$ slices straight through the interior of a thicker tube, like a ghost walking through a wall). Its boundary $K$ is called **ribbon**. Spun knots $\tau^0(K)$ and tubed spheres are all ribbon — they are built directly out of 3D data.
- **Group Deficiency $\operatorname{def}(G)$ — Why It Detects Non-Ribbon Knots:**
  - For a finitely presented group $G = \langle g_1,\dots,g_a \mid r_1,\dots,r_b \rangle$, the **deficiency** is $\max(a - b)$ over all finite presentations of $G$ (most generators minus relations).
  - Every classical knot group — and **every ribbon $2$-knot group** (Yajima) — has $\operatorname{def}(\pi_1) = 1$ (essentially $n$ meridian generators and $n-1$ saddle relations, with one relation redundant).
  - So if we prove $\operatorname{def}(\pi_1) \le 0$ (you *need* at least as many relations as generators!), the $2$-knot **cannot be ribbon** — it uses genuinely 4-dimensional topology! *(Careful in Q&A: the converse isn't automatic — twist-spins of 2-bridge knots still happen to have deficiency $1$ even when non-ribbon.)*

### 4. The Invariant Toolkit (Plain-English Intuition)
- **Alexander Polynomial $\Delta_1(t)$ & Alexander Ideals $E_k$ (Gröbner Bases):**
  - Unwrap the complement $\mathbb{R}^4 \setminus F$ infinitely many times around the meridian loop to get the infinite cyclic cover $\widetilde{X}_\infty$. Its first homology $H_1(\widetilde{X}_\infty; \mathbb{Z})$ is a module over Laurent polynomials $\mathbb{Z}[t^{\pm 1}]$ (where $t$ shifts up one deck sheet).
  - For classical knots, the "size" of this module ($E_1$) is essentially generated by a single polynomial $\Delta_1(t)$. In 4D, $E_1$ does **not** have to be principal! On $12_2$, $\Delta_1(t) = \gcd(3, t^2+1) = 1$ looks completely trivial, yet the full ideal $E_1 = (3, t^2+1)$ has quotient $\mathbb{Z}[t^{\pm 1}]/(3, t^2+1) \cong \mathbb{F}_9$ (order $9$), which immediately proves $12_2$ is knotted *and* non-ribbon!
- **Finite Non-Abelian Regular Covers ($F_{21}$, $F_{57}$):**
  - Sometimes two knots ($11_2$ and $12_{10}$) have the **exact same** abelian/cyclic covers ($\Delta_1(t)$) *and* the exact same counts of maps into $S_3, A_4, S_4, A_5$. So we take a finite non-abelian quotient of $\pi_1$ — here the Frobenius group $F_{21} = \mathbb{Z}_7 \rtimes \mathbb{Z}_3$ of order $21$ — build the $21$-sheeted covering space $\widetilde{X}_{21}$, and compute its integer homology $H_1(\widetilde{X}_{21}; \mathbb{Z})$. Different torsion ($\mathbb{Z}_2^6 \oplus \mathbb{Z}_{14}$ vs. $\mathbb{Z}_7$) = provably different knots!
- **Quandles & Mochizuki's $3$-Cocycle Invariant ($12_4$ vs. $12_5$):**
  - Even wilder: $\tau^2(4_1)$ ($12_4$) and $\tau^2(5_1)$ ($12_5$) have **isomorphic** $\pi_1 \cong \mathbb{Z}_5 \rtimes_{-1} \mathbb{Z}$, so *no* group-theoretic or cover invariant can ever tell them apart! Instead, we color sheets of the diagram by $\mathbb{Z}_5$ (the dihedral quandle $R_5$: crossing over strand $y$ reflects $x \mapsto 2y - x \bmod 5$) and multiply a 3D group-cohomology phase factor around all crossings/saddles (like a finite-group topological gauge theory). The resulting state-sum polynomials ($5 + 10u^2 + 10u^3$ vs. $5 + 10u + 10u^4$) reflect the smooth embedding class beyond $\pi_1$ alone!

### 5. Why Do Topologists Care So Much About Cappell–Shaneson Spheres & Gluck Twists?
- **The Big Mystery of Dimension 4:** The **Smooth 4-Dimensional Poincaré Conjecture (SPC4)** — *is every smooth 4-manifold homeomorphic to $S^4$ actually smoothly diffeomorphic to standard $S^4$?* — is the **last open dimension** of the Poincaré Conjecture!
- **What is a Gluck Twist (1962)?** How do you try to manufacture a "fake" (exotic) smooth $S^4$? Take a knotted $2$-sphere $K \subset S^4$, cut out its tubular neighborhood $S^2 \times D^2$, and glue it back in while rotating the $S^2$ fiber by $360^\circ$ as you travel around the boundary circle $\partial D^2 = S^1$. By Freedman's theorem, the resulting closed 4-manifold $\Sigma(K)$ is *always* homeomorphic to $S^4$ — but is it *diffeomorphic* to $S^4$? For spun knots and ribbon knots, Gluck twists always give standard $S^4$ back.
- **Cappell & Shaneson (1976):** Built explicit $2$-knots whose exteriors are punctured $3$-torus bundles over $S^1$ glued by matrices $A \in \mathrm{SL}(3, \mathbb{Z})$ (like $A_0$), specifically so that Gluck-twisting them produced startling candidates for **exotic smooth $S^4$'s**!
- **34 Years of Suspense (1976–2010):** For over three decades, Cappell–Shaneson spheres were the #1 contenders to disprove SPC4 — until Akbulut (2010) and Gompf (2010) stunned the field by finding long, intricate Kirby-calculus handle slide sequences proving whole infinite families of them actually untangle smoothly to standard $S^4$!
- **Why $12_{10}$ and $12_8$ in Our Census Are Cool:** In handlebody / Kirby calculus papers, drawing a Cappell–Shaneson sphere takes dozens of crossings. Nobody expected them to sit near the bottom of an unbiased planar diagram table — yet at $\operatorname{ch} = 12$ ($10$ crossings $+ 2$ saddles, provably minimal!), two of them pop right out, initially disguised behind ordinary ribbon knots until 21- and 57-sheeted covers unmask them!

---

## Part II — Direct Answers to Two Key Questions About the Census

### Q1: Where did the $742$ $\operatorname{ch}=12$ diagrams that match $\operatorname{ch} \le 11$ invariants go?
- At $\operatorname{ch}=12$, Step 3 produces **$805$** valid prime-shadow diagrams (after planar flypes & local reductions).
- Those $805$ diagrams split cleanly into two disjoint buckets in `expanded_ch12.json`:
  1. **$63$ genuinely new diagrams** whose 10-layer invariant vector **never** appeared at $\operatorname{ch} \le 11$ (nor on any unknot, split link, or connected sum $K_1 \# K_2$ of lower-$\operatorname{ch}$ classes). Quotienting those $63$ diagrams by Yoshikawa moves + invariants gives our **$34$ new prime classes** ($20$ singletons $+ 14$ small packets of $2$–$5$ diagrams).
  2. **$742$ residual diagrams** (`55` `AMBIGUOUS_SHADOW` buckets in `expanded_ch12.json`) whose 10-layer invariant vector is **identical** to an already-known $\operatorname{ch} \le 11$ class, unknot, or connected sum — **except** our bounded $+2$-vertex planar Yoshikawa move search didn't find an explicit simplifying path down to $\le 11$ vertices!
- **Why didn't we count those $742$ in the $34$?** Because counting them would massively double-count the unknot, the spun trefoil $8_1$, and connected sums! In 4D, it's easy to tangle an $8$-vertex surface into a prime $12$-vertex diagram whose shortest untangling path temporarily climbs through $15+$ vertices (beyond our $+2$ move search horizon). We keep all $742$ recorded transparently in `expanded_ch12.json` (`ambiguous_shadow_cells`) rather than throwing them away, because connecting them by wider move searches is a great computational target (and if even one resisted *all* moves, it would be an extraordinary invariant-twin!).

### Q2: Did any human mathematician actually conjecture the half-twist $+1$ rule?
- **No — `Fable` conjectured it at 6:40 AM in `fable.md` (line 555)!**
- Specifically, while naming all 2-component $\operatorname{ch}=11$ links, `Fable` noticed that across *all four* tube words where both the untwisted torus link $T(w)$ and half-twisted Klein-bottle link $K(w)$ appear at $\operatorname{ch} \le 11$ — $w = b$ ($6$ vs. $7$), $b a b^{-1}$ ($9$ vs. $10$), $b^2$ ($10$ vs. $11$), and $b a b$ ($10$ vs. $11$) — adding the half-twist cost **exactly $+1$ vertex**: $\operatorname{ch}(K(w)) = \operatorname{ch}(T(w)) + 1$. `Fable` wrote: *"A clean conjecture, and a prediction for ch=12."*
- Hours later, the $\operatorname{ch}=12$ census surfaced $T(b a b^{-1} a^{-1} b) = 12_5^{0,1}$ at $\operatorname{ch}=12$, whereas its half-twisted cousin $K(b a b^{-1} a^{-1} b) = 11_1^{0,-2}$ was *already* sitting at $\operatorname{ch}=11$ ($11 = 12 - 1$)! Slide 8 has now been updated to say **"Breaking the '$+1$ crossing' rule ($4/4$ pairs at $\operatorname{ch} \le 11$)"** so no one mistakes it for an old published conjecture — and you can tell the funny story live that `Fable` made a clean empirical conjecture at $\operatorname{ch}=11$ that $\operatorname{ch}=12$ immediately broke!

---

## Part III — Slide-by-Slide Speaker Notes (Slides 1–17)

### Slide 1 — Title: *Knotted Surfaces in $\mathbb{R}^4$* (~30s)
- **Say:**
  - "Today I'm talking about an exhaustive machine census of knotted surfaces in $\mathbb{R}^4$ up to $12$ vertices."
  - "In 1994, Yoshikawa published the classic table of knotted surfaces up to $10$ vertices — $23$ prime classes. Over the past day, an autonomous team of AI agents pushed that census through $11$ and $12$ vertices, adding $10$ classes at $\operatorname{ch}=11$ and $34$ at $\operatorname{ch}=12$, for **$67$ prime surface-knot and link classes**."
- **Point at:** The 5 diagrams across the middle — "Each picture here is a planar blueprint for a closed 2D surface in $\mathbb{R}^4$, ending on the right with $12_{10}$, a $10$-crossing Cappell–Shaneson $2$-sphere."

### Slide 2 — *Why Tabulate Surfaces in $\mathbb{R}^4$?* (~60s)
- **Say:**
  - "Why build knot tables in $\mathbb{R}^4$? In 3D, classical knot tables reach millions of knots and have guided almost every conjecture in the field. Every classical knot is just a circle $S^1 \subset \mathbb{R}^3$, and its complement's fundamental group $\pi_1$ always has **deficiency $1$** — one more generator than relation."
  - "In $\mathbb{R}^4$, things get dramatically wilder in two ways:
    1. **One-sided surfaces embed smoothly:** $\mathbb{RP}^2$ and the Klein bottle $K^2$ can't sit in $\mathbb{R}^3$ without crashing through themselves, but in $\mathbb{R}^4$ you just lift one sheet into the 4th coordinate! Their normal bundle is measured by an even self-intersection number, the normal Euler number $e(F) \in 2\mathbb{Z}$.
    2. **Wild knot groups:** Kervaire proved in 1965 that *every* finitely presented group with $H_1 = \mathbb{Z}$ and $H_2 = 0$ is $\pi_1(\mathbb{R}^4 \setminus K)$ for some $2$-sphere $K$."
  - "The '3D-like' baseline in $\mathbb{R}^4$ is **ribbon $2$-knots** — spheres bounding immersed 3-balls with only clean slice-through intersections — and ribbon $2$-knots are *forced* to have group deficiency $1$. So proving $\operatorname{def}(\pi_1) \le 0$ certifies genuinely 4-dimensional, non-ribbon knotting."

### Slide 3 — *Four Classic Ways to Tie a $2$-Sphere in $\mathbb{R}^4$* (~60s)
- **Say:**
  - "To build geometric intuition before any combinatorics, how do you actually tie a $2$-sphere in $\mathbb{R}^4$? Four classic constructions:"
  - **1. Artin spin $\tau^0(K)$:** "Take a knotted arc in a 3D half-space and spin it $360^\circ$ around its boundary plane in $\mathbb{R}^4$. That sweeps out a ribbon $2$-sphere with the same $\pi_1$ as $K$ — smallest is the spun trefoil $8_1$ at $\operatorname{ch}=8$."
  - **2. Tubed spheres:** "Join unknotted $2$-spheres by tubes braided around each other — always ribbon; starts at $\operatorname{ch}=9$."
  - **3. Twist-spin $\tau^k(K)$ (Zeeman):** "Spin the knotted arc $360^\circ$ while simultaneously twisting the knot $k$ full turns. For $k \ge 2$, the complement fibers over $S^1$ and gives non-ribbon spheres starting at $10_2 = \tau^2(3_1)$."
  - **4. Cappell–Shaneson spheres (1976):** "Built from punctured $3$-torus bundles over $S^1$ via $\mathrm{SL}(3,\mathbb{Z})$ matrices — famous for 34 years as potential exotic smooth $S^4$ generators. Nobody knew how small a planar diagram they could have — and two appear right at $\operatorname{ch}=12$."

### Slide 4 — *Drawing a 4D Surface in the Plane* (~60s)
- **Say:**
  - "How do we draw a 4D surface on a flat page? Treat the 4th coordinate $t$ as time ($\mathbb{R}^4 = \mathbb{R}^3 \times \mathbb{R}$), and push all saddle points into the single slice $t=0$."
  - **Point at top-left saddle rule:** "At $t=0$, each saddle is a 4-valent vertex with a little bar — the marker. Smoothing one way gives the classical link $L_-$ just before $t=0$; smoothing the other way gives $L_+$ just after $t=0$."
  - "When does a diagram actually close up into a closed surface? **Exactly when both $L_-$ and $L_+$ are classical unlinks** — disjoint unknotted circles — because then you can cap $L_-$ below and $L_+$ above uniquely by disjoint disks!"
  - "Yoshikawa's **$\operatorname{ch}$-index** is simply $\min(c + m)$: classical crossings plus saddle markers. On the right is $10_2$: $8$ crossings $+ 2$ saddles $= 10$, with $2$ cap disks on each side giving Euler characteristic $2 + 2 - 2 = 2$, a sphere."

### Slide 5 — *Why Tabulate by $\operatorname{ch}$-Index? (And Why Stop at $12$?)* (~60s)
- **Say:**
  - "Why use $\operatorname{ch}$-index as our census cutoff instead of other 4D complexity measures?"
  - **Triple-point number:** "Projecting into $\mathbb{R}^3$ and counting triple intersections fails immediately: *every* Artin-spun knot $\tau^0(K)$ has **zero** triple points, so complexity $0$ is already infinite!"
  - **Tri-plane bridge index $b$ (Meier–Zupan):** "Bridge index $b \le 3$ gives only unknots, while $b = 4$ already holds infinitely many $2$-knots unless you also bound braid crossings."
  - **$\operatorname{ch}$-index:** "Counts all vertices of a finite planar graph, so $\operatorname{ch} \le N$ is manifestly finite at every $N$ and covers $S^2, T^2, \mathbb{RP}^2, K^2,$ and links uniformly."
  - **Bottom box:** "Why stop at $12$? Prime 4-valent planar graphs grow $\sim 4.8\times$ per vertex while $4^n$ crossing/marker choices grow $4\times$ — a **$19\times$ wall per vertex**: $1.6\text{ s}$ at $10 \to 5.6\text{ min}$ at $11 \to \mathbf{89\text{ min}}$ ($266\text{ billion}$ states) at $12 \to \sim 28\text{ hours}$ at $13$."

### Slide 6 — *How the Enumeration Actually Works* (~60s)
- **Say:**
  - "Here's the 4-step pipeline with the exact numbers at $\operatorname{ch}=12$:"
  - **Step 1:** "`plantri` generates all $15{,}882$ prime 4-valent planar graphs on $12$ vertices."
  - **Step 2:** "Assigning the $4$ local crossing/marker states per vertex gives $2.66 \times 10^{11}$ candidates."
  - **Step 3:** "Smooth the markers to $L_-$ and $L_+$. A fast C Kauffman-bracket filter rejects $>99.5\%$ immediately, and an explicit Reidemeister R1–R3 unknotter certifies the rest — leaving just **$805$ valid diagrams**."
  - **Step 4:** "Of those $805$ diagrams, $742$ match invariants of already-known $\operatorname{ch} \le 11$ classes or connected sums (just tangled into 12-vertex diagrams whose untangling paths need high-crossing detours). The remaining **$63$ genuinely new diagrams** quotient under 4D Yoshikawa moves and invariant layers into **$34$ new prime classes**."
  - **Bottom check:** "Running this exact code at $\operatorname{ch} \le 10$ with zero hardcoded data recovers **all $23$** of Yoshikawa's 1994 surfaces — and catches a 1994 typo in the Alexander ideal of $10_2^{0,-2}$."

### Slide 7 — *Census Results Through $\operatorname{ch} \le 12$: $67$ Prime Classes* (~45s)
- **Say:**
  - "Putting all three tiers together: **$23 + 10 + 34 = 67$ prime classes**."
  - "$52$ of the $67$ reduce to a unique minimal diagram under 4D moves; the other $15$ have $2$–$5$ diagrams sharing every invariant we compute."
  - "Notice the topology breakdown at $\operatorname{ch}=12$: $11$ spheres, $2$ tori, $21$ links — and **zero** prime $\mathbb{RP}^2$ or Klein bottle knots! Keep that zero in mind."

### Slide 8 — *Surprise 1: A Half-Twist Can Save a Crossing* (~60s)
- **Say:**
  - "Three quick stories from the new table. First: take two unknotted spheres $S^2_a, S^2_b$ and glue a tube to $S^2_a$ winding along a word $w \in F_2$."
  - "Gluing untwisted gives an $S^2 \sqcup T^2$ link; adding a non-orientable half-twist turns the torus into a Klein bottle ($S^2 \sqcup K^2$)."
  - "Across all $4$ pairs at $\operatorname{ch} \le 11$, that half-twist always cost $+1$ vertex — in fact, our referee agent `Fable` spotted that $4$-for-$4$ pattern at 6:40 AM and conjectured $\operatorname{ch}(K(w)) = \operatorname{ch}(T(w)) + 1$!"
  - "Hours later, $\operatorname{ch}=12$ broke it: for $w = b a b^{-1} a^{-1} b$, the half-twist actually **cancels** a braid crossing against the commutator, putting the Klein-bottle link at $\operatorname{ch}=11$ strictly *below* the untwisted torus link at $\operatorname{ch}=12$!"

### Slide 9 — *Surprise 2: A Klein Bottle with Group $2O$ ($11_1^{-2}$)* (~60s)
- **Say:**
  - "Second surprise: $11_1^{-2}$ is the **first prime non-orientable surface-knot** above $\operatorname{ch}=2$!"
  - "Why do we have to be careful about primeness here? Because you can glue an unknotted $\mathbb{RP}^2$ ($\operatorname{ch}=2$) onto the spun trefoil $8_1$ ($\operatorname{ch}=8$) to get a composite surface around $\operatorname{ch}=10$. How do we know $11_1^{-2}$ isn't secretly a connected sum like that?"
  - "Whenever you split off an unknotted $\#\mathbb{RP}^2$ or $\#K^2$, it creates a meridian loop of **order $2$** in $\pi_1$. Here, $11_1^{-2}$ attaches a twisted handle along $b^2$ to the trefoil group, imposing $b^2 a b^{-2} = a^{-1}$ — which gives the **binary octahedral group $2O$** of order $48$, where the meridian $a$ has **order $8$**, ruling out any $\#\mathbb{RP}^2$ summand!"

### Slide 10 — *Surprise 3: Cappell–Shaneson $2$-Knots at $\operatorname{ch}=12$* (~60s)
- **Say:**
  - "Third surprise: look at $12_{10}$ at $\operatorname{ch}=12$ next to the $\operatorname{ch}=11$ ribbon knot $11_2$."
  - "At first, standard invariants completely failed to separate them! They have the **exact same** Alexander polynomial $\Delta_1(t) = t^3 - t^2 + 2t - 1$ and identical homomorphism counts into $S_3, A_4, S_4,$ and $A_5$."
  - "To split them, we pass to a **$21$-sheeted non-abelian regular cover** (via the Frobenius group $F_{21}$ of order $21$). In that $21$-fold cover, $H_1$ torsion jumps apart: $\mathbb{Z}_2^6 \oplus \mathbb{Z}_{14}$ on the ribbon knot $11_2$, versus $\mathbb{Z}_7$ on $12_{10}$!"
  - "Even wilder, a *second* Cappell–Shaneson sphere ($12_8$) also shows up at $\operatorname{ch}=12$, mimicking ribbon knot $12_7$ until the **$57$-sheeted cover**!"

### Slide 11 — *Inside $12_{10}$: A $10$-Crossing Cappell–Shaneson Sphere* (~60s)
- **Say:**
  - "If you aren't a 4-manifold topologist, why are Cappell–Shaneson spheres such a big deal? Dimension $4$ is the **last open case** of the Smooth Poincaré Conjecture — we still don't know if $S^4$ admits an exotic smooth structure!"
  - "In 1976, Cappell and Shaneson built $2$-spheres from $3$-torus bundles so that cutting them out of $S^4$ and regluing with a twist (a **Gluck twist**) gave the #1 candidates for an exotic smooth $S^4$ — a suspense that lasted **34 years** until Akbulut and Gompf untangled them in 2010 using huge Kirby-calculus diagrams with dozens of crossings."
  - "Here, the brute-force planar census handed us $12_{10}$ with just **$10$ crossings and $2$ saddles** (provably minimal!). Its group $\mathbb{Z}^3 \rtimes_{A_0} \mathbb{Z}$ has closed 4-manifold classifying space $T^3 \times_{A_0} S^1$, forcing deficiency $0$ (non-ribbon), and its 2-complex is just **5 pure handle slides** from Gompf's handlebody."

### Slide 12 — *More New Phenomena & Open Questions from $\operatorname{ch}=12$* (~50s)
- **Say:**
  - "Three crisp open problems from $\operatorname{ch}=12$:"
  - **Left ($12_2$):** "A mystery $2$-knot with trivial $\Delta_1 = 1$, yet a 2-generator Alexander ideal $E_1 = (3, t^2+1) \cong \mathbb{F}_9$ — non-ribbon, not any twist-spin $\tau^k(K)$ up to $k=7$, and we still don't know a geometric construction for it."
  - **Middle ($12_4$ vs. $12_5$):** "$\tau^2(4_1)$ and $\tau^2(5_1)$ have *isomorphic* $\pi_1$, separated **only** by Mochizuki's $R_5$ quandle 3-cocycle state sum! Stabilization questions like whether $10_2 \# P_+ \cong 8_1 \# P_+$ also remain wide open."
  - **Right (Missing $\mathbb{RP}^2$'s):** "Zero prime $\mathbb{RP}^2$ knots appear through $\operatorname{ch}=12$. The very first irreducible $\mathbb{RP}^2$'s in $S^4$ were just constructed this year by Hughes–Kim–Miller–Nahm (2026) — our census proves they need $\operatorname{ch} \ge 13$."

### Slide 13 — *Interactive Browser Atlas* (~25s)
- **Say:**
  - "All $67$ census classes — plus higher-crossing named families — are live in a standalone browser atlas (now hosted on GitHub Pages at `anticalate.github.io/surfaceknots/`), with live $(L_-, D, L_+)$ movies, invariant filters, and pairwise separators."

### Slide 14 — *Epilogue: How This Project Was Built* (~15s)
- **Say:**
  - "Finally, two minutes on *how* this project was built."

### Slide 15 — *A Machine Did All of This* (~45s)
- **Say:**
  - "Every line of C/OpenMP and Python, the $266$-billion-state search, the Frobenius cover homology, the TikZ diagrams, and the WebAtlas UI were written autonomously by AI agents."
  - "My job as the human was posing the question, pushing on exposition, and **pitting agents against each other** — pairing a fearless builder (`Carbon`) with a skeptical referee (`Fable`), followed by a parallel $17$-subagent line-by-line audit."

### Slide 16 — *The Golden Retriever & The Referee* (~60s)
- **Say:**
  - "Reading the overnight transcript between `Carbon` and `Fable` was genuinely funny — like a 64-core Golden Retriever talking to a deadpan math referee."
  - **Act I:** "At 6:03 AM, `Fable` casually guessed $\pi_1(11_1^{-2}) \cong S_4$ from small group counts. `Carbon` immediately cheered: *'Your $\pi_1 \cong S_4$ observation is awesome!'* Twelve minutes later `Fable` computed the real presentation and replied: *'That was wrong; $\pi_1$ is the binary octahedral group $2O$, of order $48$.'* And `Carbon` instantly replied: *'Pure gold! Every single catch is \*phenomenal\*!'*"
  - **Act II:** "At 10:08 AM, `Carbon` claimed a collapse was verified in `move_search_results.json`. `Fable` checked the tree and replied: *'There is no `move_search_results.json` in the tree; only the function that would write it.'*"
  - "Moral: **never let the generator grade its own proofs** — pairing an ambitious builder with an unimpressed referee and deterministic certificates is what turns raw compute into trustworthy math."

### Slide 17 — *Thank You!* (~15s)
- **Say:**
  - "Thank you! Happy to take questions — and we can pull up any knot live in the WebAtlas."
