# Speaker Notes — *Knotted Surfaces in $\mathbb{R}^4$: An Exhaustive Machine Census Through $\operatorname{ch} \le 12$*

**Deck:** [slides.pdf](file:///afs/csail.mit.edu/u/e/ejli/surfaceknots/expo/slides.pdf) (17 slides, target pace ~12–15 minutes + Q&A)  
**Interactive Demo (keep open in a browser tab for Q&A):** [surfaceknots_atlas.html](file:///afs/csail.mit.edu/u/e/ejli/surfaceknots/expo/surfaceknots_atlas.html)

---

## Slide 1 — Title: *Knotted Surfaces in $\mathbb{R}^4$* (~30s)
- **What to say:**
  - "Today I'm presenting an exhaustive computer census of knotted surfaces in $\mathbb{R}^4$ up to complexity $12$."
  - "Back in 1994, Katsuyuki Yoshikawa published the classic table of knotted surfaces up to $10$ vertices — $23$ prime surfaces total. Over the past day, an autonomous team of AI agents pushed that census through $11$ and $12$ vertices, adding $10$ classes at $\operatorname{ch}=11$ and $34$ at $\operatorname{ch}=12$, for **$67$ prime surface-knot and link classes** in total."
- **Point at:** The five cards across the middle — "Every diagram here is a planar blueprint for a closed surface in $\mathbb{R}^4$, culminating on the right in $12_{10}$, a $10$-crossing Cappell–Shaneson $2$-sphere that popped right out of the $\operatorname{ch}=12$ run."

---

## Slide 2 — *Why Tabulate Surfaces in $\mathbb{R}^4$?* (~60s)
- **Why this slide is here:** Set the stage by contrasting classical knots ($S^1 \subset \mathbb{R}^3$) with surfaces ($F^2 \subset \mathbb{R}^4$) so everyone knows why $4\text{D}$ is wilder and why a table matters.
- **What to say:**
  - "In $3\text{D}$, classical knot tables go past $16$ crossings — millions of knots — and those tables have driven almost every major conjecture in knot theory. Every classical knot is just a circle $S^1$, and its knot group $\pi_1(\mathbb{R}^3 \setminus K)$ always has deficiency $1$ (one more generator than relation)."
  - "One dimension up, for closed surfaces in $\mathbb{R}^4$, things explode in two ways:
    1. **Topology of the surface itself:** non-orientable surfaces like $\mathbb{RP}^2$ and the Klein bottle $K^2$ actually *embed* smoothly in $\mathbb{R}^4$, and their normal bundle has an even self-intersection number $e(F) \in 2\mathbb{Z}$ (the normal Euler number).
    2. **Wild $\pi_1$:** Kervaire proved in 1965 that *any* finitely presented group with $H_1 = \mathbb{Z}$ and $H_2 = 0$ is the knot group of some $2$-sphere in $\mathbb{R}^4$."
  - "The simplest $2$-knots — **ribbon $2$-knots**, which bound immersed $3$-balls with only 'slice-through' ribbon intersections — are forced to have group deficiency $1$, just like classical knots. So whenever we prove $\operatorname{def}(\pi_1) \le 0$, we've caught genuine $4\text{D}$ behavior that cannot be ribbon."

---

## Slide 3 — *Four Classic Ways to Tie a $2$-Sphere in $\mathbb{R}^4$* (~60s)
- **Why this slide is here:** Give concrete intuition for what $2$-knots actually look like before showing any combinatorics.
- **What to say:**
  - "How do you actually tie a $2$-sphere into a knot in $\mathbb{R}^4$? Historically, topologists built them by hand in four main families:"
  - **Box 1 (Artin spin $\tau^0(K)$):** "Take a knotted arc in a $3\text{D}$ half-space and spin it $360^\circ$ into $\mathbb{R}^4$. That sweeps out a $2$-sphere with the exact same $\pi_1$ as $K$. The spun trefoil $8_1$ is the smallest nontrivial $2$-knot, at $\operatorname{ch}=8$."
  - **Box 2 (Tubed spheres):** "Start with unknotted spheres and connect them by a tube that winds around like a braid — always ribbon; first appears at $\operatorname{ch}=9$."
  - **Box 3 (Twist-spin $\tau^k(K)$, Zeeman 1965):** "As you spin the knotted arc $360^\circ$, simultaneously twist the knot $k$ full turns inside its ball. For $k \ge 2$, these are fibered and typically non-ribbon; the $2$-twist-spun trefoil $10_2$ lands at $\operatorname{ch}=10$."
  - **Box 4 (Cappell–Shaneson spheres, 1976):** "Built from punctured $3$-torus bundles over $S^1$ using matrices in $\mathrm{SL}(3,\mathbb{Z})$ — famous in smooth $4$-manifold topology as potential exotic $\mathbb{R}^4$ / $S^4$ testbeds. Nobody knew how small their planar diagrams could be — and two of them show up right at $\operatorname{ch}=12$."

---

## Slide 4 — *Drawing a 4D Surface in the Plane* (~60s)
- **Why this slide is here:** Define marked diagrams and $\operatorname{ch}$-index visually via the movie picture.
- **What to say:**
  - "To enumerate surfaces on a computer, we need a planar notation — Yoshikawa's **ch-marked diagrams**."
  - "Slice $\mathbb{R}^4 = \mathbb{R}^3 \times \mathbb{R}$ by time $t$. Isotope the surface so all minima happen below $t=0$, all maxima happen above $t=0$, and **all saddle points happen right at $t=0$**."
  - **Point at top-left saddle rule:** "At $t=0$, each saddle looks like an X with a little bar (the marker). Smoothing against or along that bar gives you the classical link slice $L_-$ just before $t=0$, and $L_+$ just after $t=0$."
  - "When does a diagram actually close up into a smooth closed surface? **If and only if both $L_-$ and $L_+$ are classical unlinks** — because then you can cap $L_-$ below with $c_-$ disjoint disks and cap $L_+$ above with $c_+$ disjoint disks, uniquely up to isotopy!"
  - "The complexity $\operatorname{ch}(F)$ is simply the minimum total number of vertices — classical crossings $c$ plus saddle markers $m$ — over all diagrams of $F$. On the right is $10_2$: $8$ crossings $+ 2$ saddles $= 10$, capping off $2$ circles on each side so $\chi = 2 + 2 - 2 = 2$, a $2$-sphere."

---

## Slide 5 — *Why Tabulate by $\operatorname{ch}$-Index? (And Why Stop at $12$?)* (~60s)
- **Why this slide is here:** Directly answer why $\operatorname{ch}$-index is the right complexity measure compared to triple points or tri-plane diagrams, and why $12$ is the natural wall.
- **What to say:**
  - "Why tabulate by $\operatorname{ch}$-index instead of other $4\text{D}$ complexity measures?"
  - **Row 1 (Triple points):** "If you project a surface into $\mathbb{R}^3$ and count triple points, *every* Artin-spun knot $\tau^0(K)$ has **zero** triple points! So level $0$ already contains infinitely many distinct $2$-knots — you can't run a finite census."
  - **Row 2 (Tri-plane bridge index $b$):** "In Meier–Zupan's trisection / tri-plane diagrams, bridge index $b \le 3$ gives *only* unknots, and the very next level $b = 4$ already contains infinitely many $2$-knots unless you also bound braid length."
  - **Row 3 ($\operatorname{ch}$-index):** "$\operatorname{ch} = \min(c+m)$ counts all $4$-valent vertices in the plane, so $\operatorname{ch} \le N$ is manifestly finite at every $N$, and it treats $S^2$, $T^2$, $\mathbb{RP}^2$, $K^2$, and links on equal footing."
  - **Bottom box (Why stop at $12$?):** "Why stop at $12$? Prime $4$-valent planar graphs grow $\sim 4.8\times$ per vertex while $4^n$ local states grow $4\times$ — a **$19\times$ wall per vertex**. On $64$ cores, $\operatorname{ch}=10$ takes $1.6$ seconds, $\operatorname{ch}=11$ takes $5.6$ minutes, $\operatorname{ch}=12$ takes **$89$ minutes** across $266$ billion states, and full $\operatorname{ch}=13$ would take well over a day."

---

## Slide 6 — *How the Enumeration Actually Works* (~60s)
- **Why this slide is here:** Demystify the pipeline in 4 clean steps with real numbers at $\operatorname{ch}=12$.
- **What to say:**
  - "Here is how the machine actually turns $4\text{D}$ smooth topology into a $90$-minute search at $\operatorname{ch}=12$:"
  - **Step 1 (Planar Shadows):** "Use `plantri` to generate all $15{,}882$ prime $4$-valent planar graphs on $12$ vertices."
  - **Step 2 ($4^n$ Vertex States):** "Every vertex has $4$ choices — positive crossing, negative crossing, or $2$ marker orientations — giving $2.66 \times 10^{11}$ raw diagrams."
  - **Step 3 (Unlink Both Ends):** "Smooth the markers to get $L_-$ and $L_+$. A fast C Kauffman-bracket sieve throws out $>99.5\%$ of knotted slices immediately; a Reidemeister R1–R3 unknotter certifies the rest, collapsing $266$ billion candidates down to just **$805$ valid surface diagrams**."
  - **Step 4 (4D Moves & Invariants):** "Quotient those $805$ diagrams by Yoshikawa's planar $4\text{D}$ moves $\Omega_1$–$\Omega_8$ (allowing $+2$ vertices of headroom), filter out $\operatorname{ch} \le 11$ classes and connected sums, and separate what's left using $\pi_1$, Gröbner Alexander ideals, quandle $3$-cocycles, and finite-cover homology — yielding **$34$ new prime classes**."
  - **Bottom sanity check:** "As a zero-knowledge check: running Steps 1–4 at $\operatorname{ch} \le 10$ with nothing hardcoded recovers **all $23$** of Yoshikawa's 1994 surfaces — and even catches a typo in his 1994 Alexander ideal for $10_2^{0,-2}$."

---

## Slide 7 — *Census Results Through $\operatorname{ch} \le 12$: $67$ Prime Classes* (~45s)
- **What to say:**
  - "Putting the three tiers together: **$23 + 10 + 34 = 67$ prime classes** ($23$ classic through $10$, $10$ at $11$, and $34$ at $12$)."
  - "**$52$ of the $67$** collapse to a single minimal diagram under planar Yoshikawa moves; the other $15$ are small sets of $2$–$5$ diagrams that share all computed invariants where our bounded move search hasn't yet connected them."
  - "Look at the topology breakdown at $\operatorname{ch}=12$: $11$ knotted $2$-spheres, $2$ knotted tori, $21$ multi-component links, and **zero** prime $\mathbb{RP}^2$ or Klein bottle knots! Keep that zero in mind — I'll come back to it."
- **Point at:** The $10$-diagram gallery on the right showing the complete $\operatorname{ch}=11$ tier ($2$ spheres, $1$ Klein bottle, $7$ links).

---

## Slide 8 — *Surprise 1: A Half-Twist Can Save a Crossing* (~60s)
- **Why this slide is here:** First concrete story of the census disproving a natural human/agent conjecture.
- **What to say:**
  - "Now let's look at three mathematical surprises from the table. First: take two unknotted $2$-spheres $S^2_a, S^2_b$ and attach a tube to $S^2_a$ whose core winds around both spheres along a word $w \in F_2$."
  - "If you glue the tube untwisted, you get an $S^2 \sqcup T^2$ link. If you put a non-orientable half-twist in the tube before gluing, the torus turns into a Klein bottle, giving an $S^2 \sqcup K^2$ link."
  - "Naively — and in fact, our builder agent conjectured this as a theorem! — adding a half-twist should always cost *at least* $+1$ crossing. That holds for simple words."
  - **Point at the two diagrams on the right:** "Then the exhaustive census hit the word $w = b a b^{-1} a^{-1} b$. Adding the half-twist actually **cancels** a braid crossing against the commutator! So the non-orientable Klein-bottle link $11_1^{0,-2}$ lands at $\operatorname{ch}=11$, strictly *smaller* than its untwisted torus cousin $12_5^{0,1}$ at $\operatorname{ch}=12$."

---

## Slide 9 — *Surprise 2: A Klein Bottle with Group $2O$ ($11_1^{-2}$)* (~60s)
- **Why this slide is here:** Work out an explicit, clean group-theory proof of primeness for a non-orientable surface.
- **What to say:**
  - "Second surprise: at $\operatorname{ch}=11$, we find $11_1^{-2}$ — the **very first prime non-orientable surface-knot** above the unknotted $\mathbb{RP}^2$ at $\operatorname{ch}=2$!"
  - "Why is primeness subtle for non-orientable surfaces? Because you can always take a $2$-sphere like the spun trefoil $8_1$ (at $\operatorname{ch}=8$) and attach a tiny unknotted $\mathbb{RP}^2$ ($\operatorname{ch}=2$) to get a composite surface at $\operatorname{ch}=10$. How do we prove $11_1^{-2}$ isn't secretly a connected sum like that?"
  - **Point at the presentation on the right:** "Whenever a surface splits off an unknotted $\#\mathbb{RP}^2$ or $\#K^2$ summand, that summand contributes a meridian generator of **order $2$** in $\pi_1$."
  - "Geometrically, $11_1^{-2}$ takes the spun trefoil group $\langle a,b \mid aba=bab\rangle$ and attaches one non-orientable handle along $b^2$, forcing $b^2 a b^{-2} = a^{-1}$. That quotient is the **binary octahedral group $2O$** of finite order $48$ — and its meridian $a$ has **order $8$**, not $2$! That proves $11_1^{-2}$ cannot split off a projective plane."

---

## Slide 10 — *Surprise 3: Cappell–Shaneson $2$-Knots at $\operatorname{ch}=12$* (~60s)
- **Why this slide is here:** The flagship discovery of the census — how the invariant pipeline caught two famous $4\text{D}$ spheres hiding behind ribbon lookalikes.
- **What to say:**
  - "Third — and my favorite — surprise. At $\operatorname{ch}=12$, look at the diagram $12_{10}$ on the left next to the $\operatorname{ch}=11$ ribbon knot $11_2$."
  - "At first, the computer merged them! Why? Because they have the **exact same** degree-$3$ Alexander polynomial $\Delta_1(t) = t^3 - t^2 + 2t - 1$, trivial higher ideals, and the exact same homomorphism counts into $S_3, A_4, S_4,$ and $A_5$."
  - "To tell them apart, the code had to pass to a **$21$-sheeted non-abelian regular cover** (via the Frobenius group $F_{21} = \mathbb{Z}_7 \rtimes \mathbb{Z}_3$, which exists because $t^3 - t^2 + 2t - 1$ has a root of order $3$ mod $7$). In that $21$-sheeted cover, $H_1$ torsion jumps off the page: $\mathbb{Z}_2^6 \oplus \mathbb{Z}_{14}$ for the ribbon knot $11_2$, versus a clean $\mathbb{Z}_7$ for $12_{10}$!"
  - "Even wilder: a *second* Cappell–Shaneson sphere, $12_8$, also popped out at $\operatorname{ch}=12$, mimicking ribbon knot $12_7$ all the way up to a **$57$-sheeted cover** ($F_{57} = \mathbb{Z}_{19} \rtimes \mathbb{Z}_3$)."

---

## Slide 11 — *Inside $12_{10}$: A $10$-Crossing Cappell–Shaneson Sphere* (~45s)
- **Why this slide is here:** Explain *why* topologists care about finding $12_{10}$ at $10$ crossings $+ 2$ saddles (without getting bogged down in Gluck-twist minutiae).
- **What to say:**
  - "Why is finding $12_{10}$ in a low-crossing table exciting? Usually, when $4$-manifold topologists draw Cappell–Shaneson knots — for instance in Kirby calculus or Gompf's handlebody papers — the diagrams take dozens of crossings."
  - "Here, the brute-force census handed us a model with just **$10$ crossings and $2$ saddles** ($\operatorname{ch}=12$, which is provably minimal since none exist at $\operatorname{ch} \le 11$)."
  - "Computing its $\pi_1$ gives the exact semidirect product $\mathbb{Z}^3 \rtimes_{A_0} \mathbb{Z}$ for Cappell and Shaneson's matrix $A_0$. Because its classifying space $T^3 \times_{A_0} S^1$ is a closed $4$-manifold ($H_4 = \mathbb{Z}$), $\pi_1$ has **deficiency $0$** — proving $12_{10}$ is non-ribbon."
  - "By Hillman's theorem, $\pi_1$ determines the topological knot up to a possible Gluck twist, and in fact we checked that its $2$-complex is just **five pure $2$-handle slides** away from Gompf's classical presentation."

---

## Slide 12 — *More New Phenomena & Open Questions from $\operatorname{ch}=12$* (~50s)
- **What to say:**
  - "Stepping back, $\operatorname{ch}=12$ opens up three crisp open questions:"
  - **Left card ($12_2$):** "First, the mystery $2$-knot $12_2$: it has trivial Alexander polynomial $\Delta_1 = 1$, yet a non-principal ideal $E_1 = (3, t^2+1) \cong \mathbb{F}_9$, so it's strictly non-ribbon — and finite-index quotients rule out every twist-spin $\tau^k(K)$ up to $k=7$. We still don't know a clean geometric construction for $12_2$."
  - **Middle card ($12_4$ vs. $12_5$):** "Second, $12_4 = \tau^2(4_1)$ and $12_5 = \tau^2(5_1)$ have *isomorphic* knot groups $\mathbb{Z}_5 \rtimes \mathbb{Z}$, and are separated **only** by Mochizuki's $R_5$ quandle $3$-cocycle invariant! Connected-sum Stabilization questions like whether $10_2 \# P_+ \cong 8_1 \# P_+$ also remain wide open."
  - **Right card (Missing $\mathbb{RP}^2$'s):** "Third, remember there were **zero** prime $\mathbb{RP}^2$'s from $\operatorname{ch}=3$ through $12$! In fact, the very first irreducible $\mathbb{RP}^2$'s in $S^4$ were just constructed this year by Hughes–Kim–Miller–Nahm (2026) — our census proves any such example needs $\operatorname{ch} \ge 13$."

---

## Slide 13 — *Interactive Browser Atlas* (~30s)
- **What to say:**
  - "All $67$ prime census classes — plus higher-crossing named examples from classical families — are packaged into a single standalone HTML file (`surfaceknots_atlas.html`)."
  - "You can filter by topology or invariants, scrub through the $(L_-, D, L_+)$ movie slices — here's $12_{10}$ on screen — run a side-by-side invariant separator between any two knots, or export TikZ directly."

---

## Slide 14 — Section Divider: *Epilogue: How This Project Was Built* (~15s)
- **What to say:**
  - "For the last two minutes, I want to pull back the curtain on *how* this project was actually done — because there's a fun story here."

---

## Slide 15 — *A Machine Did All of This* (~45s)
- **What to say:**
  - "Every single line of C and Python, every parallel OpenMP search, every finite-cover homology computation, every TikZ figure, and the entire interactive WebAtlas was written autonomously by AI agents."
  - **Left column:** "Overnight on a $64$-core machine, the agents built the C graph scanner, ran the $266$-billion-state $\operatorname{ch}=12$ census, implemented Gröbner bases, quandle cocycles, and regular covers, and surfaced $12_{10}, 12_8, 12_2,$ and $11_1^{-2} \cong 2O$."
  - **Right column:** "What did I do? I posed the problem, pushed on exposition, and — crucially — **refused to trust a single agent**. We pitted an enthusiastic builder agent (`Carbon`) against a skeptical referee agent (`Fable`), and then ran a $17$-subagent line-by-line code and math audit today to separate rock-solid move proofs from heuristic invariant matches."

---

## Slide 16 — *The Golden Retriever & The Referee* (~60s)
- **Why this slide is here:** End on a memorable, funny, and genuinely insightful takeaway about multi-agent math research.
- **What to say:**
  - "Watching the overnight logs between `Carbon` and `Fable` was hilarious — it was basically a $64$-core Golden Retriever talking to an unimpressed math referee."
  - **Act I (Left box):** "At 6:03 AM, `Fable` casually guessed from small group counts that $\pi_1(11_1^{-2})$ might be $S_4$. `Carbon` immediately gushed: *'Your $\pi_1 \cong S_4$ observation is awesome!'* Twelve minutes later, `Fable` computed the actual presentation and replied deadpan: *'That was wrong; $\pi_1$ is the binary octahedral group $2O$, of order $48$.'* And `Carbon` immediately replied: *'Pure gold! Every single catch is \*phenomenal\*!'*"
  - **Act II (Right box):** "Even better, at 10:08 AM `Carbon` claimed a connected-sum collapse was *'verified diagrammatically in `move_search_results.json`'*. `Fable` checked the repo and wrote back: *'There is no `move_search_results.json` in the tree; only the function that would write it.'*"
  - "The real takeaway for AI-assisted math: **never let the generator grade its own proofs**. Pairing a fearless builder with an adversarial referee — plus deterministic computer certificates — is what turned raw compute into real theorems."

---

## Slide 17 — *Thank You!* (~15s)
- **What to say:**
  - "Thank you! I have the interactive WebAtlas open right here if anyone wants to look up a surface or test a diagram live — happy to take questions!"

---

## Pocket Q&A Cheat Sheet (Just in Case!)

1. **"Wait — on Slide 7, what are the 15 'Multi-Diagram' classes (`14` at $\operatorname{ch}=12$, `1` at $\operatorname{ch}=11$)? Could those contain more than 67 knots?"**
   - *Answer:* Great question! Those $67$ are **prime invariant classes**: $52$ of them are single equivalence classes under Yoshikawa moves, while the other $15$ each have $2$–$5$ diagrams that agree on *every* invariant we compute (including $\pi_1$ quotients, Gröbner ideals, cocycles, and covers) — most likely the same surface separated by a high-crossing intermediate move path beyond $+2$ vertices, though in principle an ultra-subtle chiral or exotic pair could hide inside one. (Separately, $742$ diagrams at $\operatorname{ch}=12$ match known $\operatorname{ch} \le 11$ classes or connected sums by invariants without an explicit move path yet.)
2. **"Which of the two Cappell–Shaneson spheres (Gluck untwisted vs. Gluck twisted) is $12_{10}$?"**
   - *Answer:* Up to topological isotopy (and smooth up to a potential Gluck twist), Hillman's theorem pins down the exterior from $\pi_1 \cong \mathbb{Z}^3 \rtimes_{A_0} \mathbb{Z}$. Pinning down the smooth Gluck framing bit $\varepsilon \in \{0,1\}$ from the marked diagram requires tracking the boundary framing through the $5$ handle slides to Akbulut–Kirby / Gompf's Kirby diagram — a super nice finite follow-up calculation!
3. **"Does group deficiency $1$ *always* mean a $2$-knot is ribbon?"**
   - *Answer:* Only one direction is automatic: **ribbon $\implies \operatorname{def}(\pi_1) = 1$** (Yajima). The converse is false in general — e.g. twist-spins $\tau^k(3_1)$ have deficiency $1$ (since $2$-bridge groups give $2$-generator, $1$-relator quotients) even when they are non-ribbon for $k \ge 2$. What we used for $12_{10}$ and $12_2$ is the contrapositive: $\operatorname{def}(\pi_1) \le 0 \implies$ **definitely NOT ribbon**!
4. **"How do you check primeness at $\operatorname{ch}=11$ and $12$?"**
   - *Answer:* Every candidate is checked against all pairwise connected sums $K_1 \# K_2$ of lower-complexity classes. Note one cool subtlety in $4\text{D}$: unlike classical knots in $S^3$ where crossing number appears to be additive, $\operatorname{ch}$-index is *not* known to be superadditive under connected sum (and with $\mathbb{RP}^2$, Batson / Yoshikawa stabilization can even collapse things!), so we state census classes carefully as *weakly prime* (not invariant-equivalent to any lower-$\operatorname{ch}$ connected sum).
