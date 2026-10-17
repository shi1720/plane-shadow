# Equal Zeta Functions, Unequal Rational Planes

Research manuscript by Shivam Gupta (17 October 2026).

For any finite field F_q of characteristic at least five with q = 1 mod 3 and any noncube a, the diagonal cubic fourfolds with coefficient vectors

- X: (1,1,1,1,a,a²)
- Y: (1,1,a,a,a²,a²)

have the same entire zeta function and the same zeta function of their varieties of lines, but 0 and 27 rational planes. The paper gives proofs and exact computations. Over F_7, this is the unique collision among ten diagonal projective classes.

[Read the paper](paper/main.pdf)

This is a research manuscript, not an accepted publication. Classical character-sum identities and the Galkin–Shinder relation are credited in the paper. We do not claim that diagonal point counting, the Fermat plane configuration, or the general incompleteness of zeta functions is new. The literature review did not locate this explicit comparison; it cannot establish universal absence of prior art.

## Reproduce

Python 3.9+ standard library, and a C++17 compiler:

```sh
python3 experiments/run.py
python3 experiments/classify.py
c++ -O2 -std=c++17 tests/grassmannian.cpp -o /tmp/plane-search
/tmp/plane-search results/grassmannian.csv
c++ -O2 -std=c++17 tests/lines.cpp -o /tmp/line-search
/tmp/line-search results/lines.csv
python3 experiments/verify.py
```

For figures install `matplotlib` and run `python3 experiments/figures.py`. To build the manuscript run `cd paper && tectonic main.tex`.

## Evidence

- 729 labeled cubic-class vectors.
- 150 histogram point counts, independently compared to exact Jacobi-sum formulas.
- 20 direct projective enumerations over F_7 and F_13.
- 22 complete Grassmannian searches with pure-cube pruning. Reported total chart sizes are covered spaces, not executed full-subspace tests.
- 6,534 explicit plane matrices checked by independent polynomial substitution and RREF uniqueness.
- 2 independent direct line searches over F_7: 2,943 lines on either fourfold.
- Exact degree-23 middle polynomials for ten classes at eleven primes.

All checks pass. Raw outputs and protocol are in `results/` and `experiments/protocol.json`. Wall-clock timings are descriptive, and no speedup benchmark is claimed. The all-extension result is proved, not extrapolated from finite tests.

## Scope

The examples differ in rational plane descent. They become Fermat cubics over a cubic extension. They do not establish different rationality, different geometric plane configurations, or a classification of all cubic fourfolds. The F_7 census is limited to diagonal equations. No patentability, market size, or commercial revenue claim is established.

Code is MIT licensed; manuscript copyright remains with the author.
