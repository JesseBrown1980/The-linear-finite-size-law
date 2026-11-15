# Prior-work search, November 15, 2026

Status of the law's priority: **PRIORITY_UNVERIFIED.** A scoped web search found no earlier statement of this law; that is a search result and not a proof of absence. A reader who knows
the field is the next step.

## What was searched

Four web searches by seat ACER-CLAUDE-FABLE52 on November 15, 2026:

1. rigorous explicit lower bounds on the spectral gap of finite-group (S₃) lattice gauge Hamiltonians on a torus, with a coupling threshold and its dependence on size;
2. explicit constants in the stability of topological order (toric code, quantum doubles);
3. the quantum double D(S₃) with an electric-field perturbation: finite-size gap thresholds;
4. Perron–Frobenius / Barta bounds with a product ansatz for lattice gauge gaps at finite volume.

## What was found

| Work | What it has | Relation to this law |
|---|---|---|
| Bravyi, Hastings, Michalakis, [arXiv:1001.0344](https://ar5iv.labs.arxiv.org/html/1001.0344); Bravyi, Hastings, [arXiv:1001.4363](https://www.arxiv-vanity.com/papers/1001.4363/); Michalakis, Zwolak, [arXiv:1109.1588](https://pdf.arxiv.org/pdf/1109.1588) | A perturbation threshold that does not depend on size, for commuting topologically ordered Hamiltonians, quantum doubles included | Stronger in kind (size-independent), with constants stated as existing and no number printed. This law is explicit and its threshold grows with the side |
| Lemm et al., [arXiv:2308.01405](https://arxiv.org/pdf/2308.01405), and finite-size criteria such as Kastoryano, Lucia, [arXiv:1705.09491](https://arxiv.org/pdf/1705.09491) | Explicit, size-uniform gap bounds for frustration-free models by finite-size criteria | Different models and method; the perturbed S₃ Hamiltonian here is not frustration free |
| A Rydberg-based D(S₃) construction, [Phys. Rev. B 99, 035112](https://link.aps.org/accepted/10.1103/PhysRevB.99.035112) (as summarised by the search) | A flux gap provable in the thermodynamic limit; the charge gap expected and described as hard to show rigorously | Adjacent: the same group, another Hamiltonian |
| Horn and others, [Phys. Rev. D 38, 2591](https://doi.org/10.1103%2FPhysRevD.38.2591) | Rigorous upper bounds on the mass gap for U(1), SU(2), SU(3) in 2+1 dimensions from an exactly known ground state | Upper bounds, continuous groups |
| Strong-coupling series, [Nucl. Phys. B (1981)](https://inspirehep.net/literature/164289) | Series for the mass gap of SU(2), SU(3), Z₂, Z₃ | Expansions, not certified bounds |

No result found states an explicit coupling threshold, linear in the side, for the gap above the flat cluster of the S₃ Kogut–Susskind torus.

## What this does and does not show

- The quadratic law this replaces, r ≥ 4L², is a step of this project's own proof plan, not a published result.
- Four searches are a small sample. Results in books, in theses, or phrased in other terms (relative boundedness, cluster expansions for classical ground states with degeneracy) can
  have been missed. The summaries above are the search tool's readings of the pages and were not checked against the papers line by line.
- The two papers of Bravyi, Hastings and Michalakis were read from their TeX sources for the hypotheses and the constants (see the Yang–Mills repository, draft pull request 201).
