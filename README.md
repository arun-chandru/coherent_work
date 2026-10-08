# Boundary coherence and work fluctuations under quantum confinement control

R. Arun Chandru — 8 October 2026


## Extended abstract

A warm particle can occupy many different states in a box. Near a potential
wall, however, those states can respond in almost the same direction when
the wall is changed. This paper studies an interference experiment in which
an auxiliary two-level system coherently chooses between leaving the wall
unchanged and raising it. Keeping one auxiliary outcome can retain a nearly
pure particle state even when the initial thermal state is becoming more
mixed. The successful outcome is rare, and the paper counts the preparations
needed to obtain it; this is conditional purification, not free cooling.

One boundary source connects this experiment to the work involved in changing
the wall. Its spectral measure determines both the interference signal and
the leading division of dissipated work into changed energy populations and
quantum coherence. The division depends on a precisely derived competition
between the final equilibrium weights and the boundary-generated population.
Smooth potential walls are included in this result.

Work measurements reveal another effect that the limiting interference signal
does not capture by itself. For a sharp wall, transitions between nearby
energy levels accumulate in the average absolute work. Their positive and
negative contributions nearly cancel in the mean signed work. Along specified
thermal paths, the ratio of absolute to net work grows without bound. The
relevant energy gaps are small compared with the barrier energy; they need
not be small in laboratory units. Every finite system considered has finite
first and second work moments.

The paper also supplies two concrete operational checks. A full-continuum
calculation for an analytic smooth wall certifies an increase in purity from
below 0.386060 to above 0.9139989 under the stated timing, visibility, and
detector errors, with fewer than 389 preparations per accepted output on
average. A separate reconstruction design uses finitely many Ramsey times
and exact count confidence intervals. Its very large preparation costs are
reported explicitly. These are theoretical predictions and reproducible
mathematical certificates, not results of a laboratory experiment.

- `r_arun_chandru_boundary_coherence_and_work.tex`: complete source.
- `r_arun_chandru_boundary_coherence_and_work.pdf`: verified 23-page article
  and appendices in one file.
- `verification/`: portable scripts, numerical certificates, dependency pins,
  integrity hashes, and detailed replay instructions.

The model is a one-particle Dirichlet confinement problem with coherently
controlled local potentials. Its asymptotic regimes, conditional success
probabilities, physical energy/time conversion, and specified error channels
are part of the statements.

The infinite-dimensional proofs are analytic. Numerical verification is used
only for the explicitly finite operating points and statistical budgets.
The sharp-wall absolute-work logarithm is not silently extended to smooth
walls. The smooth-wall thermodynamic and purification conclusions are proved
separately for their declared profile class.

## Compile

Use a LaTeX distribution providing REVTeX 4.2 and standard AMS, fontenc,
graphicx, booktabs, and hyperref packages. From this directory:

```text
pdflatex -interaction=nonstopmode -halt-on-error r_arun_chandru_boundary_coherence_and_work.tex
pdflatex -interaction=nonstopmode -halt-on-error r_arun_chandru_boundary_coherence_and_work.tex
```

The bibliography is embedded; no BibTeX run or external `.bib` file is needed.
Rerun LaTeX if it asks for reference resolution. The date is supplied by
`\today` and therefore follows the compilation date.

The delivered PDF was built with Tectonic 0.17.0 and REVTeX 4.2e. The usual
`pdflatex` commands above are provided for standard TeX distributions;
there was no local `pdflatex` executable for a second-engine replay.
The source includes a conditional small-size Latin Modern font fallback
for incomplete caches. Complete font installations are left unchanged.

## Reproduce the finite results

See [verification/README.md](verification/README.md) for the exact serial
commands and [verification/requirements.txt](verification/requirements.txt)
for the replayed dependency versions. It includes untrusted floating root
proposals for speed; signs, labels, norms, and spectral tails are rechecked
rigorously on every replay.

All reported confidence guarantees refer to independent prospective trials.
No synthetic outcomes are presented as experimental observations. Preparation,
ancilla reset, switching, readout, and apparatus overheads remain additional
to the explicit Hamiltonian evolution duration. Regenerate certificates after
changing their source files; source hashes are checked by the dependent runs.