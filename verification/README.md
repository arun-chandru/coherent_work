# numerical verification

This directory reproduces the manuscript's finite numerical examples. All local
Python imports resolve to the seven files here. Only Python and the three
external packages in `requirements.txt` are required; no research workspace,
external data, or bundled runtime is part of the package. The version pins record
the replay environment, with Python 3.12.

The certificates describe mathematical model calculations and preparation
budgets. The proofs and definitions are in the accompanying manuscript.

## Installation and replay

Copy the entire directory to any writable location. Create and activate a Python
virtual environment, then run `python -m pip install -r requirements.txt`.
Run the commands below from that copied directory, in their displayed order.
One CPU job at a time is sufficient. A full run can take several minutes.

### Smooth erf profile, then jitter and detector extension

```text
python HOT_ERF_SMOOTH_CHECK.py --output HOT_ERF_SMOOTH_CERTIFICATE.json
python HOT_ERF_NOISE_JITTER_CHECK.py --output HOT_ERF_NOISE_JITTER_CERTIFICATE.json
```

The first command recomputes certified root brackets, Sturm labels, signed
overlap integrals, normalization, and complete initial and final tail budgets.
The two `HOT_ERF_FLOAT_PROPOSALS_*.json` files are untrusted floating root
proposals used to accelerate this step. Every proposed bracket is independently
checked with interval arithmetic and an exact nodal label; the proposals are
not certificates. Without these files the checker regenerates them.

The baseline includes the actual analytic erf wall and its own complete Gibbs
state. Its timing extension covers each fixed time. The second command checks
the baseline's source hashes and extends it to an arbitrary unrecorded duration
law supported in the stated time interval, independent of the prepared state,
with the specified visibility and detector errors. Consult that second output
for the manuscript's jitter/readout claim.

### Sharp-wall finite reconstruction and preparation counts

```text
python check_hot_finite_count_certificate.py --output hot_finite_count_certificate.json
python check_hot_ramsey_count_statistics.py --certificate hot_finite_count_certificate.json --output hot_ramsey_count_statistics_first.json --epsilon 1/20
python check_hot_finite_count_certificate.py --inverse-r 200 --s 128 --initial-count 64 --output hot_finite_count_followup_certificate.json
python check_hot_ramsey_count_statistics.py --certificate hot_finite_count_followup_certificate.json --output hot_ramsey_count_statistics_followup_point.json --epsilon 1/20
python check_hot_ramsey_count_statistics.py --certificate hot_finite_count_followup_certificate.json --output hot_ramsey_count_statistics_followup_conservative.json --epsilon 1/200
```

The two deterministic runs use cap 1, horizon 32, and 256 midpoint times. The
first device is `(r,s)=(1/50,12)`; the followup is `(1/200,128)`. Computational
mode prefixes do not replace the full Gibbs preparation: omitted initial and
final modes and source-energy tails enter the output budgets. The statistics
checker reads those deterministic enclosures and records the exact rational
distribution used to choose the time label of each prospective trial.

The statistics commands above supply no observations. If independent trial
counts have actually been collected for this design, the checker additionally
accepts `--trials N --count k`, with integer `N>=1` and `0<=k<=N`; it then adds
an exact-binomial confidence enclosure for those user-supplied observations.
Its positive-sum inversion is intended for rare counts; it rejects more than
100,000 successes after reducing by binomial symmetry.

## Files and integrity

The four entry scripts are the two capitalized smooth-wall checks and the two
finite-count checks. The other three Python files are intact support modules:
`check_hot_ramsey_finite_certificate.py`,
`check_phy_info_finite_purity_certificate.py`, and
`check_phy_info_continuum_certificate.py`. Their older standalone entrypoints
are preserved to avoid rewriting audited mathematical functions; they are not
additional numerical claims of this package and need not be run.

`SOURCE_SHA256SUMS.txt` records SHA-256 hashes of the portable scripts and the
untrusted proposal files. The regenerated deterministic certificates additionally
record hashes of their actual source dependencies; the jitter certificate hashes
its baseline and extension; the statistics outputs hash their deterministic
input and checker. Source paths in those fields are local to this directory.
Elapsed run times are diagnostics and can differ between machines.

The portable copies change import paths and hash/provenance paths, plus the
smooth baseline's default output filename. Mathematical function bodies are
unchanged. Regenerate certificates after any source modification; a stale
certificate's hash cannot establish that the modified source was run.

The examples assume the exact model and channels stated in the manuscript.
Physical preparation, switching, reset, readout, and apparatus overheads are
additional to the reported fresh-preparation and controlled-duration budgets.
