## Documentation Quality & Reproducibility

This section explains how the card is organized and why, how to read the results, and how to check, rerun or extend every number. Everything described here refers to files in this folder.

### 1. How the card is organized, and why

The card is built so that a reader can stop at any depth and still hold a correct picture.

| If you have… | Read | You will learn |
|---|---|---|
| 2 minutes | *Summary* and its main findings | What was measured, the headline results, and what they do and do not show |
| 10 minutes | + *Novelty*, *Technical Execution*, *Results* and Figures 1–3 | How the problem is encoded, how the penalty A changes behaviour, and the trade-off |
| 30 minutes | + this section | How to verify, rerun or extend every number |

After reading, a reader should be able to (1) formulate a budgeted selection problem as a QUBO and a Pauli Hamiltonian, (2) compute the critical penalty A_crit and explain what happens below and above it, and (3) rerun the study or apply it to another instance.

| Design choice | Reason |
|---|---|
| Order: problem → mathematics → circuit → metrics → results → trade-off → guidance → limits | The derivation is written before the experiment, so each result can be traced back to a formula |
| A worked, verifiable chain: QUBO → coefficients lᵢ, qᵢⱼ → Pauli terms hᵢ, wᵢⱼ → `SparsePauliOp` | A reader can follow the algebra and the script checks it against all 2^N energies on every run |
| Every metric is defined in a table and reported next to baselines (uniform superposition, random feasible pick, greedy by confidence) | A good-looking number (for example a high approximation ratio on very few feasible shots) is never mistaken for useful behaviour |
| The instance is built so that greedy selection is a trap | The biological meaning of J_ij (synergy and redundancy) becomes visible in the answer, not just in the formulas |
| Results come with a trade-off table by region of A, a practical checklist, and an explicit "what this does not show" list | The reader leaves with guidance and with its limits |
| The same names are used in text, tables, figure legends, CSV columns and code (`A`, `p`, `feasibility_pct`, `AR_bio`, `P_opt`) | Any number can be cross-referenced without translation |

**Glossary.**

| Term | Meaning in this project |
|---|---|
| QUBO | Quadratic Unconstrained Binary Optimization: minimize a quadratic function of binary variables |
| Penalty parameter A | Weight of the term A (Σxᵢ − k)² that punishes selecting a number of proteins different from k |
| A_crit | Smallest A above which the ground state of H is the constrained optimum (1.000 for this instance) |
| Normalized gap | (energy of the best non-optimal state − energy of the optimum) / (E_max − E_min); negative if the optimum is not the ground state |
| QAOA depth p | Number of alternating cost and mixer layers (p = 1 and 2 are studied) |
| Restart | One independent COBYLA optimization from random initial angles (8 per setting) |
| Cost-selected run | The restart with the lowest final cost; it can be chosen without knowing the optimum |
| Feasibility rate | % of shots whose bitstring has Hamming weight (number of ones) equal to k |
| AR_bio | Mean objective value of the feasible shots, normalized between the worst and best feasible solution |
| AR_overall | Feasibility × AR_bio (infeasible shots count as zero) |
| P_opt | Probability of sampling the global optimum |
| Expected shots to optimum | 1 / median P_opt |

### 2. Reading the results

**Reading order for the main table** (in *Results* and in `results/exact/summary.md`).

| Column | Question it answers |
|---|---|
| A / A_crit and normalized gap | Where are we relative to the critical penalty, and how much signal is left? |
| Feasibility % | Do the shots respect the budget? |
| AR_bio | Among feasible shots, how good are they (0 = worst feasible, 1 = optimum)? |
| AR_overall | What is the quality if infeasible shots count as failures? |
| P_opt % | How often does a shot hit the optimum? |

**Figure guide.**

| Figure | Question | What to look at | Caveat |
|---|---|---|---|
| 1. Network | What is being selected? | Orange = optimum {P3, P4, P5}; black outline = greedy subset; green edges = synergy, red dashed = redundancy | Synthetic instance, not experimental data |
| 2. Feasibility and ratios vs A | What does A buy? | Blue (feasibility) rises to a plateau after A_crit while red (AR_bio) stays on the dotted random-feasible baseline | At A = 0.25 AR_bio looks high only because ~5 % of shots are feasible; read AR_overall (orange) |
| 3. Histograms | What does the sampled distribution look like? | At A = 0.25 the mass sits on infeasible states; at A = 7.5 it is spread over the 20 feasible states and the optimum (orange) is not singled out | Shows one cost-selected restart at p = 2; the two A values are chosen by median AR_overall |
| 4. Gap vs A (supplementary) | How much signal survives? | The normalized gap falls as A grows; negative below A_crit | Computed exactly by enumeration (small N only) |

**Where each number comes from.**

| Claim in the card | Source file |
|---|---|
| Main results table (exact simulation) | `results/exact/summary.md` |
| Noise comparison | `results/noisy_generic/summary.md` |
| Ablation without rescaling | `results/ablation_no_normalize/summary.md` |
| Instance, A_crit, classical optimum, greedy and random baselines, normalized gaps | `results/exact/run_metadata.json` (keys `instance`, `classical`, `normalized_gap`) |
| Every restart's metrics | `results/*/metrics.csv` |
| Measured distributions and optimal angles (used by Figure 3) | `results/*/records.jsonl` |
| Figures | `results/exact/fig1_network.png` … `fig4_gap_vs_A.png` (same four for the other runs) |

### 3. Reproducibility

**Environment (reference runs).** Python 3.12.3, Qiskit 1.4.6, Qiskit Aer 0.17.2, NetworkX 3.6.1, NumPy 2.4.4, SciPy 1.17.1, Matplotlib 3.10.8. See `requirements.txt` and `results/*/run_metadata.json`.

**Four ways to engage, from cheapest to most complete.**

| Level | Command | Time | What you get |
|---|---|---|---|
| 0. Inspect | Open `results/` | Minutes | Tables, figures and per-run data, no installation needed |
| 1. Rebuild tables and figures | `python qaoa_penalty_sensitivity.py --replot --outdir results/exact` | Seconds | Regenerates `summary.md` and the figures from `records.jsonl` |
| 2. Smoke test | `python qaoa_penalty_sensitivity.py --quick` | ~30 s | 3 values of A, 2 restarts, p = 1 and 2 |
| 3. Full reference runs | See below | ~7 min, ~2 min, ~2 min | The three result sets |

```bash
pip install -r requirements.txt

python qaoa_penalty_sensitivity.py --outdir results/exact                                # ~7 min (413 s measured), 1 CPU core
python qaoa_penalty_sensitivity.py --backend aer-noise --outdir results/noisy_generic     # ~2 min (132 s); needs qiskit-aer
python qaoa_penalty_sensitivity.py --no-normalize --outdir results/ablation_no_normalize  # ~2 min (108 s)
```

**Determinism.** All random choices are seeded (`--seed 7`, with a separate generator per A, p and restart), the exact backend has no sampling noise in the cost, and the final shot distribution uses a seeded sampler. In a check across sessions, a `--quick` run reproduced all 12 matching (A, p, restart) entries of the delivered full run exactly (evaluation counts, feasibility and final cost); timings are not reproducible.

**Checks built into the script.**
- Before every run, the Pauli Hamiltonian is compared against the QUBO energies of all 2^N bitstrings for every A in the grid; the script stops if the difference exceeds 10⁻⁹.
- The exact optimum, the greedy baseline, A_crit, the uniform-superposition feasibility and the random-feasible AR_bio are computed by enumeration and printed at the start.
- The restart used for the histogram is selected by cost, without access to the optimum.

**Output files (per run folder).**

| File | Content |
|---|---|
| `summary.md` | Results table |
| `metrics.csv` | One row per restart (`A`, `p`, `restart`, `nfev`, `time_s`, `cost_final`, `E_H`, `AR_H`, `feasibility_pct`, `AR_bio`, `P_opt`, exact-probability versions) |
| `records.jsonl` | Full records, including the measured distribution and optimal angles |
| `run_metadata.json` | Arguments, versions, wall time, instance, classical results, normalized gaps |
| `fig1_network.png`, `fig2_feasibility_ratio_vs_A.png`, `fig3_histogram.png`, `fig4_gap_vs_A.png` | Figures |
| `qollab_validation_records.jsonl` (project root) | Cost-selected restart for every (A, p), with optimal angles, for single-job validation |
| `results_qollab/validation_results.json`, `fig5_backend_validation_A<A>_p<p>.png` | Output of `--validate-on-backend` |
| `qollab_sampler.py` (project root) | Optional plug-in for `--backend custom` |

**Backends: what was tested.**

| Option | Status |
|---|---|
| `--backend statevector` (default) | Tested end-to-end; used for the exact results |
| `--backend aer-noise` | Tested end-to-end (generic depolarizing noise, **not** IonQ-calibrated) |
| `--backend ionq` | Implemented from IonQ's documentation (`qiskit-ionq` < 1.0 for Qiskit 1.x, `IONQ_API_KEY`, `--ionq-noise-model`). **Not tested** |
| `--backend qollab` | Uses the global `backend` that Qollab injects (no API key). One sampler call is allowed when the injected backend is an IonQ simulator (see below). Tested only with a *locally simulated* injected backend (Aer with generic noise), **not on Qollab itself** |
| `--backend custom --sampler-factory qollab_sampler.py:make_sampler` | Plug-in file `qollab_sampler.py` is included. Its `local` mode was tested end-to-end; its `qollab` mode behaves like `--backend qollab` (tested with the simulated injected backend); its `ionq-simulator` mode is untested |

**Running on Qollab (one job per run).** What follows is based on the platform information supplied by the author; it was not checked against Qollab's own documentation. In the Code Playground / Python projects the selected simulator is injected as a global `backend` (chosen in the *Run on Qollab / Select QPU* menu), no API key is needed, and the IonQ simulator's noise model is a string set with `backend.set_options(noise_model="forte-1")` (names reported: `forte-1`, `forte-enterprise-1`, `aria-1`, `aria-2`). Simulators are free; the remote IonQ simulators accept **one job per run**, and a single `backend.run([...])` call counts as one job even when it carries several circuits. The base libraries are `qiskit` 1.x, `numpy`, `scipy`, `matplotlib` and `networkx`.

Every COBYLA evaluation is one job, so the optimization cannot run on the Qollab backend. The study is therefore split in two:

| Step | Where | What |
|---|---|---|
| 1. Optimize | Local machine, exact statevector | `python qaoa_penalty_sensitivity.py --outdir results/exact` produces the full study |
| 2. Select | Already done | `qollab_validation_records.jsonl` keeps, for every (A, p), the restart selected by cost (22 records, with optimal angles) |
| 3. Validate | Qollab (needs only Qiskit, NumPy, SciPy and Matplotlib) | One sampler call (one circuit, one job) runs the stored angles for one value of A on the selected backend and compares the result with the noiseless distribution |

The Code Playground runs the file as a script and does not accept command-line arguments, so the script falls back to the list `QOLLAB_DEFAULT_ARGS` at the top of the file whenever it finds an injected `backend` and no arguments (A = 7.5, p = 2, 8192 shots). Without an injected `backend` (a normal local run) the same defaults reproduce the full study. To run step 3 upload `qaoa_penalty_sensitivity.py` and `qollab_validation_records.jsonl` to the project (uploaded files sit in the working-directory root, where the script looks for them), select the IonQ simulator and its noise model in the run menu and press Run. Each run validates one value of A, so edit `QOLLAB_DEFAULT_ARGS` and run again for others (for example A = 0.25). Outside the playground, or if arguments ever become available:

```bash
python qaoa_penalty_sensitivity.py --backend qollab --validate-on-backend \
    --records qollab_validation_records.jsonl --validate-A 7.5 --validate-p 2 \
    --final-shots 8192 --qollab-noise-model forte-1
```

The output is a table (feasibility, AR_bio, P_opt, AR_H: noiseless versus backend), `results_qollab/validation_results.json` and `results_qollab/fig5_backend_validation_A<A>_p<p>.png`. Before running anything, the script checks that the stored angles reproduce the stored cost for the instance and normalization in use, and, when the injected backend is an IonQ simulator (detected from its name, e.g. `ionq_simulator`), it stops with an explicit "job limit reached" error if a run attempts a second sampler call. Other injected backends are not limited unless `--max-jobs` is set.

*Status.* Tested only with a locally simulated injected backend (Aer, generic depolarizing noise of 2·10⁻⁴ per single-qubit gate and 5·10⁻³ per CX, **not** IonQ's). For A = 7.5, p = 2 the feasibility rate fell from 93.0 % (noiseless) to 79.3–79.7 % across four repeated simulated runs, P_opt from 4.2 % to about 3.4 %, and AR_bio stayed near 0.44–0.45 (0.446 noiseless), with 8192 shots. These numbers illustrate the procedure and say nothing about IonQ. **Validation results from Qollab itself are not yet included**; paste the table printed by the Qollab run here when available. Confirmed by the author, not independently checked: several circuits sent in one call count as one job (this study sends one circuit per run); the injected `backend` exposes the standard Qiskit configuration, so `transpile(circuit, backend=backend)` works (if it ever fails, the script prints a warning and falls back to the {rz, sx, x, cx} basis); the noise model is set with `backend.set_options(noise_model="forte-1")`, `forte-1` being the recommended default.

**Included plug-in.** `qollab_sampler.py` is the optional file for `--backend custom`. Its mode is set with the `QOLLAB_MODE` environment variable (`local` by default, `ionq-simulator`, `qollab`); see the comments inside the file. The built-in `--backend qollab` option needs no plug-in.

```bash
python qaoa_penalty_sensitivity.py --backend custom --sampler-factory qollab_sampler.py:make_sampler
```

For any remote backend outside Qollab, every COBYLA evaluation is one job: start small, for example `--backend ionq --restarts 2 --maxiter 60 --A-values 1 3 10 --depths 1 --shots 1024 --final-shots 1024`, and check queue time and cost first.

**Troubleshooting.**

| Symptom | Cause and fix |
|---|---|
| Qiskit 2.x got installed | `qiskit-ionq` ≥ 1.0 pulls Qiskit 2.x. Reinstall `qiskit<2` and `qiskit-ionq<1.0` |
| `ImportError` with `--backend aer-noise` | Install `qiskit-aer` (`<0.18` was tested) |
| Low feasibility with `--no-normalize` | Expected: the energy scale grows with A (see the ablation in *Results*) |
| COBYLA stops after few evaluations with a sampler backend | Observed in the generic-noise run (~22–46 evaluations); the cause was not investigated. Raising `--shots` or tightening `--tol` are untested options |
| Remote runs too slow | Reduce `--restarts`, `--A-values`, `--depths` and `--maxiter` |
| `RuntimeError: Job limit reached` | Expected on Qollab's remote IonQ simulator (one job per run). Use `--validate-on-backend` instead of an optimization loop |
| `No global backend found` | `--backend qollab` only works where the platform injects `backend`; elsewhere use `statevector`, `aer-noise` or `ionq` |
| "The stored angles do not reproduce the stored cost" | The records come from a different instance, `--n`/`--k` or normalization (for example the `--no-normalize` ablation); pass the matching options |

### 4. Limitations

- One synthetic instance with N = 6 and k = 3 studied. The N = 5, k = 2 variant (`--n 5 --k 2`, A_crit = 1.55) runs and passes the Hamiltonian check, but only as a smoke test; no sensitivity study was done for it. Conclusions about the optimum are specific to the instance's structure.
- COBYLA reaches the 400-evaluation cap in most exact runs, so p = 2 optima may be local; restarts and the cost-selected run mitigate but do not remove this.
- Medians over 8 restarts; no confidence intervals beyond the interquartile ranges.
- The normalization constant uses the exact spectrum, which is only available for small N.
- The ablation does not separate the effect of the energy scale from that of the initialization range.
- The generic noise model is illustrative; the IonQ path is untested.
- Only the standard X mixer is used; a constraint-preserving mixer is the natural next comparison and is not implemented.
- Qollab-specific paths (`--backend qollab`, `--validate-on-backend`) were tested only with a locally simulated injected backend and generic noise; no run on Qollab itself, and no IonQ-calibrated noise results, are included. One job per run on the remote IonQ simulator means noise cannot be included in the optimization loop; the noisy simulator is used only to validate stored optima.
