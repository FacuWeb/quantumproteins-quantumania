#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sensitivity of the Lagrange penalty parameter A in QAOA for protein-interactor
selection under a hard budget
===============================================================================

This effort is supported by Qollab & IonQ.

Qollab x IonQ Hackathon - Category: Quantum Impact (tag: quantum-impact).
License: MIT (see LICENSE).

Problem
-------
Choose exactly k of N candidate protein interactors (N = 5-6) maximizing

    value(x) = sum_i c_i x_i + sum_{i<j} J_ij x_i x_j        (x_i in {0, 1})

where c_i are individual confidence scores and J_ij > 0 / J_ij < 0 encode synergy /
redundancy between a pair. The hard budget sum_i x_i = k is enforced with a quadratic
(Lagrange-type) penalty, giving the QUBO

    H(x) = - sum_i c_i x_i - sum_{i<j} J_ij x_i x_j + A (sum_i x_i - k)^2.

Mapping to Ising / Pauli operators (x_i = (1 - Z_i)/2)
------------------------------------------------------
    l_i  = -c_i + A (1 - 2k)            (linear QUBO coefficient, using x_i^2 = x_i)
    q_ij = -J_ij + 2A                   (quadratic QUBO coefficient)
    H    = const + sum_i h_i Z_i + sum_{i<j} w_ij Z_i Z_j
    w_ij = q_ij / 4
    h_i  = -l_i / 2 - sum_{j != i} q_ij / 4
    const = A k^2 + sum_i l_i / 2 + sum_{i<j} q_ij / 4
The script checks this operator against brute-force energies of all 2^N states.

What the script does
--------------------
1. Builds a SYNTHETIC, illustrative interactor instance (not real experimental data).
2. Solves it exactly by enumeration (global optimum, greedy baseline, critical penalty
   A_crit above which the ground state of H is the constrained optimum).
3. For every A in a grid and QAOA depth p in {1, 2}: runs QAOA with COBYLA from several
   random restarts and measures the feasibility rate (shots with Hamming weight k),
   the approximation ratios and the probability of sampling the optimum.
4. Writes CSV/JSON/Markdown results and the figures.

Dependencies
------------
  numpy, scipy, matplotlib, qiskit (1.x), networkx
  optional: qiskit-aer (--backend aer-noise), qiskit-ionq < 1.0 (--backend ionq)

Quick use
---------
  python qaoa_penalty_sensitivity.py --quick
  python qaoa_penalty_sensitivity.py                       # reference configuration
  python qaoa_penalty_sensitivity.py --backend aer-noise   # generic noise (NOT IonQ-calibrated)
  # Qollab: optimize locally, then validate ONE stored optimum per run on the injected `backend`
  python qaoa_penalty_sensitivity.py --backend qollab --validate-on-backend --records qollab_validation_records.jsonl
  IONQ_API_KEY=... python qaoa_penalty_sensitivity.py --backend ionq --restarts 2 --maxiter 60

MIT License
Copyright (c) 2026 Facundo Ramírez
Permission is hereby granted, free of charge, to any person obtaining a copy of this
software and associated documentation files (the "Software"), to deal in the Software
without restriction, including without limitation the rights to use, copy, modify,
merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to the following
conditions: The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS
IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import os
import platform
import sys
import builtins
import time
from dataclasses import asdict, dataclass
from importlib import metadata
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import QAOAAnsatz
from qiskit.primitives import StatevectorEstimator, StatevectorSampler
from qiskit.quantum_info import SparsePauliOp, Statevector

BASIS_GATES = ["rz", "sx", "x", "cx"]
DEFAULT_A_GRID = [0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 7.5, 10.0, 15.0, 20.0]

# Used only when the script runs without command-line arguments AND a global `backend`
# has been injected by the platform (Qollab Code Playground). Edit to validate other A / p.
QOLLAB_DEFAULT_ARGS = ["--backend", "qollab", "--validate-on-backend",
                       "--records", "qollab_validation_records.jsonl",
                       "--validate-A", "7.5", "--validate-p", "2",
                       "--outdir", "results_qollab"]


# ----------------------------------------------------------------------------
# 1. Biomedical (interactomics) instance
# ----------------------------------------------------------------------------
@dataclass
class Instance:
    names: list
    c: np.ndarray      # individual confidence scores c_i
    J: np.ndarray      # symmetric pair matrix, zero diagonal; J>0 synergy, J<0 redundancy
    k: int

    @property
    def N(self) -> int:
        return len(self.names)


def default_instance(n: int = 6, k: int = 3) -> Instance:
    """Synthetic instance. P1 and P2 are individually best but redundant with each other;
    P3-P4-P5 form a synergistic complex, so picking the top-c proteins is a trap."""
    if n not in (5, 6):
        raise ValueError("This study uses N = 5 or 6 candidates")
    names = [f"P{i + 1}" for i in range(6)]
    c = np.array([0.90, 0.85, 0.80, 0.70, 0.65, 0.60])
    pairs = {(0, 1): -0.90,                       # redundancy (e.g. paralogs / same complex)
             (2, 3): 0.60, (3, 4): 0.50, (2, 4): 0.40,   # synergistic complex
             (0, 2): 0.10, (1, 5): 0.20, (0, 5): -0.30, (1, 3): 0.15}
    J = np.zeros((6, 6))
    for (i, j), v in pairs.items():
        J[i, j] = J[j, i] = v
    return Instance(names[:n], c[:n].copy(), J[:n, :n].copy(), k)


# ----------------------------------------------------------------------------
# 2. QUBO -> Ising/Pauli Hamiltonian and exact classical analysis
# ----------------------------------------------------------------------------
def all_states(N: int) -> np.ndarray:
    """Row s is the bit vector of integer s; bit i = qubit i (Qiskit little-endian)."""
    idx = np.arange(2 ** N)
    return ((idx[:, None] >> np.arange(N)) & 1).astype(float)


def bio_value(X: np.ndarray, inst: Instance) -> np.ndarray:
    return X @ inst.c + 0.5 * np.einsum("si,ij,sj->s", X, inst.J, X)


def qubo_energies(X: np.ndarray, inst: Instance, A: float) -> np.ndarray:
    return -bio_value(X, inst) + A * (X.sum(axis=1) - inst.k) ** 2


def pauli_hamiltonian(inst: Instance, A: float) -> SparsePauliOp:
    N, k = inst.N, inst.k
    l = -inst.c + A * (1 - 2 * k)
    q = -inst.J + 2 * A
    terms = []
    const = A * k ** 2 + l.sum() / 2
    for i in range(N):
        h = -l[i] / 2 - sum(q[i, j] / 4 for j in range(N) if j != i)
        terms.append(("Z", [i], h))
        for j in range(i + 1, N):
            terms.append(("ZZ", [i, j], q[i, j] / 4))
            const += q[i, j] / 4
    terms.append(("", [], const))
    return SparsePauliOp.from_sparse_list(terms, num_qubits=N).simplify()


def check_hamiltonian(inst: Instance, A: float) -> float:
    H = pauli_hamiltonian(inst, A)
    diag = np.real(np.diag(H.to_matrix()))
    err = float(np.max(np.abs(diag - qubo_energies(all_states(inst.N), inst, A))))
    if err > 1e-9:
        raise AssertionError(f"Pauli Hamiltonian disagrees with the QUBO (max error {err:.2e})")
    return err


def critical_penalty(inst: Instance) -> float:
    """Smallest A above which every infeasible state has higher energy than the optimum:
    A (s-k)^2 > value(x) - value(x*)  for every infeasible x."""
    X = all_states(inst.N)
    v, s = bio_value(X, inst), X.sum(axis=1)
    feas = s == inst.k
    v_opt = v[feas].max()
    inf = ~feas
    return float(max(0.0, np.max((v[inf] - v_opt) / (s[inf] - inst.k) ** 2)))


def classical_summary(inst: Instance) -> dict:
    X = all_states(inst.N)
    v, s = bio_value(X, inst), X.sum(axis=1)
    feas = np.where(s == inst.k)[0]
    best = feas[np.argmax(v[feas])]
    order = np.argsort(-inst.c)[: inst.k]
    greedy = int(sum(1 << int(i) for i in order))
    return {"opt_index": int(best), "opt_value": float(v[best]),
            "opt_set": [inst.names[i] for i in range(inst.N) if (best >> i) & 1],
            "worst_feasible_value": float(v[feas].min()),
            "mean_feasible_value": float(v[feas].mean()),
            "greedy_index": greedy, "greedy_value": float(v[greedy]),
            "greedy_set": [inst.names[i] for i in order],
            "n_feasible": int(len(feas)),
            "random_feasible_ARbio": float((v[feas].mean() - v[feas].min()) / (v[best] - v[feas].min())),
            "greedy_ARbio": float((v[greedy] - v[feas].min()) / (v[best] - v[feas].min())),
            "random_feasibility_pct": 100.0 * len(feas) / 2 ** inst.N,
            "A_crit": critical_penalty(inst)}


# ----------------------------------------------------------------------------
# 3. QAOA evaluators (exact simulation, noisy simulators, hardware-like backends)
# ----------------------------------------------------------------------------
class ExactEvaluator:
    """Cost from StatevectorEstimator (exact); sampling from StatevectorSampler (shots)."""

    def __init__(self, seed: int, shots: int):
        self.estimator = StatevectorEstimator(default_precision=0.0, seed=seed)
        self.sampler = StatevectorSampler(default_shots=shots, seed=seed)
        self.exact = True

    def bind(self, circ: QuantumCircuit, H_full: SparsePauliOp, energies: np.ndarray):
        self.circ, self.H = circ, H_full
        self.meas = circ.copy()
        self.meas.measure_all()
        return self

    def cost(self, theta) -> float:
        return float(np.real(self.estimator.run([(self.circ, self.H, np.asarray(theta))]).result()[0].data.evs))

    def sample_probs(self, theta, shots: int) -> np.ndarray:
        data = self.sampler.run([(self.meas, np.asarray(theta))], shots=shots).result()[0].data
        return _counts_to_probs(data.meas.get_counts(), self.circ.num_qubits)

    def exact_probs(self, theta) -> np.ndarray:
        return Statevector(self.circ.assign_parameters(np.asarray(theta))).probabilities()


class SamplerEvaluator:
    """Cost estimated from measured bitstrings (the Hamiltonian is diagonal)."""

    def __init__(self, sampler, shots: int, backend=None):
        self.sampler, self.shots, self.backend, self.exact = sampler, shots, backend, False

    def bind(self, circ: QuantumCircuit, H_full: SparsePauliOp, energies: np.ndarray):
        self.energies = energies
        circ = circ.copy()
        circ.measure_all()
        try:
            if self.backend is None:
                raise ValueError("no backend")
            self.meas = transpile(circ, backend=self.backend, optimization_level=1, seed_transpiler=7)
        except Exception as exc:  # backend without a usable target: fall back to the common basis
            if self.backend is not None:
                print(f"  [warning] transpiling for the backend failed ({exc}); using {BASIS_GATES}")
            self.meas = transpile(circ, basis_gates=BASIS_GATES + ["measure"], optimization_level=1,
                                  seed_transpiler=7)
        self.n = circ.num_qubits - 0
        return self

    def _probs(self, theta, shots):
        res = self.sampler.run([(self.meas, np.asarray(theta))], shots=shots).result()[0]
        return _counts_to_probs(res.data.meas.get_counts(), self.n)

    def cost(self, theta) -> float:
        return float(self._probs(theta, self.shots) @ self.energies)

    def sample_probs(self, theta, shots: int) -> np.ndarray:
        return self._probs(theta, shots)

    def exact_probs(self, theta):
        return None


def _counts_to_probs(counts: dict, n: int) -> np.ndarray:
    p = np.zeros(2 ** n)
    for key, v in counts.items():
        p[int(key.replace(" ", ""), 2)] += v
    return p / p.sum()


def get_injected_backend():
    """Return the global `backend` that Qollab injects into the execution namespace, or None.
    Looks in this module, __main__, builtins and the globals of every calling frame."""
    import __main__
    spaces = [globals(), vars(__main__), vars(builtins)]
    frame = sys._getframe(1)
    while frame is not None:
        spaces.append(frame.f_globals)
        frame = frame.f_back
    for ns in spaces:
        b = ns.get("backend")
        if b is not None and hasattr(b, "run"):
            return b
    return None


def default_max_jobs(backend):
    """Qollab's remote IonQ simulators accept one job per run; other injected backends are not
    limited by default. Heuristic based on the backend name (e.g. 'ionq_simulator')."""
    name = getattr(backend, "name", "")
    name = name() if callable(name) else name
    return 1 if "ionq" in str(name).lower() else None


class JobGuard:
    """Refuses to submit more than `max_jobs` primitive calls (Qollab allows 1 job per run)."""

    def __init__(self, primitive, max_jobs):
        self._primitive, self.max_jobs, self.count = primitive, max_jobs, 0

    def run(self, *a, **k):
        if self.max_jobs is not None and self.count >= self.max_jobs:
            raise RuntimeError(
                f"Job limit reached ({self.max_jobs}). The Qollab simulators accept one job per run, so "
                "optimization loops cannot run on that backend. Optimize locally and validate the stored "
                "parameters with --validate-on-backend (one job per run).")
        self.count += 1
        return self._primitive.run(*a, **k)


def make_evaluator(args):
    if args.backend == "statevector":
        return ExactEvaluator(args.seed, args.shots), "StatevectorEstimator (cost, exact) + StatevectorSampler (shots)"
    if args.backend == "aer-noise":
        from qiskit_aer.noise import NoiseModel, depolarizing_error
        from qiskit_aer.primitives import SamplerV2 as AerSampler
        nm = NoiseModel()
        nm.add_all_qubit_quantum_error(depolarizing_error(args.p1, 1), ["sx", "x"])
        nm.add_all_qubit_quantum_error(depolarizing_error(args.p2, 2), ["cx"])
        sampler = AerSampler(options={"backend_options": {"noise_model": nm, "seed_simulator": args.seed}})
        desc = (f"Aer SamplerV2 with GENERIC depolarizing noise (1q error {args.p1}, 2q error {args.p2}); "
                f"NOT the IonQ noise model; {args.shots} shots per evaluation")
        return SamplerEvaluator(sampler, args.shots), desc
    if args.backend == "ionq":
        from qiskit.primitives import BackendSamplerV2
        from qiskit_ionq import IonQProvider
        token = os.environ.get("IONQ_API_KEY")
        if not token:
            raise SystemExit("Set the IONQ_API_KEY environment variable to use --backend ionq")
        backend = IonQProvider(token).get_backend("simulator")
        backend.set_options(noise_model=args.ionq_noise_model)
        sampler = BackendSamplerV2(backend=backend, options={"default_shots": args.shots})
        return SamplerEvaluator(sampler, args.shots, backend), \
            f"IonQ cloud simulator, noise_model={args.ionq_noise_model} (UNTESTED by the script author)"
    if args.backend == "qollab":
        from qiskit.primitives import BackendSamplerV2
        backend = get_injected_backend()
        if backend is None:
            raise SystemExit("No global `backend` found. Run inside Qollab (choose the simulator in "
                             "'Run on Qollab / Select QPU') or use --backend statevector / aer-noise / ionq.")
        if args.qollab_noise_model and hasattr(backend, "set_options"):
            backend.set_options(noise_model=args.qollab_noise_model)
        name = getattr(backend, "name", "unknown")
        name = name() if callable(name) else name
        max_jobs = default_max_jobs(backend) if args.max_jobs is None else (args.max_jobs if args.max_jobs > 0 else None)
        sampler = JobGuard(BackendSamplerV2(backend=backend, options={"default_shots": args.shots}), max_jobs)
        return SamplerEvaluator(sampler, args.shots, backend), \
            f"Qollab injected backend '{name}' (job limit per run: {max_jobs if max_jobs else 'none'})"
    if args.backend == "custom":
        if not args.sampler_factory or ":" not in args.sampler_factory:
            raise SystemExit("--backend custom requires --sampler-factory file.py:function")
        path, func = args.sampler_factory.rsplit(":", 1)
        spec = importlib.util.spec_from_file_location("custom_sampler_module", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        out = getattr(module, func)(args)
        sampler, backend = out if isinstance(out, tuple) else (out, None)
        return SamplerEvaluator(sampler, args.shots, backend), f"Custom sampler factory {args.sampler_factory}"
    raise ValueError(args.backend)


# ----------------------------------------------------------------------------
# 4. QAOA runs and metrics
# ----------------------------------------------------------------------------
@dataclass
class Record:
    backend: str
    A: float
    p: int
    restart: int
    nfev: int
    time_s: float
    cost_final: float        # <H> / scale at the optimizer's end (cost used by COBYLA)
    E_H: float               # <H> in original units, from the final shots
    AR_H: float              # (E_max - <H>) / (E_max - E_min)
    feasibility_pct: float   # % of shots with Hamming weight k
    AR_bio: float            # feasible-conditioned, normalized objective ratio
    P_opt: float             # probability of sampling the global optimum
    feas_exact_pct: float    # exact probabilities (nan if unavailable)
    P_opt_exact: float
    probs: list
    theta: list


def build_circuit(H_evo: SparsePauliOp, p: int) -> QuantumCircuit:
    ansatz = QAOAAnsatz(cost_operator=H_evo, reps=p)
    return transpile(ansatz, basis_gates=BASIS_GATES, optimization_level=1, seed_transpiler=7)


def init_theta(circ, rng) -> np.ndarray:
    th = []
    for prm in circ.parameters:
        th.append(rng.uniform(0, math.pi) if prm.name.startswith("β") else rng.uniform(0, 2 * math.pi))
    return np.array(th)


def metrics_from_probs(pr, X, v, s, inst, E_orig, cls) -> dict:
    feas = s == inst.k
    pf = float(pr[feas].sum())
    E_min, E_max = float(E_orig.min()), float(E_orig.max())
    E_H = float(pr @ E_orig)
    if pf > 0:
        ev_feas = float((pr[feas] * v[feas]).sum() / pf)
        ar_bio = (ev_feas - cls["worst_feasible_value"]) / (cls["opt_value"] - cls["worst_feasible_value"])
    else:
        ar_bio = 0.0
    return {"E_H": E_H, "AR_H": (E_max - E_H) / (E_max - E_min), "feas": 100 * pf,
            "AR_bio": float(ar_bio), "P_opt": float(pr[cls["opt_index"]])}


def run_one(ev, circ, scale, E_orig, X, v, s, inst, cls, A, p, restart, args, backend_name):
    rng = np.random.default_rng([args.seed, int(round(A * 1000)), p, restart])
    theta0 = init_theta(circ, rng)
    n_eval = [0]

    def fun(th):
        n_eval[0] += 1
        return ev.cost(th)

    t0 = time.perf_counter()
    res = minimize(fun, theta0, method="COBYLA",
                   options={"maxiter": args.maxiter, "rhobeg": args.rhobeg, "tol": args.tol})
    dt = time.perf_counter() - t0
    pr = ev.sample_probs(res.x, args.final_shots)
    m = metrics_from_probs(pr, X, v, s, inst, E_orig, cls)
    pe = ev.exact_probs(res.x)
    if pe is not None:
        me = metrics_from_probs(pe, X, v, s, inst, E_orig, cls)
        fe, poe = me["feas"], me["P_opt"]
    else:
        fe, poe = float("nan"), float("nan")
    return Record(backend=backend_name, A=A, p=p, restart=restart, nfev=n_eval[0], time_s=dt,
                  cost_final=float(res.fun), E_H=m["E_H"], AR_H=m["AR_H"], feasibility_pct=m["feas"],
                  AR_bio=m["AR_bio"], P_opt=m["P_opt"], feas_exact_pct=fe, P_opt_exact=poe,
                  probs=[float(x) for x in pr], theta=[float(x) for x in res.x])


# ----------------------------------------------------------------------------
# 5. Aggregation, tables and figures
# ----------------------------------------------------------------------------
def aggregate(records, A_values, depths):
    rows = []
    for p in depths:
        for A in A_values:
            rs = [r for r in records if r.p == p and r.A == A]
            if not rs:
                continue
            sel = min(rs, key=lambda r: r.cost_final)    # restart chosen by cost (no oracle needed)
            q = lambda a, f: float(np.percentile([getattr(r, a) for r in rs], f))
            rows.append({
                "p": p, "A": A, "n_restarts": len(rs),
                "feas_med": q("feasibility_pct", 50), "feas_q25": q("feasibility_pct", 25), "feas_q75": q("feasibility_pct", 75),
                "ARbio_med": q("AR_bio", 50), "ARbio_q25": q("AR_bio", 25), "ARbio_q75": q("AR_bio", 75),
                "ARH_med": q("AR_H", 50),
                "ARo_med": float(np.median([r.feasibility_pct / 100 * r.AR_bio for r in rs])),
                "Popt_med": q("P_opt", 50), "Popt_q25": q("P_opt", 25), "Popt_q75": q("P_opt", 75),
                "nfev_med": q("nfev", 50),
                "shots_to_opt": (1.0 / q("P_opt", 50)) if q("P_opt", 50) > 0 else float("inf"),
                "sel_feas": sel.feasibility_pct, "sel_ARbio": sel.AR_bio, "sel_Popt": sel.P_opt,
                "sel_restart": sel.restart})
    return rows


def summary_markdown(rows, cls, norm_gap) -> str:
    out = ("| p | A | A / A_crit | norm. gap | feasibility % (median [IQR]) | AR_bio (median) | AR_overall (median) | "
           "P_opt % (median) | expected shots to optimum | cost-selected run: feas % / AR_bio / P_opt % | median nfev |\n"
           "|---|---|---|---|---|---|---|---|---|---|---|\n")
    for r in rows:
        out += (f"| {r['p']} | {r['A']:g} | {r['A'] / cls['A_crit']:.2f} | {norm_gap[r['A']]:.3f} | "
                f"{r['feas_med']:.1f} [{r['feas_q25']:.1f}, {r['feas_q75']:.1f}] | {r['ARbio_med']:.3f} | {r['ARo_med']:.3f} | "
                f"{100 * r['Popt_med']:.1f} | {r['shots_to_opt']:.0f} | "
                f"{r['sel_feas']:.1f} / {r['sel_ARbio']:.3f} / {100 * r['sel_Popt']:.1f} | {r['nfev_med']:.0f} |\n")
    return out


def figure_network(inst, cls, outdir: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(inst.names)
    for i in range(inst.N):
        for j in range(i + 1, inst.N):
            if inst.J[i, j] != 0:
                G.add_edge(inst.names[i], inst.names[j], J=inst.J[i, j])
    pos = nx.circular_layout(G)
    opt = set(cls["opt_set"])
    greedy = set(cls["greedy_set"])
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    colors = ["#f28e2b" if n in opt else "#a6cee3" for n in G.nodes]
    edge_col = ["#2ca02c" if d["J"] > 0 else "#d62728" for _, _, d in G.edges(data=True)]
    widths = [1 + 5 * abs(d["J"]) for _, _, d in G.edges(data=True)]
    styles = ["solid" if d["J"] > 0 else "dashed" for _, _, d in G.edges(data=True)]
    nx.draw_networkx_edges(G, pos, edge_color=edge_col, width=widths, style=styles, ax=ax, alpha=0.85)
    nx.draw_networkx_nodes(G, pos, node_color=colors, edgecolors=["black" if n in greedy else "none" for n in G.nodes],
                           linewidths=2.5, node_size=[1400 * inst.c[inst.names.index(n)] for n in G.nodes], ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=11, ax=ax)
    nx.draw_networkx_edge_labels(G, pos, edge_labels={(u, w): f"{d['J']:+.2f}" for u, w, d in G.edges(data=True)},
                                 font_size=8, ax=ax)
    for n, (x, y) in pos.items():
        ax.text(x * 1.22, y * 1.22, f"c={inst.c[inst.names.index(n)]:.2f}", ha="center", va="center", fontsize=8)
    ax.set_title("Figure 1. Interactor network (synthetic instance)\n"
                 f"Orange: optimal subset {sorted(opt)} (value {cls['opt_value']:.2f}) | "
                 f"black outline: greedy-by-c subset (value {cls['greedy_value']:.2f})\n"
                 "green solid = synergy (J>0), red dashed = redundancy (J<0), node size ~ c", fontsize=8.5)
    ax.axis("off")
    ax.set_xlim(-1.4, 1.4); ax.set_ylim(-1.4, 1.4)
    fig.tight_layout(); fig.savefig(outdir / "fig1_network.png", dpi=170); plt.close(fig)


def figure_feasibility_ratio(rows, cls, depths, outdir: Path, backend_label: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, len(depths), figsize=(6.4 * len(depths), 4.8), squeeze=False)
    for ax, p in zip(axes[0], depths):
        rs = [r for r in rows if r["p"] == p]
        A = np.array([r["A"] for r in rs])
        ax.plot(A, [r["feas_med"] for r in rs], "o-", color="tab:blue", label="Feasibility rate (%) - median")
        ax.fill_between(A, [r["feas_q25"] for r in rs], [r["feas_q75"] for r in rs], color="tab:blue", alpha=0.15)
        ax.axhline(cls["random_feasibility_pct"], color="tab:blue", ls=":", lw=1,
                   label=f"Uniform-superposition baseline ({cls['random_feasibility_pct']:.1f} %)")
        ax.set_xscale("log"); ax.set_xlabel("Penalty parameter A"); ax.set_ylabel("Feasibility rate (%)", color="tab:blue")
        ax.set_ylim(0, 105); ax.grid(alpha=0.3, which="both")
        ax2 = ax.twinx()
        ax2.plot(A, [r["ARbio_med"] for r in rs], "s-", color="tab:red", label="Approx. ratio AR_bio - median")
        ax2.fill_between(A, [r["ARbio_q25"] for r in rs], [r["ARbio_q75"] for r in rs], color="tab:red", alpha=0.12)
        ax2.plot(A, [r["ARo_med"] for r in rs], "D-", color="tab:orange", label="AR_overall = feasibility x AR_bio - median")
        ax2.plot(A, [r["Popt_med"] for r in rs], "^--", color="tab:green", label="P(optimum) - median")
        ax2.axhline(cls["random_feasible_ARbio"], color="tab:red", ls=":", lw=1,
                    label=f"AR_bio of a random feasible pick ({cls['random_feasible_ARbio']:.2f})")
        ax2.set_ylim(0, 1.05); ax2.set_ylabel("Ratio / probability", color="tab:red")
        ax.axvline(cls["A_crit"], color="k", ls="--", lw=1)
        ax.text(cls["A_crit"], 4, " A_crit", fontsize=8, rotation=90, va="bottom")
        h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=6.5, loc="lower right")
        ax.set_title(f"QAOA p = {p}")
    fig.suptitle(f"Figure 2. Feasibility rate and approximation ratio vs penalty A ({backend_label})", fontsize=10)
    fig.tight_layout(); fig.savefig(outdir / "fig2_feasibility_ratio_vs_A.png", dpi=170); plt.close(fig)


def figure_histogram(records, rows, inst, cls, p, outdir: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rs = [r for r in rows if r["p"] == p]
    A_best = max(rs, key=lambda r: r["ARo_med"])["A"]      # highest median AR_overall
    A_sub = min(rs, key=lambda r: r["ARo_med"])["A"]       # lowest median AR_overall
    N, k = inst.N, inst.k
    pop = [bin(i).count("1") for i in range(2 ** N)]
    order = sorted(range(2 ** N), key=lambda i: (abs(pop[i] - k), i))
    fig, axes = plt.subplots(2, 1, figsize=(13, 7.6), sharey=True)
    for ax, A, tag in zip(axes, (A_sub, A_best), ("least favourable", "most favourable")):
        rr = [r for r in records if r.p == p and r.A == A]
        sel = min(rr, key=lambda r: r.cost_final)
        pr = np.array(sel.probs)
        cols = ["#f28e2b" if i == cls["opt_index"] else ("#1f77b4" if pop[i] == k else "#b0b0b0") for i in order]
        ax.bar(range(2 ** N), pr[order], color=cols)
        ax.set_xticks(range(2 ** N))
        ax.set_xticklabels([format(i, f"0{N}b") for i in order], rotation=90, fontsize=5.5)
        ax.set_ylabel("Measurement probability")
        ax.set_title(f"A = {A:g} ({tag} A at p = {p} by median AR_overall; A_crit = {cls['A_crit']:.2f}) - cost-selected restart: "
                     f"feasibility {sel.feasibility_pct:.1f} %, P(optimum) = {100 * sel.P_opt:.1f} %", fontsize=9)
        ax.grid(alpha=0.3, axis="y")
    from matplotlib.patches import Patch
    axes[0].legend(handles=[Patch(color="#f28e2b", label="global optimum"),
                            Patch(color="#1f77b4", label=f"feasible (Hamming weight = {k})"),
                            Patch(color="#b0b0b0", label="infeasible")], fontsize=8)
    axes[1].set_xlabel("Bitstring (q_{N-1} ... q_0), feasible states first")
    fig.suptitle("Figure 3. Measurement histograms: least vs most favourable penalty A", fontsize=10)
    fig.tight_layout(); fig.savefig(outdir / "fig3_histogram.png", dpi=170); plt.close(fig)
    return A_sub, A_best


def figure_diagnostics(A_values, norm_gap, cls, outdir: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4.2))
    ax.plot(A_values, [norm_gap[a] for a in A_values], "o-")
    ax.axvline(cls["A_crit"], color="k", ls="--", lw=1, label=f"A_crit = {cls['A_crit']:.2f}")
    ax.set_xscale("log"); ax.set_xlabel("Penalty parameter A")
    ax.set_ylabel("Normalized gap  (E_2nd - E_0) / (E_max - E_min)")
    ax.set_title("Figure 4 (supplementary). Signal dilution: normalized gap vs A\n(negative gap = ground state of H is not the constrained optimum)", fontsize=8.5)
    ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(outdir / "fig4_gap_vs_A.png", dpi=170); plt.close(fig)


def normalized_gap(inst, cls, A) -> float:
    """(E of best non-optimal state - E of optimum) / (E_max - E_min); negative if the
    constrained optimum is not the ground state of H."""
    X = all_states(inst.N)
    E = qubo_energies(X, inst, A)
    opt = cls["opt_index"]
    others = np.delete(E, opt)
    return float((others.min() - E[opt]) / (E.max() - E.min()))


def figure_validation(pr_ref, pr_bk, inst, cls, A, p, backend_desc, outdir: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    N, k = inst.N, inst.k
    pop = [bin(i).count("1") for i in range(2 ** N)]
    order = sorted(range(2 ** N), key=lambda i: (abs(pop[i] - k), i))
    fig, ax = plt.subplots(figsize=(13, 4.4))
    x = np.arange(2 ** N)
    ax.bar(x - 0.2, np.asarray(pr_ref)[order], width=0.4, color="#1f77b4", label="noiseless (exact)")
    ax.bar(x + 0.2, np.asarray(pr_bk)[order], width=0.4, color="#ff7f0e", label="backend (shots)")
    opt_pos = order.index(cls["opt_index"])
    ax.axvline(opt_pos, color="k", ls=":", lw=1)
    ax.text(opt_pos, ax.get_ylim()[1] * 0.95, " optimum", fontsize=8, va="top")
    ax.set_xticks(x)
    ax.set_xticklabels([format(i, f"0{N}b") for i in order], rotation=90, fontsize=5.5)
    ax.set_ylabel("Measurement probability"); ax.grid(alpha=0.3, axis="y"); ax.legend(fontsize=8)
    ax.set_title(f"Figure 5. Validation of stored optimal angles (A = {A:g}, p = {p}); feasible states first\n"
                 f"{backend_desc}", fontsize=8.5)
    fig.tight_layout()
    path = outdir / f"fig5_backend_validation_A{A:g}_p{p}.png"
    fig.savefig(path, dpi=170); plt.close(fig)
    return path


def validate_on_backend(args):
    """Run the stored optimal angles of the selected restarts on the chosen backend.
    One sampler call (one job) per value of A. No optimization is performed."""
    args.outdir.mkdir(parents=True, exist_ok=True)
    inst = default_instance(args.n, args.k)
    cls = classical_summary(inst)
    X = all_states(inst.N)
    v, s = bio_value(X, inst), X.sum(axis=1)
    rec_path = args.records if args.records is not None else args.outdir / "records.jsonl"
    records = [Record(**json.loads(l)) for l in Path(rec_path).read_text().splitlines() if l.strip()]
    ev, desc = make_evaluator(args)
    print("Backend:", desc)
    out = []
    for A in args.validate_A:
        rr = [r for r in records if r.p == args.validate_p and abs(r.A - A) < 1e-9]
        if not rr:
            raise SystemExit(f"No records for A={A}, p={args.validate_p} in {rec_path}")
        sel = min(rr, key=lambda r: r.cost_final)           # restart selected by cost
        E_orig = qubo_energies(X, inst, A)
        scale = 1.0 if args.no_normalize else float(E_orig.max() - E_orig.min())
        H_full = pauli_hamiltonian(inst, A) / scale if scale != 1.0 else pauli_hamiltonian(inst, A)
        H_evo = SparsePauliOp.from_list([(l, c) for l, c in zip(H_full.paulis.to_labels(), H_full.coeffs)
                                         if set(l) != {"I"}])
        circ = build_circuit(H_evo, args.validate_p)
        theta = np.array(sel.theta)
        sv = Statevector(circ.assign_parameters(theta))
        cost_check = float(np.real(sv.expectation_value(H_full)))
        if abs(cost_check - sel.cost_final) > 1e-6:
            raise SystemExit("The stored angles do not reproduce the stored cost for this instance and "
                             f"normalization (difference {abs(cost_check - sel.cost_final):.2e}); check "
                             "--no-normalize, --n and --k against the run that produced the records.")
        pr_ref = sv.probabilities()
        ev.bind(circ, H_full, E_orig / scale)
        pr_bk = ev.sample_probs(theta, args.final_shots)    # ONE sampler call = one job
        m_ref = metrics_from_probs(pr_ref, X, v, s, inst, E_orig, cls)
        m_bk = metrics_from_probs(pr_bk, X, v, s, inst, E_orig, cls)
        print(f"\nValidation: A = {A:g} (A/A_crit = {A / cls['A_crit']:.2f}), p = {args.validate_p}, "
              f"restart {sel.restart}, {args.final_shots} shots")
        print("| Metric | Noiseless (exact) | Backend (shots) |\n|---|---|---|")
        for label, key, f in (("Feasibility %", "feas", "{:.1f}"), ("AR_bio", "AR_bio", "{:.3f}"),
                              ("P_opt %", "P_opt", "{:.1f}"), ("AR_H", "AR_H", "{:.3f}")):
            sc = 100 if key == "P_opt" else 1
            print(f"| {label} | {f.format(sc * m_ref[key])} | {f.format(sc * m_bk[key])} |")
        fig = figure_validation(pr_ref, pr_bk, inst, cls, A, args.validate_p, desc, args.outdir)
        out.append({"A": A, "p": args.validate_p, "restart": sel.restart, "shots": args.final_shots,
                    "noiseless": m_ref, "backend": m_bk, "probs_backend": [float(x) for x in pr_bk],
                    "figure": str(fig)})
    path = args.outdir / "validation_results.json"
    path.write_text(json.dumps({"backend": desc, "results": out}, indent=2))
    print(f"\nSaved {path} and the validation figure(s) in {args.outdir.resolve()}")


# ----------------------------------------------------------------------------
# 6. Command line
# ----------------------------------------------------------------------------
def parse_args(argv=None):
    ap = argparse.ArgumentParser(prog="qaoa_penalty_sensitivity.py",   # explicit: sys.argv may be empty in a playground
                                 description="Lagrange-penalty sensitivity of QAOA for budgeted interactor selection")
    ap.add_argument("--n", type=int, default=6, choices=[5, 6], help="number of candidate proteins")
    ap.add_argument("--k", type=int, default=3, choices=[2, 3], help="budget: proteins to select")
    ap.add_argument("--A-values", type=float, nargs="+", default=DEFAULT_A_GRID)
    ap.add_argument("--depths", type=int, nargs="+", default=[1, 2])
    ap.add_argument("--restarts", type=int, default=8, help="random COBYLA restarts per (A, p)")
    ap.add_argument("--no-normalize", action="store_true",
                    help="do not divide H by its spectral range (exposes the energy-scale effect of A)")
    ap.add_argument("--maxiter", type=int, default=400)
    ap.add_argument("--rhobeg", type=float, default=0.5)
    ap.add_argument("--tol", type=float, default=1e-3)
    ap.add_argument("--shots", type=int, default=4096, help="shots per cost evaluation (sampler backends)")
    ap.add_argument("--final-shots", type=int, default=8192, help="shots for the final distribution")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--backend", choices=["statevector", "aer-noise", "ionq", "qollab", "custom"], default="statevector")
    ap.add_argument("--p1", type=float, default=2e-4, help="1q depolarizing error (aer-noise, illustrative)")
    ap.add_argument("--p2", type=float, default=5e-3, help="2q depolarizing error (aer-noise, illustrative)")
    ap.add_argument("--ionq-noise-model", default="forte-1", help="check IonQ docs for current model names")
    ap.add_argument("--sampler-factory", default=None, help="file.py:function (with --backend custom)")
    ap.add_argument("--outdir", type=Path, default=Path("results"))
    ap.add_argument("--quick", action="store_true", help="short test: 3 values of A, 2 restarts")
    ap.add_argument("--replot", action="store_true", help="regenerate tables and figures from records.jsonl")
    ap.add_argument("--qollab-noise-model", default=None,
                    help="optional noise model set on the injected backend, e.g. forte-1 (leave unset to keep "
                         "the selection made in the Qollab run menu)")
    ap.add_argument("--max-jobs", type=int, default=None,
                    help="maximum sampler calls; default: 1 for an injected IonQ backend (Qollab: one job per run), "
                         "unlimited otherwise; 0 = unlimited")
    ap.add_argument("--validate-on-backend", action="store_true",
                    help="skip optimization: run stored optimal angles on the chosen backend")
    ap.add_argument("--records", type=Path, default=None, help="records file with optimal angles (validation)")
    ap.add_argument("--validate-A", type=float, nargs="+", default=[7.5], help="penalty values to validate")
    ap.add_argument("--validate-p", type=int, default=2, help="QAOA depth to validate")
    if argv is None and len(sys.argv) <= 1 and get_injected_backend() is not None:
        argv = QOLLAB_DEFAULT_ARGS          # running inside Qollab without command-line arguments
        print("Qollab backend detected; using QOLLAB_DEFAULT_ARGS:", " ".join(argv))
    args = ap.parse_args(argv)
    if args.quick:
        args.A_values, args.restarts = [0.5, 2.0, 10.0], 2
    return args


def _versions():
    out = {}
    for pkg in ["qiskit", "qiskit-aer", "qiskit-ionq", "numpy", "scipy", "matplotlib", "networkx"]:
        try:
            out[pkg] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            out[pkg] = None
    out["python"] = platform.python_version()
    return out


def main(argv=None):
    args = parse_args(argv)
    if args.validate_on_backend:
        validate_on_backend(args)
        return
    args.outdir.mkdir(parents=True, exist_ok=True)
    inst = default_instance(args.n, args.k)
    N = inst.N
    cls = classical_summary(inst)
    X = all_states(N)
    v, s = bio_value(X, inst), X.sum(axis=1)
    for A in args.A_values:
        check_hamiltonian(inst, A)
    print(f"Instance N={N}, k={inst.k}: optimum {cls['opt_set']} value {cls['opt_value']:.3f} | "
          f"greedy-by-c {cls['greedy_set']} value {cls['greedy_value']:.3f} "
          f"(ratio {cls['greedy_value'] / cls['opt_value']:.2f})")
    print(f"A_crit = {cls['A_crit']:.3f} | uniform-superposition feasibility = {cls['random_feasibility_pct']:.2f} % "
          f"| Pauli Hamiltonian verified against all {2 ** N} QUBO energies")
    print(f"Baselines: AR_bio of greedy-by-c = {cls['greedy_ARbio']:.3f} | AR_bio of a random feasible pick = "
          f"{cls['random_feasible_ARbio']:.3f} | uniform P(optimum) = {100 / 2 ** N:.2f} %")

    norm_gap = {A: normalized_gap(inst, cls, A) for A in args.A_values}
    records, t_start = [], time.perf_counter()
    rec_path = args.outdir / "records.jsonl"
    if args.replot:
        records = [Record(**json.loads(l)) for l in rec_path.read_text().splitlines() if l.strip()]
        backend_desc = "(regenerated from records.jsonl)"
        loop_depths = []
    else:
        ev, backend_desc = make_evaluator(args)
        print("Backend:", backend_desc)
        loop_depths = args.depths
    for p in loop_depths:
        for A in args.A_values:
            E_orig = qubo_energies(X, inst, A)
            scale = 1.0 if args.no_normalize else float(E_orig.max() - E_orig.min())
            H_full = pauli_hamiltonian(inst, A) / scale if scale != 1.0 else pauli_hamiltonian(inst, A)
            H_evo = SparsePauliOp.from_list([(l, c) for l, c in zip(H_full.paulis.to_labels(), H_full.coeffs)
                                             if set(l) != {"I"}])
            circ = build_circuit(H_evo, p)
            ev.bind(circ, H_full, E_orig / scale)
            for r in range(args.restarts):
                records.append(run_one(ev, circ, scale, E_orig, X, v, s, inst, cls, A, p, r, args, args.backend))
            last = [x for x in records if x.p == p and x.A == A]
            print(f"  p={p} A={A:5g} | feasibility med {np.median([x.feasibility_pct for x in last]):5.1f} % | "
                  f"AR_bio med {np.median([x.AR_bio for x in last]):.3f} | P_opt med "
                  f"{100 * np.median([x.P_opt for x in last]):5.1f} % | nfev med {np.median([x.nfev for x in last]):.0f}",
                  flush=True)
    wall = time.perf_counter() - t_start
    if not args.replot:
        with open(rec_path, "w") as fh:
            for r in records:
                fh.write(json.dumps(asdict(r)) + "\n")

    rows = aggregate(records, args.A_values, args.depths)
    md = summary_markdown(rows, cls, norm_gap)
    (args.outdir / "summary.md").write_text(md)
    with open(args.outdir / "metrics.csv", "w", newline="") as fh:
        fields = [f for f in asdict(records[0]) if f not in ("probs", "theta")]
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in records:
            d = asdict(r); d.pop("probs"); d.pop("theta")
            w.writerow(d)
    figure_network(inst, cls, args.outdir)
    label = "exact statevector" if args.backend == "statevector" else args.backend
    figure_feasibility_ratio(rows, cls, args.depths, args.outdir, label)
    hist_p = max(args.depths)
    A_sub, A_best = figure_histogram(records, rows, inst, cls, hist_p, args.outdir)
    figure_diagnostics(args.A_values, norm_gap, cls, args.outdir)
    meta = {"args": {k_: str(v_) for k_, v_ in vars(args).items()}, "backend": backend_desc,
            "versions": _versions(), "wall_time_s": wall, "classical": cls,
            "instance": {"names": inst.names, "c": inst.c.tolist(), "J": inst.J.tolist(), "k": inst.k},
            "normalized_gap": {str(a): g for a, g in norm_gap.items()},
            "histogram_A": {"least_favourable": A_sub, "most_favourable": A_best, "p": hist_p},
            "summary_rows": rows}
    (args.outdir / "run_metadata.json").write_text(json.dumps(meta, indent=2))
    print("\n" + md)
    print(f"Histogram uses A_least = {A_sub:g}, A_best = {A_best:g} (p = {hist_p})")
    print(f"Results written to {args.outdir.resolve()} (wall time {wall:.0f} s)")


if __name__ == "__main__":
    main()
