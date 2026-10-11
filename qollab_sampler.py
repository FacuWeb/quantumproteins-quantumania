#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
qollab_sampler.py - plug-in for `--backend custom` (QAOA penalty-sensitivity project)

This effort is supported by Qollab & IonQ.
License: MIT (see LICENSE).

Usage
-----
    python main.py.py --backend custom \
        --sampler-factory qollab_sampler.py:make_sampler

The script calls make_sampler(args) and accepts either `sampler` or
`(sampler, backend)`:
  * `sampler` must follow the Qiskit primitives V2 interface:
        sampler.run([(circuit, parameter_values)], shots=N).result()[0].data.meas.get_counts()
    The circuits sent to it end with `measure_all()`, so the classical register is
    named "meas" and bitstrings use Qiskit's little-endian order (qubit 0 rightmost).
  * If a `backend` is returned, the circuits are transpiled for it. If not, they are
    transpiled to {rz, sx, x, cx} plus measurement only.

Mode selection (environment variable QOLLAB_MODE)
-------------------------------------------------
    local           (default) seeded statevector sampling. TESTED.
    ionq-simulator  IonQ cloud simulator through qiskit-ionq (< 1.0 for Qiskit 1.x).
                    Needs IONQ_API_KEY. Written from IonQ's documentation. UNTESTED.
    qollab          the global `backend` injected by Qollab (Code Playground). No API key.
                    One job per run: use it with --validate-on-backend, not for optimization
                    loops (a guard stops after QOLLAB_MAX_JOBS primitive calls, default 1).
                    Written from the platform information provided by the author. TESTED only
                    with a locally simulated injected backend, not on Qollab itself.

Optional environment variables
------------------------------
    QOLLAB_NOISE_MODEL optional noise model set on the injected backend (e.g. forte-1); if unset,
                       the choice made in the Qollab run menu is kept
    QOLLAB_MAX_JOBS    primitive calls allowed in qollab mode (default: 1 for an IonQ backend,
                       unlimited otherwise; 0 = unlimited)
    IONQ_NOISE_MODEL   noise model name for the IonQ simulator (default "forte-1";
                       check IonQ's documentation for the models currently available)

Remember that every COBYLA evaluation on a sampler backend is one job: for remote
runs reduce --restarts, --A-values, --depths and --maxiter first.
"""
import builtins
import os
import sys

IONQ_NOISE_MODEL = os.environ.get("IONQ_NOISE_MODEL", "forte-1")



def _find_injected_backend():
    """Global `backend` injected by Qollab (Code Playground). Looks in __main__, builtins and
    the globals of every calling frame."""
    import __main__
    spaces = [vars(__main__), vars(builtins)]
    frame = sys._getframe(1)
    while frame is not None:
        spaces.append(frame.f_globals)
        frame = frame.f_back
    for ns in spaces:
        b = ns.get("backend")
        if b is not None and hasattr(b, "run"):
            return b
    return None


class _JobGuard:
    """Stops after `max_jobs` primitive calls (Qollab accepts one job per run)."""

    def __init__(self, primitive, max_jobs):
        self._primitive, self.max_jobs, self.count = primitive, max_jobs, 0

    def run(self, *a, **k):
        if self.max_jobs is not None and self.count >= self.max_jobs:
            raise RuntimeError("Job limit reached. Qollab simulators accept one job per run: optimize "
                               "locally and validate stored parameters with --validate-on-backend.")
        self.count += 1
        return self._primitive.run(*a, **k)


def _max_jobs(backend):
    """QOLLAB_MAX_JOBS if set (0 = unlimited); otherwise 1 for an IonQ backend (one job per run
    on Qollab's remote IonQ simulators) and unlimited for any other injected backend."""
    env = os.environ.get("QOLLAB_MAX_JOBS")
    if env is not None:
        return int(env) if int(env) > 0 else None
    name = getattr(backend, "name", "")
    name = name() if callable(name) else name
    return 1 if "ionq" in str(name).lower() else None


def make_sampler(args):
    mode = os.environ.get("QOLLAB_MODE", "local").lower()

    if mode == "local":
        from qiskit.primitives import StatevectorSampler
        return StatevectorSampler(default_shots=args.shots, seed=args.seed)

    if mode == "ionq-simulator":
        from qiskit.primitives import BackendSamplerV2
        from qiskit_ionq import IonQProvider
        token = os.environ.get("IONQ_API_KEY")
        if not token:
            raise SystemExit("Set IONQ_API_KEY to use QOLLAB_MODE=ionq-simulator")
        backend = IonQProvider(token).get_backend("simulator")
        backend.set_options(noise_model=IONQ_NOISE_MODEL)
        sampler = BackendSamplerV2(backend=backend, options={"default_shots": args.shots})
        return sampler, backend

    if mode == "qollab":
        from qiskit.primitives import BackendSamplerV2
        backend = _find_injected_backend()
        if backend is None:
            raise SystemExit("No global `backend` found. Run inside Qollab (choose the simulator in "
                             "'Run on Qollab / Select QPU') or use QOLLAB_MODE=local.")
        noise = os.environ.get("QOLLAB_NOISE_MODEL")
        if noise and hasattr(backend, "set_options"):
            backend.set_options(noise_model=noise)
        primitive = _JobGuard(BackendSamplerV2(backend=backend, options={"default_shots": args.shots}), _max_jobs(backend))
        return primitive, backend

    raise SystemExit(f"Unknown QOLLAB_MODE='{mode}'. Use local, ionq-simulator or qollab.")
