# QEC Unitary Circuits

A benchmark set of **purely unitary Quantum Error Correction circuits** in OpenQASM 2.0, for testing
quantum-circuit **equivalence checkers** such as
[Manjushri](https://arxiv.org/html/2601.22372v1).

All circuits are **Clifford+T** with **no measurement, reset, or classical control**, so they are
valid inputs to equivalence-checking tools that compare two unitaries up to global phase. Each
circuit is a single *origin* circuit: a checker's harness can run an optimizer (e.g. PyZX) to
produce an equivalent "optimized" twin and verify the two match (and a gate-deleted mutant to
confirm inequivalence is caught).

The seed example is Chamberland & Cross, *Fault-tolerant magic state preparation with flag qubits*
([arXiv:1811.00566](https://arxiv.org/abs/1811.00566)), **Fig. 4** — a non-fault-tolerant logical
|H⟩ preparation for the [[7,1,3]] Steane code.

## Contents

| Path | Description |
|------|-------------|
| [`qec_circuits/`](qec_circuits/) | 19 verified `.qasm` circuits |
| [`qec_circuits/MANIFEST.md`](qec_circuits/MANIFEST.md) | per-circuit table: code, qubits, gate count, verification, reference |
| [`generate_qec_circuits.py`](generate_qec_circuits.py) | pure-numpy generator + verifier + stabilizer synthesizer |
| [`QEC_unitary_circuits_catalog.md`](QEC_unitary_circuits_catalog.md) | literature catalog of 45 cited circuits (shortlist, equivalence pairs, gaps) |
| `_verified_circuits.json` | machine-readable backing data for the catalog |

## The circuits (19)

**Encoders / logical state preparation** — Steane [[7,1,3]] (|0_L⟩, |+_L⟩, |H⟩), Shor [[9,1,3]],
[[4,2,2]], [[5,1,3]] perfect (non-CSS), [[8,3,3]] (non-CSS, k=3), [[15,1,3]] Reed–Muller, 3-qubit
bit-/phase-flip repetition.

**Stabilizer / cat / GHZ states** — cat/GHZ ladders (n = 3,4,5,7) and a GHZ fan-out construction.

**Transversal / logical gates** — Steane transversal H, S, and CNOT (between two blocks); [[15,1,3]]
transversal T† (= logical T).

Coverage: CSS and non-CSS codes; k = 1, 2, 3; n up to 15.

## Verification

Every circuit is checked **by direct statevector simulation before it is written** — nothing is
taken on trust:

- **Encoders / state prep:** the prepared state is confirmed to be a +1 eigenstate of every
  stabilizer generator, with the intended logical-operator value (e.g. ⟨Z_L⟩ = +1 for |0_L⟩).
- **Transversal / logical gates:** verified on *encoded* states (the transversal layer applied to
  |0…0⟩ is not a codeword). E.g. transversal H maps |0_L⟩ → |+_L⟩; transversal T† on [[15,1,3]]
  fixes |0_L⟩, consistent with all |0_L⟩ basis-state weights being ≡ 0 mod 8 (triorthogonality).

Non-CSS encoders ([[5,1,3]], [[8,3,3]]) are produced by a general **graph-state stabilizer
synthesizer** (Hadamard subset → X-block = I → CZ/S for the symmetric Z-block → Pauli sign-fix),
cross-checked against known CSS codes first.

## Reproduce

```bash
python3 generate_qec_circuits.py    # requires only numpy; regenerates qec_circuits/ + MANIFEST.md
```

## References & license

The circuits are standard constructions from the quantum error correction literature; see
[`qec_circuits/MANIFEST.md`](qec_circuits/MANIFEST.md) and
[`QEC_unitary_circuits_catalog.md`](QEC_unitary_circuits_catalog.md) for per-circuit references.

No license file is included yet — add one (e.g. MIT, Apache-2.0, or a public-domain dedication) to
set reuse terms.
