# QEC purely-unitary origin circuits -- benchmark set for Manjushri

OpenQASM 2.0, Clifford+T, no measurement/reset/classical control. Each circuit is a single
"origin" input: the Manjushri harness applies PyZX to produce the equivalent "optimized" twin
(and a gate-deleted mutant) for equivalence/inequivalence checking. Every circuit below was
verified by direct statevector simulation (stabilizer eigenstate / logical-action checks).

| file | code | qubits | gates | verify | reference |
|------|------|:-----:|:----:|:----:|-----------|
| steane_713_H_nf_fig4.qasm | [[7,1,3]] Steane | 7 | 18 | PASS | Chamberland & Cross, arXiv:1811.00566 Fig. 4 (logical |H>) |
| steane_713_encode_0L.qasm | [[7,1,3]] Steane | 7 | 12 | PASS | Steane, PRL 77, 793 (1996); standard CSS encoder (ECZoo c/steane) |
| steane_713_encode_plusL.qasm | [[7,1,3]] Steane | 7 | 19 | PASS | Steane code; |+_L> = transversal H on |0_L> |
| shor_913_encode.qasm | [[9,1,3]] Shor | 9 | 11 | PASS | Shor, PRA 52, R2493 (1995); concatenated phase/bit-flip encoder |
| rep3_bitflip_encode.qasm | [[3,1,1]] bit-flip repetition | 3 | 2 | PASS | 3-qubit bit-flip repetition code (textbook) |
| rep3_phaseflip_encode.qasm | [[3,1,1]] phase-flip repetition | 3 | 5 | PASS | 3-qubit phase-flip repetition code (textbook) |
| code_422_encode_00L.qasm | [[4,2,2]] error-detecting | 4 | 4 | PASS | [[4,2,2]] code (Vaidman-Goldenberg-Wiesner 1996; Grassl 1997) |
| cat3_prep.qasm | 3-qubit cat / GHZ state | 3 | 3 | PASS | GHZ/cat state; Shor-style ancilla & flag verification prep |
| cat4_prep.qasm | 4-qubit cat / GHZ state | 4 | 4 | PASS | GHZ/cat state (ladder construction) |
| cat5_prep.qasm | 5-qubit cat / GHZ state | 5 | 5 | PASS | GHZ/cat state (ladder construction) |
| cat7_prep.qasm | 7-qubit cat / GHZ state | 7 | 7 | PASS | GHZ/cat state (ladder construction) |
| ghz5_fanout.qasm | 5-qubit GHZ state | 5 | 5 | PASS | GHZ state via H + CNOT fan-out (alternate construction to cat ladder) |
| steane_713_transversal_H.qasm | [[7,1,3]] Steane | 7 | 7 | PASS | Steane self-dual CSS: transversal H = logical H |
| steane_713_transversal_S.qasm | [[7,1,3]] Steane | 7 | 7 | PASS | Steane: transversal S^dagger = logical S (verified codespace-preserving) |
| steane_713_transversal_CNOT.qasm | [[7,1,3]] Steane x2 | 14 | 7 | PASS | Transversal logical CNOT between two Steane blocks (CSS) |
| rm15_1531_encode_0L.qasm | [[15,1,3]] Reed-Muller (RM(1,4)) | 15 | 32 | PASS | [[15,1,3]] punctured RM(1,4); Knill-Laflamme-Zurek; CSS encoder |
| rm15_1531_transversal_Tdg.qasm | [[15,1,3]] Reed-Muller | 15 | 15 | PASS | [[15,1,3]] RM: transversal T^dagger = logical T (triorthogonal; weights in {0,8}) |
| code_513_encode_0L.qasm | [[5,1,3]] perfect (non-CSS) | 5 | 31 | PASS | Laflamme-Miquel-Paz-Zurek, quant-ph/9602019; [[5,1,3]] perfect code (non-CSS) |
| code_833_encode_basis.qasm | [[8,3,3]] Gottesman (non-CSS, k=3) | 8 | 52 | PASS | Gottesman thesis quant-ph/9705052, Table 3.3; [[8,3,3]] (non-CSS, k=3) |

Total: 19 circuits, all verified. Regenerate with `python3 generate_qec_circuits.py`.
