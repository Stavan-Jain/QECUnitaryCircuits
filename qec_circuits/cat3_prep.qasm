OPENQASM 2.0;
include "qelib1.inc";
// cat3_prep: 3-qubit cat / GHZ state -- cat-state ladder (H then CNOT chain) -- used in Shor-style ancilla/flag prep
qreg q[3];
h q[0];
cx q[0],q[1];
cx q[1],q[2];
