OPENQASM 2.0;
include "qelib1.inc";
// cat4_prep: 4-qubit cat / GHZ state -- cat-state ladder (H then CNOT chain) -- used in Shor-style ancilla/flag prep
qreg q[4];
h q[0];
cx q[0],q[1];
cx q[1],q[2];
cx q[2],q[3];
