OPENQASM 2.0;
include "qelib1.inc";
// cat5_prep: 5-qubit cat / GHZ state -- cat-state ladder (H then CNOT chain) -- used in Shor-style ancilla/flag prep
qreg q[5];
h q[0];
cx q[0],q[1];
cx q[1],q[2];
cx q[2],q[3];
cx q[3],q[4];
