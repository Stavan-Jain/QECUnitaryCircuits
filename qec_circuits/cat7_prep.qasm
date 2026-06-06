OPENQASM 2.0;
include "qelib1.inc";
// cat7_prep: 7-qubit cat / GHZ state -- cat-state ladder (H then CNOT chain) -- used in Shor-style ancilla/flag prep
qreg q[7];
h q[0];
cx q[0],q[1];
cx q[1],q[2];
cx q[2],q[3];
cx q[3],q[4];
cx q[4],q[5];
cx q[5],q[6];
