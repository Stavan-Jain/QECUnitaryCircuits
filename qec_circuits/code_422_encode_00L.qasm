OPENQASM 2.0;
include "qelib1.inc";
// code_422_encode_00L: [[4,2,2]] error-detecting -- logical |00_L> state-prep (GHZ-like): H on q0 then CNOT chain
qreg q[4];
h q[0];
cx q[0],q[1];
cx q[0],q[2];
cx q[0],q[3];
