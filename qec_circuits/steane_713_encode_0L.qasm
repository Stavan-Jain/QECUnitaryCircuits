OPENQASM 2.0;
include "qelib1.inc";
// steane_713_encode_0L: [[7,1,3]] Steane -- logical |0_L> state-prep (H-seed + CNOT-spread over the 3 X-stabilizers)
qreg q[7];
h q[0];
h q[1];
h q[3];
cx q[0],q[2];
cx q[0],q[4];
cx q[0],q[6];
cx q[1],q[2];
cx q[1],q[5];
cx q[1],q[6];
cx q[3],q[4];
cx q[3],q[5];
cx q[3],q[6];
