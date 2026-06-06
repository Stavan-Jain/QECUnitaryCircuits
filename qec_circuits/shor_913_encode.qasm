OPENQASM 2.0;
include "qelib1.inc";
// shor_913_encode: [[9,1,3]] Shor -- full encoder (phase-flip outer then bit-flip inner); on |0> prepares |0_L>
qreg q[9];
cx q[0],q[3];
cx q[0],q[6];
h q[0];
h q[3];
h q[6];
cx q[0],q[1];
cx q[0],q[2];
cx q[3],q[4];
cx q[3],q[5];
cx q[6],q[7];
cx q[6],q[8];
