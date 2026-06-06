OPENQASM 2.0;
include "qelib1.inc";
// rep3_bitflip_encode: [[3,1,1]] bit-flip repetition -- encoder: CNOT fan-out; on |0> prepares |000>, on |1> prepares |111>
qreg q[3];
cx q[0],q[1];
cx q[0],q[2];
