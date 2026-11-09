# W1: a closed parent for the calibrated record-action family

This is a new mathematical construction. The original Construction1 has an eight-level clock and six qubits, dimension512. W1 inserts an explicit environment register E and a parity coupling into that history. Its dimension is9 times2^7=1152. The original experiment and its verdict are unchanged.

Use work registers (S,M1,M2,A,E,D1,D2), P=Z_M1 Z_M2, hbar=1, and supplied angular frequency Omega=1. The initial state is

    |0>_clock |+>_S |Phi+>_M1,M2 |0>_A
      (sqrt((1+eta)/2)|0>_E + sqrt((1-eta)/2)|1>_E) |d1,d2>,

with0<=eta<=1 and selector bits d1,d2. The parameters theta1,theta2,beta,eta are physical model inputs. Eta is a coherence factor; it is distinct from the historical within-parity hold angle lambda=pi/2.

The eight unitary work links are:

1. CNOT(S,M2), the writer.
2. exp(-i*pi*Z_M1/4), the within-parity hold.
3. Rx_M2(theta1), controlled by D1=1.
4. U_EP=|0><0|_E tensor I + |1><1|_E tensor P.
5. exp(-i*beta*P/2), the relative parity phase.
6. Rx_M2(theta2), controlled by D2=1.
7. CNOT(M1,A) CNOT(M2,A), the parity-controlled actuator.
8. Identity, the readout link.

No later work link touches E. At each subsequent supported clock tick, tracing E implements

    Lambda_eta(rho)=(1+eta)rho/2+(1-eta)P rho P/2

at the inserted slot. This is an exact finite dilation, not an uncontrolled dissipative approximation. Eta1 leaves the work unchanged and eta0 implements nonselective parity dephasing. The finite parent can reverse this coupling on reverse clock motion, so no irreversible bath or universal Markov law is derived.

Define the time-independent joint Hamiltonian

    H=sum_{k=0}^7 sqrt((k+1)(8-k))/2
      [|k+1><k| tensor U_k + |k><k+1| tensor U_k^dagger].

Every link remains present throughout the unmeasured evolution. Put V0=I, V_{k+1}=U_k V_k and T=sum_k |k><k| tensor V_k. Direct multiplication gives H=T(H_clock tensor I)T^dagger, where H_clock is J_x for spin4. Thus at every nonzero-weight tick k the work state is exactly V_k rho_work V_k^dagger. The clock probability is C(8,k)p^k(1-p)^(8-k), p=sin^2(t/2); at t=pi the last tick has unit probability. This is the standard engineered spin-chain clock construction, here verified directly; see [Christandl et al.](https://arxiv.org/abs/quant-ph/0309131).

Change the memory basis from (m1,m2) to (g,p)=(m1,m1 xor m2). Then the parity operators tau_x=X_M2, tau_y=Z_M1 Y_M2, tau_z=P act on p, while the hold acts on g. The Bell blank becomes |+>_g|0>_p. On the declared clock-conditioned Z_S,P,Z_A readouts, the effective work model is therefore three logical qubits (S,p,A) plus the explicitly traced dephasing ancilla. This reduction does not preserve arbitrary complementary-memory observables or the parent's full clock coherence.

For the primary preparation, writer, two enabled pulses and actuator, the final prediction is

    Pr(Z_A=Z_S) = [1+cos(theta1)cos(theta2)
                     -eta sin(theta1)sin(theta2)cos(beta)]/2.

A disabled pulse sets its theta to zero. At theta1=theta2=pi/3 and beta=2pi/3 this is5/8+3eta/16. The memoryless population comparator is the eta0 value. The ability to distinguish that comparator comes from a fixed phase-sensitive prediction, not from post-hoc parameter fitting.

The production consistency test independently assembles H as a sparse1152-dimensional matrix, propagates it at the registered times, and compares its clock and conditional readout distributions to the ordered work circuit. The ordered final probability is also compared to the analytic formula. The numerical energy check covers mean energy; conservation of the full energy distribution follows exactly from time-independent unitary dynamics.

The physical access table of this parent must admit the stated clock-work links and readouts. Their engineering, locality and calibration are not supplied by the unitary equivalence. The proposed four-qubit hardware experiment tests the externally scheduled work protocol; it does not by itself validate this larger autonomous Hamiltonian. Its pulse, preparation and readout errors remain subject to the M1 calibration and held-out acceptance rules.
