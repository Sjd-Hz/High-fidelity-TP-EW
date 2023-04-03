% Original Controlled teleportation with W state through ADC
% Derive the final expressions of total teleporttaion success probability and teleporttaion fidelity 

clear, clc
close all

% input state parameters
syms alpha beta real 
assume((0<=alpha)&(alpha<=1))
assume((0<=beta)&(beta<=1))
assume(alpha^2+beta^2==1)
rho_in = [alpha^2, alpha*beta; alpha*beta, beta^2];

% The decaying rate of ADC (r) and the weak measurement strength (q)
syms r q  real 
assume(r,'real')
assume(q,'real')
assume((0<=r)&(r<=1))
assume((0<=q)&(q<=1))

%Pauli operators
sigmax = [0 1; 1 0];
sigmaz = [1 0; 0 -1];

Id2 = eye(2,2); %Identity operator 
H = [1;0]; %|0>
V = [0;1]; %|1>

%W entangled state
psi_W = 1/2*(kron(V, kron(H, H)) + kron(H, kron(V, H)) + sqrt(2)*kron(H, kron(H, V))); 
rho_W = psi_W*psi_W';

%ADC Kraus operators
e0 = [1 0; 0 sqrt(1 - r)];
e1 = [0 sqrt(r); 0 0];

%ADC Kraus operators in case of teleportation with W state
E0 = kron(e0, kron(e0, e0));
E1 = kron(e0, kron(e0, e1));
E2 = kron(e0, kron(e1, e0));
E3 = kron(e0, kron(e1, e1));
E4 = kron(e1, kron(e0, e0));
E5 = kron(e1, kron(e0, e1));
E6 = kron(e1, kron(e1, e0));
E7 = kron(e1, kron(e1, e1));

% Normalized W entangled state after passing through ADC
rho_W_damp = (E0*rho_W*E0'+E1*rho_W*E1'+E2*rho_W*E2'+E3*rho_W*E3'+E4*rho_W*E4'+E5*rho_W*E5'+E6*rho_W*E6'+E7*rho_W*E7')/trace(E0*rho_W*E0'+E1*rho_W*E1'+E2*rho_W*E2'+E3*rho_W*E3'++E4*rho_W*E4'+E5*rho_W*E5'+E6*rho_W*E6'+E7*rho_W*E7');

%interacting the input state with W shared entangled state at Alice's end
rho_com = kron(rho_in, rho_W_damp); 

%joint measurement operators at Alice's in case of W state
psi_P1 = 1/2*(kron(H, kron(V, H)) + kron(H, kron(H, V)) + sqrt(2)*kron(V, kron(H, H)));
psi_P2 = 1/2*(kron(H, kron(V, H)) + kron(H, kron(H, V)) - sqrt(2)*kron(V, kron(H, H)));
psi_P3 = 1/2*(kron(V, kron(V, H)) + kron(V, kron(H, V)) + sqrt(2)*kron(H, kron(H, H)));
psi_P4 = 1/2*(kron(V, kron(V, H)) + kron(V, kron(H, V)) - sqrt(2)*kron(H, kron(H, H)));
rho_P1 = psi_P1*psi_P1';
rho_P2 = psi_P2*psi_P2';
rho_P3 = psi_P3*psi_P3';
rho_P4 = psi_P4*psi_P4';

% Alice applies measurements on the input state and her qubit of the shared W entangled state
Proj1 = kron(rho_P1, Id2);
Proj2 = kron(rho_P2, Id2);
Proj3 = kron(rho_P3, Id2);
Proj4 = kron(rho_P4, Id2);

% Probability of occurance of each measurement operator
P1 = simplify(trace(Proj1*rho_com*Proj1'/trace(rho_com)));
P2 = simplify(trace(Proj2*rho_com*Proj2'/trace(rho_com)));
P3 = simplify(trace(Proj3*rho_com*Proj3'/trace(rho_com)));
P4 = simplify(trace(Proj4*rho_com*Proj4'/trace(rho_com)));

% Bob's state corresponidng to each measurement outcome of Alice
rho_Proj1 = PartialTrace(Proj1*rho_com*Proj1', [1,2,3]);
rho_Proj2 = PartialTrace(Proj2*rho_com*Proj2', [1,2,3]);
rho_Proj3 = PartialTrace(Proj3*rho_com*Proj3', [1,2,3]);
rho_Proj4 = PartialTrace(Proj4*rho_com*Proj4', [1,2,3]);

% Unitary operations corresponding to each measurement outcome of Alice
U1 = Id2;
U2 = sigmaz;
U3 = sigmax;
U4 = sigmax*sigmaz;

% Non normalized output states at Bob's end corresponding to each measurement outcome of Alice
rho_U1 = U1*rho_Proj1*U1';
rho_U2 = U2*rho_Proj2*U2';
rho_U3 = U3*rho_Proj3*U3';
rho_U4 = U4*rho_Proj4*U4';

% Success probability of obtaining each final state 
g1 = trace(rho_U1);
g2 = trace(rho_U2);
g3 = trace(rho_U3);
g4 = trace(rho_U4);

%%%%%%%% Total teleportation success probability  %%%%%%%
g_tot = simplify(2*(g1 + g3))

% Normalized output states at Bob's end corresponding to each measurement outcome of Alice
rho_U1_nor = rho_U1/g1;
rho_U2_nor = rho_U2/g2;
rho_U3_nor = rho_U3/g3;
rho_U4_nor = rho_U4/g4;

% Fidelity between input state and each output state 
fid1 = simplify(trace(rho_in*rho_U1_nor));
fid2 = simplify(trace(rho_in*rho_U2_nor));
fid3 = simplify(trace(rho_in*rho_U3_nor));
fid4 = simplify(trace(rho_in*rho_U4_nor));

%%%%%%%%%  average teleportation fidlelity %%%%%%%%%
fid_tot=simplify(2*(P1*fid1+P3*fid3));
Fid_av = int(fid_tot,beta, 0, 1)

