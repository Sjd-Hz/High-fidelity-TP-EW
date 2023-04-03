% Teleportation with Bell state via weak measurement and environment-assisted measurement
% Derive the final expressions of total teleporttaion success probability and teleporttaion fidelity 

clear, clc
close all

% input state parameters
syms alpha beta real 
assume((0<=alpha)&(alpha<=1))
assume((0<=beta)&(beta<=1))
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

%Bell entangled state
psi_Bell = 1/sqrt(2)*(kron(H, H) + kron(V, V));
rho_Bell = psi_Bell*psi_Bell';

%ADC Kraus operators
e0 = [1 0; 0 sqrt(1 - r)];
e1 = [0 sqrt(r); 0 0];

%ADC Kraus operators in case of teleportation with Bell state
E0 = kron(Id2, e0);
E1 = kron(Id2, e1);

% Bell entangled state after passing through ADC by employing EAM
rho_Bell_damp = (E0*rho_Bell*E0')/trace(E0*rho_Bell*E0');

%Success probability of EAM in entanglement distribution process
g_EAM = trace(rho_Bell_damp);

%interacting the input state with Alice's qubit of shared Bell entangled state
rho_com = kron(rho_in, rho_Bell_damp);

%joint measurement operators at Alice's in case of Bell state
b1 = 1/sqrt(2)*(kron(H, H) + kron(V, V));
b2 = 1/sqrt(2)*(kron(H, H) - kron(V, V));
b3 = 1/sqrt(2)*(kron(H, V) + kron(V, H));
b4 = 1/sqrt(2)*(kron(H, V) - kron(V, H));


% Alice applies measurements on the input state and her qubit of the shared Bell entangled state
B1 = kron(b1*b1', Id2);
B2 = kron(b2*b2', Id2);
B3 = kron(b3*b3', Id2);
B4 = kron(b4*b4', Id2);

% Probability of occurance of each measurement operator
P1 = simplify(trace(B1*rho_com/trace(rho_com)*B1'))
P2 = simplify(trace(B2*rho_com/trace(rho_com)*B2'));
P3 = simplify(trace(B3*rho_com/trace(rho_com)*B3'))
P4 = simplify(trace(B2*rho_com/trace(rho_com)*B4'));

% Bob's state corresponidng to each measurement outcome of Alice
rho_B1 = simplify(PartialTrace(B1*rho_com*B1', [1,2]));
rho_B2 = simplify(PartialTrace(B2*rho_com*B2', [1,2]));
rho_B3 = simplify(PartialTrace(B3*rho_com*B3', [1,2]));
rho_B4 = simplify(PartialTrace(B4*rho_com*B4', [1,2]));

% Designed weak measurement operator
M1 = [sqrt(1 - q), 0; 0, 1];
M2 = [sqrt(q), 0; 0, 0];

% Weak measurement operator corresponding to each measurement outcome of Alice
M1_1 = M1;
M1_2 = sigmaz*M1;
M1_3 = sigmax*M1;
M1_4 = sigmax*sigmaz*M1;

% Non normalized output states at Bob's end corresponding to each measurement outcome of Alice
rho_WM1 = M1_1*rho_B1*M1_1';
rho_WM2 = M1_2*rho_B2*M1_2';
rho_WM3 = M1_3*rho_B3*M1_3';
rho_WM4 = M1_4*rho_B4*M1_4';

% Success probability of obtaining each final state 
g1 = trace(rho_WM1);
g2 = trace(rho_WM2);
g3 = trace(rho_WM3);
g4 = trace(rho_WM4);

%%%%%%%% Total teleportation success probability  %%%%%%%
g_tot = simplify(2*(g1 + g3))

% Normalized output states at Bob's end corresponding to each measurement outcome of Alice
rho_WM1_nor = rho_WM1/g1;
rho_WM2_nor = rho_WM2/g2;
rho_WM3_nor = rho_WM3/g3;
rho_WM4_nor = rho_WM4/g4;

%%%%%%%%% Fidelity between input state and each output state %%%%%%%%
fid1 = simplify(trace(rho_in*rho_WM1_nor))
fid2 = simplify(trace(rho_in*rho_WM2_nor))
fid3 = simplify(trace(rho_in*rho_WM3_nor))
fid4 = simplify(trace(rho_in*rho_WM4_nor))
