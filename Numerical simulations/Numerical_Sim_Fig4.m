% Numerical simulation,(controlled teleportation) average teleportation fidelity and teleportation
% success probability of: TP-EAM via W state, TP-EW via Bell state with q=0 and
% q=r, the MR framework and the original controlled teleportation protocols
% without any protection 

clear, clc
close all
syms alpha real
beta = sqrt(1-alpha^2);
jj = 0;
dd = 0.01;

for r = 0.01: dd: 1
    jj = jj + 1;
    
    %CTP-EW (Bell) q=r: Average teleportation fidelity
    q=r;
	fid_qr(jj)=1;
    %CTP-EW (Bell) q=r: Total teleportation success probability
    g_tot_qr(jj) =(q^2 - 2*q + r^2 - 2*r + 2)/(r^2 - 2*r + 2);
    
    
    %CTP-EW (Bell) q=0: Average teleporttaion fidelity
    q=0;
    fid1 =(alpha^2*q + beta^2*r - 1)^2/(alpha^2*(q^2 - 2*q)+ beta^2*(r^2 - 2*r) + 1);
	fid3 =(beta^2*q + alpha^2*r - 1)^2/(beta^2*(q^2 - 2*q)+ alpha^2*(r^2 - 2*r) + 1);
	P1 = (beta^2*r^2 - 2*beta^2*r + 1)/(2*(r^2 - 2*r + 2));
	P3 = (alpha^2*r^2 - 2*alpha^2*r + 1)/(2*(r^2 - 2*r + 2));
	ff = 2*(P1*fid1 + P3*fid3);
	FF = matlabFunction(ff);
	fid_r(jj)= integral(FF, 0, 1);
    %CTP-EW (Bell) q=0: Total teleportation success probability
    g_tot_r(jj) =(q^2 - 2*q + r^2 - 2*r + 2)/(r^2 - 2*r + 2);
    
    
    %CTP-EAM (W state) : Average teleporttaion fidelity
    fid_W(jj) = 1;
    %CTP-EAM (W state) : Total teleportation success probability
    g_tot_W(jj) = 1 ;
    
    %Original controlled teleportation without any protection (Bell state)
    fid_Bell_ori(jj) = (7*r^2)/15 - (11*r)/15 + 1;
    
    %Original controlled teleportation without any protection (W state)
    fid_W_ori(jj) = 1 - 11*r/15;
end

figure(1)
r = 0.01: dd: 1;
b1=plot(r, fid_qr, 'r--', 'LineWidth', 3); hold on
b2=plot(r, fid_r, 'k:', 'LineWidth', 2.5); hold on
b3=plot(r, fid_W, 'g', 'LineWidth', 1); hold on
b4=plot(r,  fid_Bell_ori, 'b-.', 'LineWidth', 2); hold on
b5=plot(r, fid_W_ori,'m--', 'LineWidth', 2); hold on

legend([b1 b2 b3 b4 b5],'Fid_{av}^{CTP-EW}(Bell)(q=r)', 'Fid_{av}^{CTP-EW}(Bell)(q=0)', 'Fid_{av}^{CTP-EAM}(W)(q=0)', 'Fid_{av}^{C-ori}(Bell)', 'Fid_{av}^{C-ori}(W)')
xlabel('\itr')
ylabel('Average teleportation fidelity')
set(gca,'fontsize',12,'fontname','Times');

figure(2)
r = 0.01: dd: 1;
b1=plot(r, g_tot_qr, 'r--', 'LineWidth', 2); hold on
b2=plot(r, g_tot_r, 'k:', 'LineWidth', 3.5); hold on
b3=plot(r, g_tot_W, 'g', 'LineWidth', 1.5); hold on


legend([b1 b2 b3],'g_{tot}^{CTP-EW}(Bell)(q=r)', 'g_{tot}^{CTP-EW}(Bell)(q=0)', 'g_{tot}^{CTP-EAM}(W)(q=0)')
xlabel('\itr')
ylabel('Total teleportation success probability')
set(gca,'fontsize',12,'fontname','Times');