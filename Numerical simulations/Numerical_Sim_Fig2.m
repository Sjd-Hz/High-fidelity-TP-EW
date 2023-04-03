% Numerical simulation, average teleportation fidelity and teleportation
% success probability of TP-EW, original teleportation without protection
% and the MR framework

clear, clc
close all
syms alpha real
beta = sqrt(1-alpha^2);
jj = 0;
dd = 0.02;

for r = 0.001: dd: 0.999 % y axis
    jj = jj + 1;
    kk = 0;
    for q = 0.001: dd: 0.999 % x axis
        kk = kk + 1;
        
        %TP-EW: Average teleporttaion fidelity
        fid1 = -(alpha^4 - beta^4*r - alpha^4*q + beta^4 + 2*alpha^2*beta^2*(1 - q)^(1/2)*(1 - r)^(1/2))/((q - 1)*alpha^2 + (r - 1)*beta^2);
        fid3 = -(alpha^4 - alpha^4*r - beta^4*q + beta^4 + 2*alpha^2*beta^2*(1 - q)^(1/2)*(1 - r)^(1/2))/((r - 1)*alpha^2 + (q - 1)*beta^2);
        P1 = (1-beta^2*r )/(2*(2-r));
        P3 = (1-alpha^2*r )/(2*(2-r));
        ff = 2*(P1*fid1 + P3*fid3);
        FF = matlabFunction(ff);
        fid(jj, kk) = integral(FF, 0, 1);
        %TP-EW: Total teleportation success probability
        g_tot(jj, kk) = 1 - q/(2 - r);
        
        
        %MR framework: Average teleporttaion fidelity
        ff_MR = 1/2*((1 + r*alpha^2*beta^2)/(1 + r*alpha^2) + (1 + r*alpha^2*beta^2)/(1 + r*beta^2));
        FF_MR = matlabFunction(ff_MR);
        fid_MR(jj, kk) = integral(FF_MR, 0, 1);
        %MR framework: Total teleportation success probability
        g_tot_MR(jj, kk) = 1 - (r + r^2)/2;
        
        
        %Original teleportation without any protection
        fid_ori(jj, kk) = 1/15*(4*sqrt(1 - r) - 3.5*r + 11);
    end
end

q = 0.01: dd: 1;
r = 0.01: dd: 1;
%Plot average teleportation fidelity 
figure(1)
surf(q, r, fid); hold on 
surf(q, r, fid_MR); hold on
surf(q, r, fid_ori);

%The line corresponidng to maximum fildelity by setting q=r
z(1:100) = 1;
p = plot3(q, r, z(1: length(r)), 'w');
p.LineWidth = 3; 

xlabel('\itq')
ylabel('\itr')
zlabel('Average teleportation fidelity')   
colormap(jet)
colorbar
set(gca,'fontsize',12,'fontname','Times');
hold off


%Plot total teleporttaion success probability
figure(2)
surf(q, r, g_tot); hold on
surf(q, r, g_tot_MR);

xlabel('\itq')
ylabel('\itr')
zlabel('Total teleportation success probability')  
colormap(jet)
colorbar
set(gca,'fontsize',12,'fontname','Times');
hold off

