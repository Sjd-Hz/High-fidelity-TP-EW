% Numerical simulation, average teleportation fidelity and teleportation
% success probability of TP-EW for q=0 and q=r and the MR framework

clear, clc
close all
syms alpha real
beta = sqrt(1-alpha^2);
jj=0;
dd=0.02;
for r=0.01:dd:1
        jj=jj+1;
    
        %TP-EW (q=r): Average teleporttaion fidelity
        q=r;
        fid_qr(jj) =1;
        %TP-EW (q=r): Total teleportation success probability
        g_tot_qr(jj) = 1 - q/(2 - r);
  
        
        %TP-EW (q=0): Average teleporttaion fidelity
        q=0;
        if q==r 
            fid_r(jj)=1;
      % fidw(jj,kk)=1;
        else 
        fid1 = -(alpha^4 - beta^4*r - alpha^4*q + beta^4 + 2*alpha^2*beta^2*(1 - q)^(1/2)*(1 - r)^(1/2))/((q - 1)*alpha^2 + (r - 1)*beta^2);
        fid3 = -(alpha^4 - alpha^4*r - beta^4*q + beta^4 + 2*alpha^2*beta^2*(1 - q)^(1/2)*(1 - r)^(1/2))/((r - 1)*alpha^2 + (q - 1)*beta^2);
        P1 = (1-beta^2*r )/(2*(2-r));
        P3 = (1-alpha^2*r )/(2*(2-r));
        ff = 2*(P1*fid1 + P3*fid3);
        FF = matlabFunction(ff);
        fid_r(jj) = integral(FF, 0, 1);
        end
        %TP-EW (q=0): Total teleportation success probability
        g_tot_r(jj) = 1 - q/(2 - r);
        
        
        %MR framework: Average teleporttaion fidelity
        ff_MR = 1/2*((1 + r*alpha^2*beta^2)/(1 + r*alpha^2) + (1 + r*alpha^2*beta^2)/(1 + r*beta^2));
        FF_MR = matlabFunction(ff_MR);
        fid_MR(jj) = integral(FF_MR, 0, 1);
        %MR framework: Total teleportation success probability
        g_tot_MR(jj) = 1 - (r + r^2)/2;
end

r = 0.01: dd: 1;
figure(1)

b1=plot(r,fid_qr,'r','LineWidth',1.5','LineStyle','--');hold on
b2=plot(r,fid_r,'b','LineWidth',2.5,'LineStyle',':');hold on
b3=plot(r,fid_MR,'k','LineWidth',1.5);hold on

legend([b1 b2 b3],'Fid_{av}^{TP-EW}(q=r)','Fid_{av}^{TP-EW}(q=0)','Fid_{av}^{MR}')
axis tight
xlabel('\itr')
ylabel('Average teleportation fidelity')
set(gca,'fontsize',12,'fontname','Times');


figure(2)
b1=plot(r,g_tot_qr,'r','LineWidth',1.5,'LineStyle','--'); hold on
b2=plot(r,g_tot_r,'b','LineWidth',2.5,'LineStyle',':'); hold on
b3=plot(r,g_tot_MR,'k','LineWidth',1.5);hold on

legend([b1 b2 b3],'g_{tot}^{TP-EW}(q=r)','g_{tot}^{TP-EW}(q=0)','g_{tot}^{MR}')
axis tight
xlabel('\itr')
ylabel('Total teleportation success probability')
set(gca,'fontsize',12,'fontname','Times');


%}
