%Test of potential curves 

r = linspace(0,2,100);

%Sydora potential 
C0111=1
C1=-3.309;
C2=-.690;
C3=5.712;
C4=0.397;
C5=-75.41;
PHI_BP = C1 + C2.*exp(-C3.*(r-C4).^2)+ C5.*(r + 3).^-4;
%Sydora feild 
Er_Sy=C0111.*(-2.*C3.*C2.*(r-C4).*exp(-C3.*(r-C4).^2)-((4.*C5)./((r+3).^5)));

%Test potential 
C011=1
Acon=.25;
ConPHI=2;
C=1;
PHI_BP1=ConPHI*(C-Acon.*(r.^2));
%Test feild
Er_tp=C011.*(2.*ConPHI.*Acon.*r);

%Plots of potentials 
figure(1)
plot(r, PHI_BP)
title("Potentials")
hold on 
plot(r, PHI_BP1)
hold off 
legend("Sydora potential", "Test potential")

%Plots of feilds 
figure(2)
plot(r,Er_Sy)
title("Electric Feilds")
hold on 
plot(r,Er_tp)
hold off
legend("Sydora feild", "Test feild")



