%2D Sydora Potential and E-Feild 
%Last edited by M. H. Kent Dec 26. 2021

%Set background of plots to white 
set (0, 'DefaultFigureColor', 'white' ); 

%Sets array for "r" value 
r = linspace(0,2,10000);

%Sydora Potential 
C0111=1; %Set to 1 to plot. Set to 0 to nullify curve 
C1=-3.309;
C2=-.690;
C3=5.712;
C4=0.397;
C5=-75.41;
PHI_BP = C1 + C2.*exp(-C3.*(r-C4).^2)+ C5.*(r + 3).^-4;
%Sydora E-Feild 
Er_Sy=C0111.*(-2.*C3.*C2.*(r-C4).*exp(-C3.*(r-C4).^2)-((4.*C5)./((r+3).^5)));

%Plot Sydora Potential 
figure(1)
plot(r, PHI_BP, 'p' ,col="Red")
%title("Sydora Potential")
xlabel('r (cm)')
ylabel('\delta (V)')

%Plots Sydora E-Feild 
figure(2)
plot(r, Er_Sy, 'p' ,col="Red")
%title("Sydora E-Feild ")%
xlabel('r (cm)')
ylabel('E (eV)')




