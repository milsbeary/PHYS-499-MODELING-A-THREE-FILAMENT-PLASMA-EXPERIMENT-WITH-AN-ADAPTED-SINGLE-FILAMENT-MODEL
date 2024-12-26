%Playing with force between magnetic dipoles 
%Last edited by M. H. Kent Oct 3, 2021

%We will assume that all dipoles will have the same moment in this case.
%We will investigate how the force will occilate as the moments rotate with
%a constant few radii. 

%Three dipoles with two thetas changeing
r=2;
r2=4;
p=10;
th1=0; %linspace(-pi,pi,100);
th2=linspace(0,4*pi,100);
th3=linspace(pi/4,6*pi/4,100);
th4=pi/4;
F1=((3./(r.^4)).*((p.^2)*cos(th1+th2)))+(((p).*cos(th1))+((p).*cos(th2))).*(15./(r.^6));
F2=((3./(r.^4)).*((p.^2)*cos(th3+th4)))+(((p).*cos(th3))+((p).*cos(th4))).*(15./(r.^6));
F=F1+F2;
plot(th2,F)
hold on 
plot(th3,F)
hold off
xlabel("theta")
ylabel("Magnitude")
title("Angle of dipole vs Force")
legend("P1 rotating [0,4pi] at 4pi(radians)/100(points)", "P2 rotating [pi/4,6pi/4], at 5pi/4(radians)/100(points)")