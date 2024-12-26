%Mapping a magnetic dipole 

u = 4*pi*10.^-7; %Permativity of free space 
m = 10; %         %Magnetic dipole moment 

%Will be mapping in 3D plane but focousing on facing cordinates x and y

%Feild equations 
%By = (m./(((x.^2)+(y.^2))^(3./2))).*(3.*((y.^2)./((x.^2)+(y.^2)))-1);
%Bx = 3.*m.*((y.*x)./((x.^2)+(y.^2))).*(1./(((x.^2)+(y.^2))^(3./2)));

%Normilization vector/Feild density 
%B = sqrt((Bx.^2)+(By.^2));

%Plot feild density  
fig1 = figure(1)
[x,y] = meshgrid(-10:.01:10, -10:.01:10);
By = (1/4*pi).*abs(m./(((x.^2)+(y.^2))^(1.5))).*(3.*((y.^2)./((x.^2)+(y.^2)))-1);
Bx = (1/4*pi).*3.*m.*((y.*x)./((x.^2)+(y.^2))).*abs(1./(((x.^2)+(y.^2))^(1.5)));
B = sqrt((Bx.^2)+(By.^2));
feilddensity = pcolor(x,y,B);
feilddensity.FaceColor = 'interp';
set(feilddensity, 'EdgeColor', 'none');

hold on

%Plot vector feild 
%fig3 = figure(3)
[x,y] = meshgrid(-10:1:10, -10:1:10);
By = (1/4*pi).*abs(m./(((x.^2)+(y.^2))^(3./2))).*(3.*((y.^2)./((x.^2)+(y.^2)))-1);
Bx = (1/4*pi).*3.*m.*((y.*x)./((x.^2)+(y.^2))).*abs(1./(((x.^2)+(y.^2))^(3./2)));
feildvectors = quiver(x,y,Bx,By);
set(feildvectors,'LineWidth',.1,'Color','white');

hold off

%Plot streamline 
%fig2 = figure(2)
%Bottom of plot
%[x,y] = meshgrid(-10:.1:10, -10:.1:10);
%By = (1/4*pi).*abs(m./(((x.^2)+(y.^2))^(3./2))).*(3.*((y.^2)./((x.^2)+(y.^2)))-1);
%Bx = (1/4*pi).*3.*m.*((y.*x)./((x.^2)+(y.^2))).*abs(1./(((x.^2)+(y.^2))^(3./2)));
%startx = -10:.5:10;
%starty = zeros(size(startx));
%feildlines=streamline(x,y,Bx,By,x,y)
%feildlines=streamline(x,y,Bx,By,startx,starty);
%set(feildlines,'LineWidth',.2,'Color','red');
%hold on 
%Top of plot 
%[x,y] = meshgrid(-10:.1:10, -10:.1:10);
%By = -(1/4*pi).*abs(m./(((x.^2)+(y.^2))^(3./2))).*(3.*((y.^2)./((x.^2)+(y.^2)))-1);
%Bx = -(1/4*pi).*3.*m.*((y.*x)./((x.^2)+(y.^2))).*abs(1./(((x.^2)+(y.^2))^(3./2)));
%startx = -10:.5:10;
%starty = zeros(size(startx));
%feildlines=streamline(x,y,Bx,By,x,y)
%feildlines=streamline(x,y,Bx,By,startx,starty);
%set(feildlines,'LineWidth',.2,'Color','red');
%hold off 












