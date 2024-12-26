%% Single Mode Model 15 mm 
%Window set V=[0,1]
%Flow feild not valid 
%Last edited by M. H. Kent Dec 26th, 2021

clc
clear

load('modeturbo.mat');

make_fig2 = 0;
make_movie = 1;
speed = 3;
%% 
% Grid
x = linspace(-2.25,2.25,256);
y = linspace(-2.25,2.25,256);
r = linspace(0,2,256);
[X,Y] = meshgrid(x,y);

%3 Filament Rotation and Distance   
th=pi/3.5; %Rotation from x=0
d = 1.5; %Distance of filaments from each other (in cm)  
rr= d./sqrt(3);
%rr = .5;   %Distance of filaments from center of grid

% Time and Freqency of exparament
f1 = 35e3; %Will calculate base off of SHI filament 
dt = .5*3.2e-7; % Common time step in LAPD experiments
N = 5;
t = 0:dt:(N/f1); % Time array of N m1 periods
freq = (0:length(t)-1)/length(t)/dt/1e3;

%Magnetic Field: used to calculate force on particals  
B = .1; % Tesla (1000 G)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SHI FILAMENT (Not being used but can be added)

% Polar Coordinates
R = (X.^2 + Y.^2).^.5;
TH = atan2(Y,X);

% Azimuthal Mode Numbers
m1 = 0;
m2 = 6;
M = [m1, m2];

% Mode Amplitudes
A1 = 0;
A2 = 0;
A = [A1, A2];
% Radial Wavenumbers
k1 = 10; 
k2 = 15;
K = [k1, k2];

% Decay Parameter
alpha1 = 1;
alpha2 = 1;
Alpha = [alpha1, alpha2];

% Mode Frequencies
f1 = 35e3;
f2 = 35e3;
w1 = 2*pi*f1;
w2 = 2*pi*f2;
W = [w1, w2];

% Magnetic Field: used to calculate force on particals 
%B = .1; % Tesla (1000 G)

% Time and Freqency
%dt = .5*3.2e-7; % Common time step in LAPD experiments
%N = 5;
%t = 0:dt:(N/f1); % Time array of N m1 periods
%freq = (0:length(t)-1)/length(t)/dt/1e3;
    
% Potentials (still complex valued)
PHI_1t = A1.*besselj(m1,k1.*R).*exp(m1*1i.*TH).*exp(-alpha1.*R);
PHI_2t = A2.*besselj(m2,k2.*R).*exp(m2*1i.*TH).*exp(-alpha2.*R);

% E-Fields
Er_1  = -A1.*exp(m1*1i.*TH).*exp(-alpha1.*R).*(...
            besselj_p(m1,R,k1) - alpha1.*besselj(m1,k1.*R));
Er_2  = -A2.*exp(m2*1i.*TH).*exp(-alpha2.*R).*(...
            besselj_p(m2,R,k2) - alpha2.*besselj(m2,k2.*R));
        
Eth_1 = -1i.*m1.*A1.*exp(m1.*1i.*TH).*exp(-alpha1.*R).*besselj(m1,k1.*R)./R;
Eth_2 = -1i.*m2.*A2.*exp(m2.*1i.*TH).*exp(-alpha2.*R).*besselj(m2,k2.*R)./R;

Ex_1 = cos(TH).*Er_1-sin(TH).*Eth_1;
Ex_2 = cos(TH).*Er_2-sin(TH).*Eth_2;

Ey_1 = sin(TH).*Er_1+cos(TH).*Eth_1;
Ey_2 = sin(TH).*Er_2+cos(TH).*Eth_2;

% Flow Field
%U =  (real(Ey_1) + real(Ey_2))/B*1e3; %cm/s
%V = -(real(Ex_1) + real(Ex_2))/B*1e3; %cm/s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FILAMENT 1

% Polar Coordinates
X1=rr.*cos(th); 
Y1=rr.*sin(th);
R1 = ((X-X1).^2 + (Y-Y1).^2).^.5;
TH1 = atan2((Y-Y1),(X-X1));

% Azimuthal Mode Numbers
m3 = 0;
m4 = 6;
M1 = [m3, m4];

% Mode Amplitudes
A3 = 1;
A4 = 0;
A1 = [A3, A4];
% Radial Wavenumbers
k3 = 8;
k4 = 15;
K1 = [k3, k4];

% Decay Parameter
alpha3 = 3;
alpha4 = 1;
Alpha1 = [alpha3, alpha4];

% Mode Frequencies
f3 = 35e3;
f4 = 35e3;
w3 = 2*pi*f3;
w4 = 2*pi*f4;
W1 = [w3, w4];

% Magnetic Field: used to calculate force on particals 
%B = .1; % Tesla (1000 G)

% Time and Freqency
%dt = .5*3.2e-7; % Common time step in LAPD experiments
%N = 5;
%t = 0:dt:(N/f1); % Time array of N m1 periods
%freq = (0:length(t)-1)/length(t)/dt/1e3;
    
% Potentials (still complex valued)
PHI_3 = A3.*besselj(m3,k3.*R1).*exp(m3*1i.*TH1).*exp(-alpha3.*R1);
PHI_4 = A4.*besselj(m4,k4.*R1).*exp(m4*1i.*TH1).*exp(-alpha4.*R1);

% E-Fields
Er_3  = -A3.*exp(m3*1i.*TH1).*exp(-alpha3.*R1).*(...
            besselj_p(m3,R1,k3) - alpha3.*besselj(m3,k3.*R1));
Er_4  = -A4.*exp(m4*1i.*TH1).*exp(-alpha4.*R1).*(...
            besselj_p(m4,R1,k4) - alpha4.*besselj(m4,k4.*R1));
        
Eth_3 = -1i.*m3.*A3.*exp(m3.*1i.*TH1).*exp(-alpha3.*R1).*besselj(m3,k3.*R1)./R1;
Eth_4 = -1i.*m4.*A4.*exp(m4.*1i.*TH1).*exp(-alpha4.*R1).*besselj(m4,k4.*R1)./R1;

Ex_3 = cos(TH1).*Er_3-sin(TH1).*Eth_3;
Ex_4 = cos(TH1).*Er_4-sin(TH1).*Eth_4;

Ey_3 = sin(TH1).*Er_3+cos(TH1).*Eth_3;
Ey_4 = sin(TH1).*Er_4+cos(TH1).*Eth_4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FILAMENT 2

% Polar Coordinates
X2=rr.*cos(th + ((2*pi)/3)); 
Y2=rr.*sin(th + ((2*pi)/3));
R2 = ((X-X2).^2 + (Y-Y2).^2).^.5;
TH2 = atan2((Y-Y2),(X-X2));

% Azimuthal Mode Numbers
m5 = 0;
m6 = 6;
M2 = [m5, m6];

% Mode Amplitudes
A5 = 1;
A6 = 0;
A2 = [A5, A6];
% Radial Wavenumbers
k5 = 8; 
k6 = 15;
K2 = [k5, k6];

% Decay Parameter
alpha5 = 3;
alpha6 = 1;
Alpha2 = [alpha5, alpha6];

% Mode Frequencies
f5 = 35e3;
f6 = 35e3;
w5 = 2*pi*f5;
w6 = 2*pi*f6;
W2 = [w5, w6];

% Magnetic Field: used to calculate force on particals 
%B = .1; % Tesla (1000 G)

% Time and Freqency
%dt = .5*3.2e-7; % Common time step in LAPD experiments
%N = 5;
%t = 0:dt:(N/f1); % Time array of N m1 periods
%freq = (0:length(t)-1)/length(t)/dt/1e3;
    
% Potentials (still complex valued)
PHI_5 = A5.*besselj(m5,k5.*R2).*exp(m5*1i.*TH2).*exp(-alpha5.*R2);
PHI_6 = A6.*besselj(m6,k6.*R2).*exp(m6*1i.*TH1).*exp(-alpha6.*R2);

% E-Fields
Er_5  = -A5.*exp(m5*1i.*TH2).*exp(-alpha5.*R2).*(...
            besselj_p(m5,R2,k5) - alpha5.*besselj(m5,k5.*R2));
Er_6  = -A6.*exp(m6*1i.*TH2).*exp(-alpha6.*R2).*(...
            besselj_p(m6,R2,k6) - alpha6.*besselj(m6,k6.*R2));
        
Eth_5 = -1i.*m5.*A5.*exp(m5.*1i.*TH2).*exp(-alpha5.*R2).*besselj(m5,k5.*R2)./R2;
Eth_6 = -1i.*m6.*A6.*exp(m6.*1i.*TH2).*exp(-alpha6.*R2).*besselj(m6,k6.*R2)./R2;

Ex_5 = cos(TH2).*Er_5-sin(TH2).*Eth_5;
Ex_6 = cos(TH2).*Er_6-sin(TH2).*Eth_6;

Ey_5 = sin(TH2).*Er_5+cos(TH2).*Eth_5;
Ey_6 = sin(TH2).*Er_6+cos(TH2).*Eth_6;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILAMENT 3

% Polar Coordinates
X3=rr.*cos(th + ((4*pi)/3)); 
Y3=rr.*sin(th + ((4*pi)/3));
R3 = ((X-X3).^2 + (Y-Y3).^2).^.5;
TH3 = atan2((Y-Y3),(X-X3));

% Azimuthal Mode Numbers
m7 = 0;
m8 = 6;
M3 = [m7, m8];

% Mode Amplitudes
A7 = 1;
A8 = 0;
A3 = [A7, A8];
% Radial Wavenumbers
k7 = 8; 
k8 = 15;
K3 = [k7, k8];

% Decay Parameter
alpha7 = 3;
alpha8 = 1;
Alpha3 = [alpha7, alpha8];

% Mode Frequencies
f7 = 35e3;
f8 = 35e3;
w7 = 2*pi*f7;
w8 = 2*pi*f8;
W3 = [w7, w8];

% Magnetic Field: used to calculate force on particals 
%B = .1; % Tesla (1000 G)

% Time and Freqency
%dt = .5*3.2e-7; % Common time step in LAPD experiments
%N = 5;
%t = 0:dt:(N/f1); % Time array of N m1 periods
%freq = (0:length(t)-1)/length(t)/dt/1e3;
    
% Potentials (still complex valued)
PHI_7 = A7.*besselj(m7,k7.*R3).*exp(m7*1i.*TH3).*exp(-alpha7.*R3);
PHI_8 = A8.*besselj(m8,k8.*R3).*exp(m8*1i.*TH3).*exp(-alpha8.*R3);

% E-Fields
Er_7  = -A7.*exp(m7*1i.*TH3).*exp(-alpha7.*R3).*(...
            besselj_p(m7,R3,k7) - alpha7.*besselj(m7,k7.*R3));
Er_8  = -A8.*exp(m8*1i.*TH3).*exp(-alpha8.*R3).*(...
            besselj_p(m8,R3,k8) - alpha8.*besselj(m8,k8.*R3));
        
Eth_7 = -1i.*m7.*A7.*exp(m7.*1i.*TH3).*exp(-alpha7.*R3).*besselj(m7,k7.*R3)./R3;
Eth_8 = -1i.*m8.*A8.*exp(m8.*1i.*TH3).*exp(-alpha8.*R3).*besselj(m8,k8.*R3)./R3;

Ex_7 = cos(TH3).*Er_7-sin(TH3).*Eth_7;
Ex_8 = cos(TH3).*Er_8-sin(TH3).*Eth_8;

Ey_7 = sin(TH3).*Er_7+cos(TH3).*Eth_7;
Ey_8 = sin(TH3).*Er_8+cos(TH3).*Eth_8;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%POTENTIAL WELL (4th Filament)

%Set constants for background potential (Sydora et al. 2019 EQ:4.1)
C0=0;  %Set to '1' to use this potential or zero to negate
C1=-3.309;
C2=-.690;
C3=5.712;
C4=.397;
C5=-75.41;

%Background potential well (Sydora et al. 2019 EQ:4.1)
%Centerd at r=0
PHI_PW=C0*(C1+(C2.*exp(-C3.*(R-C4).^2))+(C5.*(R+3).^-4));

%Sydora et al. 2019 electric feild
Er_Sy=C0.*((2.*C3.*C2.*(-R+C4).*exp(-C3.*(R-C4).^2)-((4.*C5)./((R+3).^5))));

Ex_Sy=-C0.*(cos(TH).*Er_Sy);
Ey_Sy=-C0.*(sin(TH).*Er_Sy);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FEILDS

%Total potential 
PHI_1 = PHI_1t + PHI_3 + PHI_5 + PHI_7;
PHI_2 = PHI_2t + PHI_4 + PHI_6 + PHI_8;

% Flow Field
U =  ((real(Ey_1) + real(Ey_2) + real(Ey_3) + real(Ey_4)...
    + real(Ey_5) + real(Ey_6) + real(Ey_7) + real(Ey_8)...
    )/B*1e3) + (Ey_Sy/B); %cm/s

V = -((real(Ex_1) + real(Ex_2) + real(Ex_3) + real(Ex_4)...
    + real(Ex_5) + real(Ex_6) + real(Ex_7) + real(Ex_8)...
    )/B*1e3) + (Ex_Sy/B); %cm/s 

% Numerical Flow Field (checks out)
% [Vn,Un] = gradient(-real(PHI_1+PHI_2),x,y);
% U = Un/B*1e3;
% V = -Vn/B*1e3;

%% Solve Particle Motion

initC1 = [-0.2 0.2]; % Cartesian;
initC2 = initC1+[randn*1e-2 0]; % Cartesian;
initC1 = initC1+[randn*1e-2 0]; % Cartesian;
% options = odeset('NonNegative',1);
[T1,X1] = ode45(@(T,X) F(T,X,A,M,K,Alpha,W,B), t, initC1);
[T2,X2] = ode45(@(T,X) F(T,X,A,M,K,Alpha,W,B), t, initC2);

x1 = X1(:,1);
y1 = X1(:,2);
r1 = sqrt(x1.^2+y1.^2);
th1 = atan2(y1,x1);

x2 = X2(:,1);
y2 = X2(:,2);
r2 = sqrt(x2.^2+y2.^2);
th2 = atan2(y2,x2);

r1 = r1;
r2 = r2;

FFTr1 = log(abs(fft(r1)/length(r1)));
FFTr2 = log(abs(fft(r2)/length(r2)));

%% Generate Figure 1
fig1 = figure(2);
%fig1.Position = [50  100   1575   420*2];
fig1.Position = [50  100   875   420];
clf
fig1.Color = 'w';
ax1 = subplot(1,1,1);
hold on
axis square
colormap(modeturbo); % requires download of turbo map from GitHub
% colormap(turbo()); % requires download of turbo map from GitHub
% colormap(redblue()); % requires download of redblue map from GitHub
% colormap('jet');

% Plot the potential map
% S = pcolor(X,Y,real(PHI_1)); 
% S = pcolor(X,Y,real(PHI_2)); 
S = pcolor(X,Y,real(PHI_1 + PHI_2 + PHI_PW)); 
shading interp

% Add contour lines at specific phi;
% Phi_Contour = 0.0;
% [tmp,C] = contour3(X,Y,S.CData,[1,1]*Phi_Contour,'k-','LineWidth',2);

% Add ExB flow arrows;
%sc = 8;

%Q = quiver(X(1:sc:end,1:sc:end),Y(1:sc:end,1:sc:end),...
%    U(1:sc:end,1:sc:end),V(1:sc:end,1:sc:end),'color','k','LineWidth',1.5);

% Add particle trajectory
%L0 = plot(x1(1),y1(1),'ro','LineWidth',1,'MarkerSize',15,'MarkerFaceColor','r');
%L0 = plot(x2(1),y2(1),'bo','LineWidth',1,'MarkerSize',15,'MarkerFaceColor','b');
if make_movie ==0
%L1 = plot(x1,y1,'r-','LineWidth',3);
%L2 = plot(x2,y2,'b-','LineWidth',3);
else
%L1 = plot(x1(1:1),y1(1:1),'r-','LineWidth',3);
%L2 = plot(x2(1:1),y2(1:1),'b-','LineWidth',3);
end

% Subplot 2
%ax2 = subplot(1,3,2);
%hold on
%LR1 =% plot(t*1e3,r1,'k-','LineWidth',3);
%LRC1 = plot(t(1:1)*1e3,r1(1:1),'r-','LineWidth',3);
%LR2 = plot(t*1e3,r2,'k-','LineWidth',3);
%LRC2 = plot(t(1:1)*1e3,r2(1:1),'b-','LineWidth',3);
%if make_movie ==0
%LR1.Color = 'r';
%LR2.Color = 'b';
%end

% Subplot 3
%ax3 = subplot(1,3,3);
%hold on
%LFR1 = plot(freq,FFTr1,'r-','LineWidth',3);
%LFR2 = plot(freq,FFTr2,'b-','LineWidth',3);

% Figure Aesthetics
view(ax1,2)       
ax1.CLim = [0,1];
%ax1.CLim = [-5 -3.4];
%ax1.CLim= [-max([A1,A2])/2. max([A1,A2])/2.];
%ax1.CLim = [-.5 .5];
ax1.XLim = [min(x) max(x)];
ax1.YLim = [min(y) max(y)];
ax1.LineWidth = 2;
ax1.XMinorTick = 'on';
ax1.YMinorTick = 'on';
ax1.TickLength = ax1.TickLength*2;
ax1.TickDir = 'in';
ax1.TickLabelInterpreter = 'latex';
ax1.Layer = 'top';
ax1.Box = 'on';
ax1.FontSize = 18;
%ax1.Title.String  = strcat('$\omega_*t =$ ',sprintf(' %0.2f',f1*t(1)));
ax1.XLabel.String = '$x$ (cm)';
ax1.YLabel.String = '$y$ (cm)';
%ax1.Title.Interpreter = 'latex';
ax1.XLabel.Interpreter = 'latex';
ax1.YLabel.Interpreter = 'latex';

hc = colorbar(ax1);
hc.LineWidth = 2;
hc.TickLabelInterpreter = 'latex';
hc.FontSize = 18;
hc.Ticks = linspace(ax1.CLim(1),ax1.CLim(2),11);
hc.Title.String = '$\delta \phi$ (V)';
hc.Title.Interpreter = 'latex';

%ax2.XLim = [0 max(LR1.XData)];
% ax2.YLim = [min(LR1.YData)-.1 max(LR1.YData)+.1];
%ax2.LineWidth = 2;
%ax2.XMinorTick = 'on';
%ax2.YMinorTick = 'on';
%ax2.TickLength = ax2.TickLength*2;
%ax2.TickDir = 'in';
%ax2.TickLabelInterpreter = 'latex';
%ax2.Layer = 'top';
%ax2.Box = 'on';
%ax2.FontSize = 18;
%ax2.XLabel.String = '$t$ (ms)';
%ax2.YLabel.String = '$r$ (cm)';
%ax2.XLabel.Interpreter = 'latex';
%ax2.YLabel.Interpreter = 'latex';

%ax3.XLim = [0 1000];
% ax3.YLim = [min(LFR1.YData)-.25 max(LFR1.YData)+.25];
%ax3.LineWidth = 2;
%ax3.XMinorTick = 'on';
%ax3.YMinorTick = 'on';
%ax3.TickLength = ax3.TickLength*2;
%ax3.TickDir = 'in';
%ax3.TickLabelInterpreter = 'latex';
%ax3.Layer = 'top';
%ax3.Box = 'on';
%ax3.FontSize = 18;
%ax3.XLabel.String = '$f$ (kHz)';
%ax3.YLabel.String = 'Power';
%ax3.XLabel.Interpreter = 'latex';
%ax3.YLabel.Interpreter = 'latex';

%ax1.Units = 'normalized';
%ax2.Units = 'normalized';
%ax3.Units = 'normalized';
%ax1.Position = [.5 .1 .45 .8];
%ax2.Position = [.05 .55 .4 .35];
%ax3.Position = [.05 .1 .4 .35];

%% Update
if make_movie == 1;
for it = 1:speed:length(t);
% Update Potential
% S.CData = real(PHI_1.*exp(-1i.*w1.*t(it))); % m1 only
% S.CData = real(PHI_2.*exp(-1i.*w2.*t(it))); % m2 only      
S.CData = real(PHI_1.*exp(-1i.*w1.*t(it)) + PHI_2.*exp(-1i.*w2.*t(it))...
               +PHI_PW); 

% Update Contour
C.ZData = S.CData;

% Update Flows
%U =  real(Ey_1.*exp(-1i.*w1.*t(it))) + real(Ey_2.*exp(-1i.*w2.*t(it)))...
%    +real(Ey_3.*exp(-1i.*w3.*t(it))) + real(Ey_4.*exp(-1i.*w4.*t(it)))...
%    +real(Ey_5.*exp(-1i.*w5.*t(it))) + real(Ey_6.*exp(-1i.*w6.*t(it)))...
%    +real(Ey_7.*exp(-1i.*w7.*t(it))) + real(Ey_8.*exp(-1i.*w8.*t(it)))...
%    +real(Ey_Sy);

%V = -real(Ex_1.*exp(-1i.*w1.*t(it))) - real(Ex_2.*exp(-1i.*w2.*t(it)))...
%    -real(Ex_3.*exp(-1i.*w3.*t(it))) - real(Ex_4.*exp(-1i.*w4.*t(it)))...
%    -real(Ex_5.*exp(-1i.*w5.*t(it))) - real(Ex_6.*exp(-1i.*w6.*t(it)))...
%    -real(Ex_7.*exp(-1i.*w7.*t(it))) - real(Ex_8.*exp(-1i.*w8.*t(it)))...
%    -real(Ex_Sy);

%Q.UData = U(1:sc:end,1:sc:end);
%Q.VData = V(1:sc:end,1:sc:end);

% Update Particle
%L1.XData = x1(1:it);
%L1.YData = y1(1:it);
%L2.XData = x2(1:it);
%L2.YData = y2(1:it);

%LRC1.XData = t(1:it)*1e3;
%LRC1.YData = r1(1:it);
%LRC2.XData = t(1:it)*1e3;
%LRC2.YData = r2(1:it);
% Update Title
%ax1.Title.String  = strcat('$\omega_*t =$ ',sprintf(' %0.2f',f1*t(it)));

drawnow
end
end
%% Generate Figure 2
if make_fig2 == 1
% A1 = 7;
% A2 = 7;

phi_1 = A1.*besselj(m1,k1.*r).*exp(-alpha1.*r);
phi_2 = A2.*besselj(m2,k2.*r).*exp(-alpha2.*r);

fig2 = figure(2);
clf
fig2.Color = 'w';
ax1 = gca;
hold on
L0 = plot(r,phi_1,'k-.','LineWidth',4);
L1 = plot(r,phi_2,'k:','LineWidth',4);

ax1.XLim = [0,.4];
ax1.YLim = [min([phi_1,phi_2])*1.5, max([phi_1,phi_2])*1.5];
% ax.YLim = [-2,5];
ax1.LineWidth = 2;
ax1.XMinorTick = 'on';
ax1.YMinorTick = 'on';
ax1.TickLength = ax1.TickLength*4;
ax1.TickDir = 'in';
ax1.TickLabelInterpreter = 'latex';
ax1.Layer = 'top';
grid off
shading interp
box on
ax1.FontSize = 18;
% ax.Title.String  = strcat('$\omega_*t =$ ',sprintf(' %0.2f',f1*t(1)));
ax1.XLabel.String = '$r$ (cm)';
ax1.YLabel.String = '$\delta \phi$ (V)';
% ax.Title.Interpreter = 'latex';
ax1.XLabel.Interpreter = 'latex';
ax1.YLabel.Interpreter = 'latex';

leg = legend(ax1,'$\delta \phi_1$ ($m=1$)','$\delta \phi_6$ ($m=6$)');
leg.Interpreter = 'latex';
leg.FontSize = 18;
leg.Location = 'Northeast';
LP = leg.Position;
leg.Position = [LP(1)-.1, LP(2)-.1, LP(3), LP(4)];
leg.Box = 'off';
end

%% Flow Field 
function dX = F(T,X,A,M,K,Alpha,W,B)
    % Returns the velocity vectors at X according to DW potential.
    % T == time (s)
    % dX == velocity vector
    % X == position vector 
    % A == amplitude vector (V)
    % M == mode vector 
    % K == wavenumber vector (1/cm)
    % Alpha == decay vector (1/cm)
    % B == magnetic field (T)

    dX   = zeros(2,1); 
    TH = atan2(X(2),X(1));
    R = sqrt(X(1).^2+X(2).^2);
    % Cylindrical E-fields
    Er  = -exp(-1i.*W(1).*T).*A(1).*exp(M(1)*1i.*TH).*exp(-Alpha(1).*R).*(...
         besselj_p(M(1),R,K(1)) - Alpha(1).*besselj(M(1),K(1).*R))...
          -exp(-1i.*W(2).*T).*A(2).*exp(M(2)*1i.*TH).*exp(-Alpha(2).*R).*(...
         besselj_p(M(2),R,K(2)) - Alpha(2).*besselj(M(2),K(2).*R));
    Eth =  -1i.*M(1).*exp(-1i.*W(1).*T).*A(1).*exp(M(1).*1i.*TH).*...
            exp(-Alpha(1).*R).*besselj(M(1),K(1).*R)./R...       
           -1i.*M(2).*exp(-1i.*W(2).*T).*A(2).*exp(M(2).*1i.*TH).*...
           exp(-Alpha(2).*R).*besselj(M(2),K(2).*R)./R;
       
    Ex = cos(TH).*Er-sin(TH).*Eth;
    Ey = sin(TH).*Er+cos(TH).*Eth;
    
    % Velocity Fields
    dX(1) =  real(Ey)/B*1e3; % cm/s
    dX(2) = -real(Ex)/B*1e3;   % cm/s
    return
end

%% Bessel Derivative Function
function J_prime = besselj_p(nu,x,a)
J_prime = 0.5*a*(besselj(nu-1,a.*x)-besselj(nu+1,a.*x));
end



