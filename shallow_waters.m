close all
clear all
 
%% PURPOSE: Replicate the results obtained in:
%
% The Seasonal Upwelling in the Gulf of Guinea Due to Remote Forcing
% Journal of Physical Oceanography, Vol. 8, 1050-1060
% 1978
% by David Adamec and James J. O'Brien
%
% for the course Numerical Modeling of the Atmosphere of the bachelor in
% Atmospheric Sciences - University of the Republic (Uruguay).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% It consists of a linear model on an equatorial beta plane integrated over
% a 100-days period in a basin that approximates the tropical Atlantic
% Ocean.
%
% The motion is described by the shallow water equations, discretized in a
% staggered grid in space (type C). For the time derivative a leapfrog 
% scheme (leapfrog.m) is considered. Every 30 iterations, the Matsuno 
% scheme (matsuno.m) is used. 
%
%Figuresare made to show the results.
%
% A stability analysis is made
%
% This script calls the functions leapfrog.m and matsuno.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1) Settings

% Grid
m = 120; % latitude
n = 200; % longitude
t = 2400; % time

dx = 25*10^3; % grid resolution of 25 km in x direction
dy = 25*10^3; % grid resolution of 25 km in y direction
dt = 3600*1; % time step of 1 hour

h = zeros(n,   m,   t); % local height field change from undisturbed depth in a layer
u = zeros(n+1, m,   t); % x directed component of velocity
v = zeros(n,   m+1, t); % y directed component of velocity

%% Remote forcing

Tau = zeros(n+1, m); % wind stress

% x dependence of wind stress
Tau(1:60, :) = -0.025; % constant on the west coast

for i = 61:90
    Tau(i, :) = ( sin( (i-60)*pi/(30*2) )^2 - 1 ) * 0.025;
end

figure
plot(Tau(:,1))

%% Beta and Coriolis parameter

beta=2.11e-11;

fv=(-1500:25:1500)*1e3*beta; % Coriolis parameter (meridional component)

fu=(fv(1:120)+fv(2:121))/2; % Coriolis parameter (zonal component)

%% 2) Simulation of u, v and h in the stable case

for kk=1:80, 
    t=1+(kk-1)*30; % 100-days period of integration
    [u, v, h] = matsuno(u, v, h, fu, fv, Tau, t); % Matsuno scheme
    for t=1+(kk-1)*30+1:1+(kk-1)*30+29,
        [u, v, h] = leapfrog(u, v, h, fu, fv, Tau, t); % Leapfrog scheme
    end
end

%% 3) Figures of the results

% labels
v = -50:10:50;
v2 =-0.60:0.1:0.60;

% Contours of h day 10 (t=240)
figure
hold on
[c,hc] = contour(h(:,:,240)', 5:5:50,'ShowText','on',...
    'color', 'k','linewidth', 1.7); % positive values
[cneg,hcneg] = contour(h(:,:,240)','LevelList', -50:5:-5,'ShowText','on',...
    'color','k','linewidth', 1.7, 'lineStyle', ':'); % negative values
clabel(c,hc,v,'FontSize',15)
clabel(cneg,hcneg,v,'FontSize',15)
set(gca,'XTick',0:40:200)
set(gca,'XTickLabel',['0';'1';'2';'3';'4';'5'],'FontSize',15)
set(gca,'YTick',0:40:120)
set(gca,'YTickLabel',['0';'1';'2';'3'],'FontSize',15)
xlabel('km*10^3','FontSize',15)
ylabel('km*10^3','FontSize',15)
caxis([-65 50])
title('h(day 10)','FontSize',15)


% Contours of h day 50 (t=1200)
figure
hold on
[c,hc] = contour(h(:,:,1200)', 5:5:50,'ShowText','on',...
    'color', 'k','linewidth', 1.7); % positive values 
[cneg,hcneg] = contour(h(:,:,1200)','LevelList', -50:5:-5,'ShowText','on', ...
'color', 'k','linewidth', 1.7, 'lineStyle', ':'); % negative values
clabel(c,hc,v,'FontSize',15)
clabel(cneg,hcneg,v,'FontSize',15)
set(gca,'XTick',0:40:200)
set(gca,'XTickLabel',['0';'1';'2';'3';'4';'5'],'FontSize',15)
set(gca,'YTick',0:40:120)
set(gca,'YTickLabel',['0';'1';'2';'3'],'FontSize',15)
xlabel('km*10^3','FontSize',15)
ylabel('km*10^3','FontSize',15)
caxis([-65 50])
title('h(day 50)','FontSize',15)


% Contours of h day 100 (t=2400)
figure
hold on
[c,hc] = contour(h(:,:,2400)', 5:5:50,'ShowText','on',...
    'color', 'k','linewidth', 1.7); % positive values
[cneg,hcneg] = contour(h(:,:,2400)','LevelList', -50:5:-5,'ShowText','on', ...
'color', 'k','linewidth', 1.7, 'lineStyle', ':'); % negative values
clabel(c,hc,v,'FontSize',15)
clabel(cneg,hcneg,v,'FontSize',15)
set(gca,'XTick',0:40:200)
set(gca,'XTickLabel',['0';'1';'2';'3';'4';'5'],'FontSize',15)
set(gca,'YTick',0:40:120)
set(gca,'YTickLabel',['0';'1';'2';'3'],'FontSize',15)
xlabel('km*10^3','FontSize',15)
ylabel('km*10^3','FontSize',15)
caxis([-65 50])
title('h(d�a 100)','FontSize',15)


% Hovmoeller diagram (time-dependent contours) of the u component of the 
% flow 1000 km from the western boundary running 500 km at either side of
% the equator
figure
hold on
[c,hc] = contour(squeeze(u(40,40:80,:)), 0.05:0.05:0.15,'ShowText','on',...
    'color', 'k','linewidth', 1.7); % positive velues
[cneg,hcneg] = contour(squeeze(u(40,40:80,:)),'LevelList', -0.6:0.05:-0.05,'ShowText','on', ...
'color', 'k','linewidth', 1.7, 'lineStyle', ':'); % negative values
clabel(c,hc,v2,'FontSize',15)
clabel(cneg,hcneg,v2,'FontSize',15)
set(gca,'XTick',0:240:2400)
set(gca,'XTickLabel',{'0';'10';'20';'30';'40';'50';'60';'70';'80';'90';'100'},'FontSize',15)
set(gca,'YTick',0:20:40)
set(gca,'YTickLabel',{'1';'1.5';'2'},'FontSize',15)
ylabel('km*10^3','FontSize',15)
xlabel('D�as','FontSize',15)
title('u 1000 km from the western boundary and 500 km either side of the equator','FontSize',15)
caxis([-0.5 0.3])

%% 3) Maximum and minimum values of h

min_d10 = min(min(h(:,:,240))); % the minimum value of h day 10
max_d10 = max(max(h(:,:,240))); % the maximum value of h day 10

min_d50 = min(min(h(:,:,1200))); % the minimum value of h day 50
max_d50 = max(max(h(:,:,1200))); % the maximum value of h day 10

min_d100 = min(min(h(:,:,2400))); % the minimum value of h day 100
max_d100 = max(max(h(:,:,2400))); % % the maximum value of h day 100


%% 4) Stability analysis of the leapfrog scheme

% Settings
d = 50000;
y = 0; % y can be zero or the -1500*1e3 (southern boundary)
f = beta*y; % Coriolis parameter
H = 50; % undisturbed thickness of the layer
g = 9.8*0.001; % reduced gravity
vel_ref = sqrt(9.8*0.001*H); % speed of the reference wave 

lambda = NaN(10, 20); % maximum eigenvalues of the Von Neumann matrix

for i = 1:10
    k = (2*pi)/(i*d); % wavenumber
    j=1;
    for c = 0.05:0.05:1 % Courant number
        dt = c*d/vel_ref; % time
        % define Von Neumann matrix (VN)
        alpha = 2*dt*f*cos(k*d/2);
        gamma = 4*dt*1i*sin(k*d/2)/d;
        VN = [0, alpha, -g*gamma, 1, 0, 0;
            -alpha, 0, 0, 0, 1, 0;
            -H*gamma, 0, 0, 0, 0, 1;
            1, 0, 0, 0, 0, 0;
            0, 1, 0, 0, 0, 0;
            0, 0, 1, 0, 0, 0];
        % eigenvalues of VN
        emax = max(abs(eig(VN)));
        lambda(i,j) = emax;
        j = j+1;
    end
end

% the stability condition is lambda < 1


