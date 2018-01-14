close all
clear all

m = 120;
n = 200;
t = 2400;

dx = 25*10^3;
dy = 25*10^3;
% dt = 3600*3; % 
dt = 3600*1;

h = zeros(n,   m,   t);
u = zeros(n+1, m,   t);
v = zeros(n,   m+1, t);

%% Matriz Tau

Tau = zeros(n+1, m);

Tau(1:60, :) = -0.025;

for i = 61:90
    Tau(i, :) = ( sin( (i-60)*pi/(30*2) )^2 - 1 ) * 0.025;
end

beta=2.11e-11;

fv=(-1500:25:1500)*1e3*beta;

fu=(fv(1:120)+fv(2:121))/2;



figure
plot(Tau(:,1))

%% 
% 
% for kk=1:80, %80
%     t=1+(kk-1)*30;
%     [u, v, h] = matsuno(u, v, h, fu, fv, Tau, t); t,
%     for t=1+(kk-1)*30+1:1+(kk-1)*30+29,
%         [u, v, h] = leapfrog(u, v, h, fu, fv, Tau, t); t,
%     end
% end

%%
[u, v, h] = matsuno(u, v, h, fu, fv, Tau, 1);
for t=2:10;
    [u, v, h] = leapfrogSinDy(u, v, h, fu, fv, Tau, t); 
   
end