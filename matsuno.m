function[u, v, h] = matsuno(u, v, h, fu, fv, tau, t)

% Matsuno scheme applied to the shallow water equations

dx = 25*10^3;
dy = 25*10^3;
dt = 3600*1;

imax=200; % extent in x direction
jmax=120; % extent in y direction

% shape of the basin
itope(1:80)=ones(1,80)*imax;
itope(81:121)=ones(1,41)*120;

H = 50; % undisturbed thickness of the layer
gred = 9.8 * 2/1000; % reduced gravity


hfg = zeros(imax,   jmax);
ufg = zeros(imax+1, jmax);
vfg = zeros(imax,   jmax+1);


% first guess for instant t+1

for j=1:jmax;
for i=2:itope(j)
ufg(i,j) = u(i,j,t) + (dt*fu(j)/4) * (v(i-1,j+1,t) + v(i,j+1,t) + v(i-1,j,t) + v(i,j,t)) - ...
    (gred*dt/dx) * (h(i,j,t) -h(i-1,j,t)) + dt*tau(i,j)/(50*1000);
end
end

for j=2:jmax
for i=1:itope(j)
vfg(i,j) = v(i,j,t) - (dt*fv(j)/4) * (u(i,j,t) + u(i+1,j,t) + u(i,j-1,t) + u(i+1,j-1,t)) - ...
    (gred*dt/dy) * (h(i,j,t) -h(i,j-1,t));
end;
end;

for j=1:jmax
for i=1:itope(j)
hfg(i,j) = h(i,j,t) - (dt*H) * ( (u(i+1,j,t) - u(i,j,t))/dx + (v(i,j+1,t) - v(i,j,t))/dy);
end;
end


% second approximation using first guess to evaluate the operator of the 
% differential equations

for j=1:jmax;
for i=2:itope(j)
u(i,j,t+1) = u(i,j,t) + (dt*fu(j)/4) * (vfg(i-1,j+1) + vfg(i,j+1) + vfg(i-1,j) + vfg(i,j)) - ...
    (gred*dt/dx) * (hfg(i,j) -hfg(i-1,j))  + dt*tau(i,j)/(50*1000);
end
end

for j=2:jmax
for i=1:itope(j)
v(i,j,t+1) = v(i,j,t) - (dt*fv(j)/4) * (ufg(i,j) + ufg(i+1,j) + ufg(i,j-1) + ufg(i+1,j-1)) - ...
    (gred*dt/dy) * (hfg(i,j) -hfg(i,j-1));
end;
end;

for j=1:jmax
for i=1:itope(j)
h(i,j,t+1) = h(i,j,t) - (dt*H) * ( (ufg(i+1,j) - ufg(i,j))/dx + (vfg(i,j+1) - vfg(i,j))/dy);
end;
end

