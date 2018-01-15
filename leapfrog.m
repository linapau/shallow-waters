function[u, v, h] = leapfrog(u, v, h, fu, fv, tau, t)

% leapfrog scheme applied to the shallow water equations

dx = 25*10^3;
dy = 25*10^3;
dt = 3600*1;

imax=200; % extent in x direction
jmax=120; % extent in y direction

% shape of the basin
itop(1:80)=ones(1,80)*imax; 
itop(81:121)=ones(1,41)*120;

H = 50; % undisturbed thickness of the layer
gred = 9.8 * 2/1000; % reduced gravity

for j=1:jmax;
for i=2:itop(j)
u(i,j,t+1) = u(i,j,t-1) + (2*dt*fu(j)/4) * (v(i-1,j+1,t) + v(i,j+1,t) + v(i-1,j,t) + v(i,j,t)) - ...
    (gred*2*dt/dx) * (h(i,j,t) -h(i-1,j,t)) +  2*dt*tau(i,j)/(50*1000);
end
end


for j=2:jmax
for i=1:itop(j)
v(i,j,t+1) = v(i,j,t-1) - (2*dt*fv(j)/4) * (u(i,j,t) + u(i+1,j,t) + u(i,j-1,t) + u(i+1,j-1,t)) - ...
    (gred*2*dt/dy) * (h(i,j,t) -h(i,j-1,t));
end;
end;

for j=1:jmax
for i=1:itop(j)
h(i,j,t+1) = h(i,j,t-1) - (2*dt*H) * ( (u(i+1,j,t) - u(i,j,t))/dx + (v(i,j+1,t) - v(i,j,t))/dy);
end;
end
