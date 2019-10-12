function [data, rdata, lrdata, ardata, rpos] = Lorenz_Data(IC,p,tmax,tstep)
tspan = 0:tstep:tmax;
N = length(tspan);
M = length(tspan)+1;
rpos = zeros(1,3);
noise = 2;

[t,Y] = ode45(@Lorenz_eq, tspan, IC,[], p(1),p(2),p(3));
Y = Y';
data = zeros(M,N);
ardata = data + noise*(rand(M,N)-rand(M,N));
lrdata = data + noise*(rand(M,N)-rand(M,N));


for k = 0:2
    m = mod(k,3);
    rpos(k+1) = round(M/2) + k;
    data(rpos(k+1),:) = Y(m+1,:);
    lrdata(rpos(k+1),:) = Y(m+1,:);
end

rdata = data + noise*(rand(M,N)-rand(M,N));
end

