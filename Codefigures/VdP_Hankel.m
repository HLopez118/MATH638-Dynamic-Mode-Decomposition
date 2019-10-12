function [H1,H2,rH1,rH2] = VdP_Hankel(IC,mu,tcolm,trowm,tstep)

tcol = 0:tstep:tcolm;
trow = 0:tstep:trowm;
m = tcolm+trowm;

per = .1;

tspan = 0:tstep:m;
p = length(tcol);
q = length(trow);

[t,Y] = ode45(@VdP_eq, tspan, IC,[], mu);
Y = Y';

[r,c] = size(Y);
rY = per*randn(r,c) + Y;
% rY = Y;

H1 = zeros(q,p);
H2 = zeros(q,p);
rH1 = zeros(q,p);
rH2 = zeros(q,p);

for i = 1:q
    H1(i,:) = Y(1,i:(p+(i-1)));
end

for ii = 1:q
    H2(ii,:) = Y(2,ii:(p+(ii-1)));
end

for j = 1:q
    rH1(j,:) = rY(1,j:(p+(j-1)));
end

for jj = 1:q
    rH2(jj,:) = rY(2,jj:(p+(jj-1)));
end
end