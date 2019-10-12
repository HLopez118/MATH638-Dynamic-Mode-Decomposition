function [H1,H2,H3,rH1,rH2,rH3] = Lorenz_Hankel(IC,par,tcolm,trowm,tstep)

tcol = 0:tstep:tcolm;
trow = 0:tstep:trowm;
m = tcolm+trowm;

per = 4;

tspan = 0:tstep:m;
p = length(tcol);
q = length(trow);

[t,Y] = ode45(@Lorenz_eq, tspan, IC,[], par(1),par(2),par(3));
Y = Y';

[r,c] = size(Y);
rY = per*(randn(r,c)) + Y;

H1 = zeros(q,p);
H2 = zeros(q,p);
H3 = zeros(q,p);
rH1 = zeros(q,p);
rH2 = zeros(q,p);
rH3 = zeros(q,p);

for i = 1:q
    H1(i,:) = Y(1,i:(p+(i-1)));
end

for ii = 1:q
    rH1(ii,:) = rY(1,ii:(p+(ii-1)));
end

for j = 1:q
    H2(j,:) = Y(2,j:(p+(j-1)));
end

for jj = 1:q
    rH2(jj,:) = rY(2,jj:(p+(jj-1)));
end

for k = 1:q
    H3(k,:) = Y(3,k:(p+(k-1)));
end

for kk = 1:q
    rH3(kk,:) = rY(3,kk:(p+(kk-1)));
end
end