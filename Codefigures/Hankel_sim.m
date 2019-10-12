% space
clc
clear all
%par = [sigma r b]
IC = [5; 5; 5];

par = [10 28 8/3];
tcolm = 10;
trowm = 15;
dt = 0.01;
tproj = tcolm+trowm;
tpspan = 0:dt:tproj;
Np = length(tpspan);
% opts.tol = 1e-8;
% opts.maxit = 1500;


tcol = 0:dt:tcolm;
trow = 0:dt:trowm;
p = length(tcol);
q = length(trow);

pmax = tcolm+2*trowm;
pspan = 0:dt:pmax;

[t,Y1] = ode45(@Lorenz_eq, pspan, IC,[], par(1),par(2),par(3));
Y1 = Y1';

[X,Y,Z] = Lorenz_Hankel(IC,par,tcolm,trowm,dt);

X1 = X(:,1:end-1); X2 = X(:,2:end);

[U,Sigma,V] = svd(X1,'econ');
% [U,Sigma,V] = lmsvd(X1,12,opts);
S = U'*X2*V*diag(1./diag(Sigma));
[eV,D] = eig(S);
mu = diag(D);
omega = log(mu)/(dt);
Phi = U*eV;

y0 = Phi\X(:,1); % pseudo-inverse initial conditions
u_modes = zeros(size(V,2),Np);
    for iter = 1:Np
    u_modes(:,iter) =(y0.*exp(omega*tpspan(iter)));
    end
u_dmd = real(Phi*u_modes);


plot(pspan,Y1(1,:),'g')
hold on
plot(tpspan(1:p),u_dmd(1,1:p),'b--')
plot(tpspan(p:q),u_dmd(p,1:(q-p+1)),'b--')
plot(pspan(q:length(pspan)),u_dmd(q,:),'r--')
ylim([-50,50])
grid on
hold off



figure()
xx = -1:0.01:1;
c1 = sqrt(1-xx.^2);
c2 = -sqrt(1-xx.^2);
plot(xx,c1,'r',xx,c2,'r')
hold on
plot(mu,'o')
grid on
% figure()
% SSE = sum((X-u_dmd).^2);
% 
% plot(tspan,SSE)



