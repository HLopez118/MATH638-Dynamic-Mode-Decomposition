clc
clear all
IC = [1; 1];
mu = 3;
tcolm = 10;
trowm = 15;
dt = 0.01;
tproj = tcolm+trowm;
tpspan = 0:dt:tproj;
Np = length(tpspan);
opts.tol = 1e-8;
opts.maxit = 1500;


tcol = 0:dt:tcolm;
trow = 0:dt:trowm;
p = length(tcol);
q = length(trow);

pmax = tcolm+2*trowm;
pspan = 0:dt:pmax;

[t,YY] = ode45(@VdP_eq, pspan, IC,[], mu);
YY = YY';

[Y,rX,X,rY] = VdP_Hankel(IC,mu,tcolm,trowm,dt);

X1 = X(:,1:end-1); X2 = X(:,2:end);

r = rank(X1);

% [U,Sigma,V] = svd(X1,'econ');
[U,Sigma,V] = lmsvd(X1,100,opts);
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

usol_dmd = [u_dmd(1,1:p-1) u_dmd(p,1:(q-p)) u_dmd(q,:)];
trusol = usol_dmd(1:2001);
rData = [X(1,1:p-1) X(p,1:(q-p)) X(q,:)];
trData = rData(1:2001);
plot(pspan,YY(1,:),'b','linewidth',3)
hold on
plot(pspan(1:length(rData)),rData,'k')
plot(pspan,usol_dmd,'r','linewidth',1.5)
xlabel t; ylabel x;
legend('Analytical Solution','Data','DMD Reconstruction/Prediction')
% legend('Analytical Solution','DMD Reconstruction/Prediction')
title('Dynamic Mode Decomposition')
grid on
hold off

figure()
plot(pspan(1:2001),YY(1,(1:2001)),'b','linewidth',3)
hold on
plot(pspan(1:2001),trData,'k')
plot(pspan(1:2001), trusol,'r','linewidth',1.5)
xlabel t; ylabel x;
legend('Analytical Solution','Data','DMD Reconstruction')
% legend('Analytical Solution','DMD Reconstruction')
title('Dynamic Mode Decomposition')
grid on





figure()
xx = -1:0.01:1;
c1 = sqrt(1-xx.^2);
c2 = -sqrt(1-xx.^2);
plot(xx,c1,'r',xx,c2,'r')
hold on
plot(mu,'o')
title('Spectral Distribution of Eigenvalues, mu_{k}')
xlabel('Re(\mu)') 
ylabel('IM(\mu)')
grid on


figure()
SE = (YY(1,1:2001)-trusol).^2;

plot(pspan(1:2001),SE)

nrmd = norm(YY(1,1:2001)-trData); %r = 100 == 44.6230;
nrmx = norm(YY(1,1:2001)-trusol); %r = 100 == 17.2918;
SSE = sum(YY(1,1:2001)-trusol); %33.7382

% usol_dmd = [u_dmd(1,1:p-1) u_dmd(p,1:(q-p)) u_dmd(q,:)];

% plot(pspan,YY(1,:),'g')
% hold on
% plot(tpspan(1:p-1),u_dmd(1,1:p-1),'b--')
% plot(tpspan(p:q-1),u_dmd(p,1:(q-p)),'b--')
% plot(pspan(q:length(pspan)),u_dmd(q,:),'r--')
% % plot(pspan,usol_dmd,'m')
% ylim([-10,10])
% grid on
% hold off
% 
% 
% 
% figure()
% xx = -1:0.01:1;
% c1 = sqrt(1-xx.^2);
% c2 = -sqrt(1-xx.^2);
% plot(xx,c1,'r',xx,c2,'r')
% hold on
% plot(mu,'o')
% grid on
% 
% 
% 
% figure()
% SE = (YY(1,:)-usol_dmd).^2;
% 
% plot(pspan,SE)
% 
% nrm = norm(YY(1,:)-usol_dmd)
% SSE = sum((YY(1,:)-usol_dmd).^2)




