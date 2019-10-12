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
opts.tol = 1e-8; %these options are for the truncated SVD function
opts.maxit = 1500;


tcol = 0:dt:tcolm;
trow = 0:dt:trowm;
p = length(tcol);
q = length(trow);

pmax = tcolm+2*trowm;
pspan = 0:dt:pmax;

[t,Y1] = ode45(@Lorenz_eq, pspan, IC,[], par(1),par(2),par(3));
Y1 = Y1';

[Y,Z,rX,X,rY,rZ] = Lorenz_Hankel(IC,par,tcolm,trowm,dt); % move 'X' 
% to the spot of the matrix you want. [x,y,z,rx,ry,rz]. r = matrix with noise
% corresponding to which lorenz solution you want.

X1 = X(:,1:end-1); X2 = X(:,2:end);

% [U,Sigma,V] = svd(X1,'econ'); %normal matlab SVD function.
[U,Sigma,V] = lmsvd(X1,100,opts); %SVD truncated at 100 eigenvalues. 
%comment out either depending on what you want to use. 
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
u_dmd = real(Phi*u_modes); %recreated matrix and projected solutions

usol_dmd = [u_dmd(1,1:p-1) u_dmd(p,1:(q-p)) u_dmd(q,:)];
trusol = usol_dmd(1:2001);
rData = [X(1,1:p-1) X(p,1:(q-p)) X(q,:)];
trData = rData(1:2001);
plot(pspan,Y1(1,:),'b','linewidth',3)
hold on
plot(pspan(1:length(rData)),rData,'ko')
plot(pspan,usol_dmd,'r','linewidth',1.5)
xlabel t; ylabel x;
legend('Analytical Solution','Data','DMD Reconstruction/Prediction')
% legend('Analytical Solution','DMD Reconstruction/Prediction')
title('Dynamic Mode Decomposition')
grid on
hold off

figure()
plot(pspan(1:2001),Y1(1,(1:2001)),'b','linewidth',3)
hold on
plot(pspan(1:2001),trData,'ko')
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
xlabel('Re($\mu$)') 
ylabel('IM($\mu$)')
grid on


figure()
SE = (Y1(1,(1:2001))-usol_dmd(1:2001)).^2;

plot(pspan(1:2001),SE)

nrmd = norm(Y1(1,(1:2001))-trData); %r = 100 == 179.8262;
nrmx = norm(Y1(1,(1:2001))-trusol); %r = 100 == 77.9739;
SSE = sum(Y1(1,(1:2001))-trusol); %33.7382



