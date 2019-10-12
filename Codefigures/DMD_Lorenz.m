% space
clc
clear all
%p = [sigma r b]
IC = [5; 5; 5];
p = [10 28 8/3];
tproj = 30;
tmax = 20;
dt = 0.1;
tspan = 0:dt:tmax;
tpspan = 0:dt:tproj;
N = length(tspan);
Np = length(tpspan);
M = length(tspan)+1;
opts.tol = 1e-8;
opts.maxit = 1500;

[Y,X,Z,B,rpos] = Lorenz_Data(IC,p,tmax,dt);

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
[t,Y] = ode45(@Lorenz_eq, tpspan, IC,[], p(1),p(2),p(3));
Y = Y';
%upos = find(u_dmd(:,1)>4,3);
for i = 1:3
    subplot(3,1,i)
    plot(tspan,X(rpos(i),:),'r')
    hold on
    plot(tpspan,u_dmd(rpos(i),:),'b--')
    plot(tpspan,Y(i,:),'g')
    ylim([-50,50])
    grid on
    hold off
end
figure()
% SSE = sum((X-u_dmd).^2);
% 
% plot(tspan,SSE)

xx = -1:0.01:1;
c1 = sqrt(1-xx.^2);
c2 = -sqrt(1-xx.^2);
plot(xx,c1,'r',xx,c2,'r')
hold on
plot(mu,'o')
grid on



% for i = 1:length(t)
%     plot(x,usol(i,:),'b-o',x,tu_dmd(i,:),'r-');
%     ylim([-4,4])
%     xlim([-20,20])
%     title('Shrodinger eq.')
%     xlabel('x')
%     ylabel('u')
%     grid on
%     M(i) = getframe;
% end

