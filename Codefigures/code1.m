clc
clear all
% space
L=40; n=512;
x2=linspace(-L/2,L/2,n+1); x=x2(1:n);
k=(2*pi/L)*[0:n/2-1 -n/2:-1].';
% time
slices=100;
tmin = 6*pi;
t=linspace(0,tmin,slices+1); dt=t(2)-t(1);
tmax = 16*pi;
xslices = round((tmax)/dt);
t1=linspace(0,tmax,xslices+1); 


u=2*sech(x).'; % initial conditions
ut=fft(u);
[t2,utsol]=ode45('nls_rhs',t1,ut,[],k);
for j=1:length(t2)
usol(j,:)=ifft(utsol(j,:)); % bring back to space
end

tsol = (length(usol(:,1)) - length(t));

X = usol(1:end-tsol,:).'; X1 = X(:,1:end-1); X2 = X(:,2:end);

[U,Sigma,V] = svd(X1, 'econ');
S = U'*X2*V*diag(1./diag(Sigma));
[eV,D] = eig(S);
mu = diag(D);
omega = log(mu)/(dt);
Phi = U*eV;

y0 = Phi\u; % pseudo-inverse initial conditions
u_modes = zeros(size(V,2),length(t));
for iter = 1:length(t1) %this is extended to t1 in order to calculate predictions. 
u_modes(:,iter) =(y0.*exp(omega*t1(iter)));
end
u_dmd = Phi*u_modes;
tu_dmd = u_dmd';


for i = length(t):length(t1) %plot of just predicted values
    plot(x,usol(i,:),'b-o',x,tu_dmd(i,:),'r-');
    ylim([-4,4])
    xlim([-20,20])
    title(['Shrodinger eq.', num2str(t1(i))])
    xlabel('x')
    ylabel('u')
    grid on
    M(i) = getframe;
end

