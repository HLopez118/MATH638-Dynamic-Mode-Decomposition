%noisy observation every t=0.5
ntd = 1;
tdata=t(1:ntd:end);
n=length(tdata);

xn=randn(n,1);
yn=randn(n,1);

x=xsol(:,1); 
y=xsol(:,2); 
x1=xsol1(:,1); 
y1=xsol1(:,2); 

s3=[0.1 1]; %error variance in data

%perturbed observations
xdata=x(1:ntd:end)+s3(:,1).*xn;
ydata=y(1:ntd:end)+s3(:,1).*yn;
x2data=x1(1:ntd:end)+s3(:,2).*xn;
y2data=y1(1:ntd:end)+s3(:,2).*yn;


[t,xsol]=ode45('vandpol',t,x0,[],mu); %unperturbed IC

%plot the x-dynamics
figure(1)
subplot(2,1,1)
plot(t,x_true,'b'), hold on
plot(tdata,xdata,'r--'),hold off
pbaspect([2 1 1])%makes x-axis twice the y-axis
xlabel t; ylabel x;
title('Perturbed observations for $\sigma_{3}=0.01$','interpreter','latex')
legend('Analytical Solution','Perturbed Observations','location','BestOutside')
subplot(2,1,2)
plot(t,x_true,'b'), hold on
plot(tdata,x2data,'r--'),hold off
pbaspect([2 1 1]) %makes x-axis twice the y-axis
xlabel t; ylabel x;
title('Perturbed observations for $\sigma_{3}=1$','interpreter','latex')
legend('Analytical Solution','Perturbed Observations','location','BestOutside')

%plot the y-dynamics
figure(2)
subplot(2,1,1)
plot(t,y_true,'b'), hold on
plot(tdata,ydata,'r--'),hold off
% pbaspect([2 1 1])%makes x-axis twice the y-axis
xlabel t; ylabel y;
title('Perturbed observations for $\sigma_{3}=0.01$','interpreter','latex')
legend('Analytical Solution','Perturbed Observations','location','BestOutside')
subplot(2,1,2)
plot(t,y_true,'b'), hold on
plot(tdata,y2data,'r--'),hold off
% pbaspect([2 1 1]) %makes x-axis twice the y-axis
xlabel t; ylabel y;
title('Perturbed observations for $\sigma_{3}=1$','interpreter','latex')
legend('Analytical Solution','Perturbed Observations','location','BestOutside')