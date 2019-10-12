clear all
clc

t=0:0.01:20;
mu=3;

x0=[1 1];
[t,xsol]=ode45('vandpol',t,x0,[],mu);

x_true=xsol(:,1); 
y_true=xsol(:,2); 


figure(1)
subplot(2,1,1)
grid on
plot(t,x_true,'r'), hold on

xlabel('t')
ylabel('X')
title('Solution to the Vanderpol Equation for $\mu=3$','Interpreter','Latex')

subplot(2,1,2)
plot(t,y_true,'b'), 
xlabel('t')
ylabel('Y')
hold off