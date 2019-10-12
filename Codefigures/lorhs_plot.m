clear all
clc

t=0:0.01:20;
s=10;
b=8/3;
r=28;

x0=[5 5 5];
[t,xsol]=ode45('lor_rhs',t,x0,[],s,b,r);

x_true=xsol(:,1); 
y_true=xsol(:,2); 
z_true=xsol(:,3);


figure(1)
rotate3d on
plot3(x_true,y_true,z_true,'b');
grid on
xlabel x; ylabel y; zlabel z

figure(2)
subplot(3,1,1)
plot(t,x_true,'b')
xlabel t; ylabel x;
subplot(3,1,2)
plot(t,y_true,'r')
xlabel t; ylabel y;
subplot(3,1,3)
plot(t,z_true,'k')
xlabel t; ylabel z;