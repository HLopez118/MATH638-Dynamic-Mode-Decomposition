s2=[0.1 1]; %error variance

%x-dynamics
for jj = 1:length(s2)   
    if jj == 1
        x_IC = x0+s2(jj).*randn(1,2); %perturbed initial condition
        [t,xsol]=ode45('vandpol',t,x_IC,[],mu);
        x=xsol(:,1);%projected x values
        plot(t,x,'r'), hold on
    else
        x_IC = x0+s2(jj).*randn(1,2); %perturbed initial condition
        [t,xsol1]=ode45('vandpol',t,x_IC,[],mu);
        x1=xsol1(:,1);%projected x values
        plot(t,x1,'b'), hold on
    end
end

figure(1)
plot(t,x_true,'k','Linewidth',1), hold off
xlabel t; ylabel x;
legend('\sigma_{2}=0.1', '\sigma_{2}=1','Analytical Solution','Interpreter',...
    'latex','location','best')
title('Van der Pol Equation with perturbed initial conditions - x solution')

%uncomment for y-dynamics
% %y-dynamics
% for j = 1:length(s2)   
%         x_IC = x0+s2(j).*randn(1,2); %perturbed initial condition
%         [t,xsol]=ode45('vandpol',t,x_IC,[],mu);
%         y=xsol(:,2);%projected y values
%         plot(t,y,'o'), hold on
% end
% 
% figure(2)
% plot(t,y_true,'k','Linewidth',3), hold off
% xlabel t; ylabel y;
% legend('\sigma_{2}=0.1', '\sigma_{2}=1','Analytical Solution','Interpreter',...
%     'latex','location','best')
% title('Van der Pol Equation with perturbed initial conditions - x solution')