%s2=s3=0.01
x_da=[]; %data assimilation solution

for j=1:length(tdata)-1 %step through every t=0.5;
    tspan=0:0.001:0.01; %time b/w data collection
    [tspan,xsol]=ode45('vandpol',tspan,x_IC,[],mu);
    
    x_IC0=[xsol(end,1); xsol(end,2)]; %model estimate
    x_dat=[xdata(j+1); ydata(j+1)]; %data estimate
    K=s2(:,1)./(s2(:,1)+s3(:,1)); %Kalman gain %s2=s3=0.01
    x_IC=x_IC0+K*(x_dat-x_IC0); %adjusted state vector
    
    x_da=[x_da; xsol(1:end-1,:)]; %store the data

end
x_dax=[x_da; xsol(end,:)]; %store last data time
tnew = 0:0.001:20;
% clf;
% figure(1)
% plot(t,x,'b'), hold on
% plot(tdata,xdata,'k'), hold on
% x_dax=x_da(:,1);
% plot(tnew,x_dax,'r'),hold off
% xlabel t; ylabel x;
figure(1)
plot(t,x,'b','linewidth',3), hold on
plot(tdata,xdata,'k'), hold on
plot(tnew,x_dax(:,1),'r','linewidth',1.5)
xlabel t; ylabel x;
legend('Analytical Solution','Data','Data Assimilation Solution')
title('Data Assimilation')
grid on

figure()


%error b/w solution and model prediction
%error b/w solution and data assimilation technique
for j = 1:length(t)
    Edat(j,:) = abs((x(j,:)-x_dax((10*(j-1)+1),:)));
    Ervec1(j,:) = (x(j,:)-x_dax((10*(j-1)+1),:));
    errdat(j) = (sum(Edat(j,:))).^2;
end

nrmx = norm(Ervec1); %6.1928

figure(2)
plot(t,errdat,'b-')
pbaspect([2 1 1])
xlabel('Time')
ylabel('SSE')
title('Error Dynamics for X solution - Data Assimilation')
grid on


figure(3)
plot(t,y,'b'), hold on
plot(tdata,ydata,'k'), hold on
x_day=x_da(:,2);
plot(tnew,x_day,'r'),hold off
pbaspect([2 1 1]) %makes x-axis twice the y-axis
xlabel t; ylabel y;

for j = 1:length(t)
    Erdat(j,:) = abs((y(j,:)-x_day((10*(j-1)+1),:)));
    Ervec2(j,:) = (y(j,:)-x_day((10*(j-1)+1),:));
    err_dat(j) = (sum(Erdat(j,:))).^2;
end

nrmy = norm(Ervec2); %12.8875

figure(4)
plot(t,err_dat,'b-')
pbaspect([2 1 1])
xlabel('Time')
ylabel('SSE')
title('Error Dynamics for Y solution - Data Assimilation')
grid on