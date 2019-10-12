%s2=s3=1
x_da=[]; %data assimilation solution

for j=1:length(tdata)-1 %step through every t=0.5;
    tspan=0:0.01:0.5; %time b/w data collection
    [tspan,xsol]=ode45('vandpol',tspan,x_IC,[],mu);
    
    x_IC0=[xsol(end,1); xsol(end,2)]; %model estimate
    x_dat=[xdata(j+1); ydata(j+1)]; %data estimate
    K=s2(:,2)./(s2(:,2)+s3(:,2)); %Kalman gain %s2=s3=1
    x_IC=x_IC0+K*(x_dat-x_IC0); %adjusted state vector
    
    x_da=[x_da; xsol(1:end-1,:)]; %store the data

end
x_da=[x_da; xsol(end,:)]; %store last data time

clf;
figure(1)
plot(t,x_true,'b'), hold on
plot(tdata,xdata,'ko'), hold on
x_dax=x_da(:,1);
plot(t,x_dax,'r'),hold off
xlabel t; ylabel x;


%error b/w solution and model prediction
%error b/w solution and data assimilation technique
for j = 1:length(t)
    Edat(j,:) = abs((x_true(j,:)-x_dax(j,:)));
    errdat(j) = (sum(Edat(j,:))).^2;
end

figure(2)
plot(t,errdat,'b-')
pbaspect([2 1 1])
xlabel('Time')
ylabel('SSE')
title('Error Dynamics for X solution - Data Assimilation')
grid on


figure(3)
plot(t,y_true,'b'), hold on
plot(tdata,ydata,'ko'), hold on
x_day=x_da(:,2);
plot(t,x_day,'r'),hold off
pbaspect([2 1 1]) %makes x-axis twice the y-axis
xlabel t; ylabel y;

for j = 1:length(t)
    Erdat(j,:) = abs((y_true(j,:)-x_day(j,:)));
    err_dat(j) = (sum(Erdat(j,:))).^2;
end

figure(4)
plot(t,err_dat,'b-')
pbaspect([2 1 1])
xlabel('Time')
ylabel('SSE')
title('Error Dynamics for Y solution - Data Assimilation')
grid on