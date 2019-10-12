s2=1; %error variance


x_IC = x0+s2*randn(1,3); %perturbed initial condition
[t,xsol]=ode45('lor_rhs',t,x_IC,[],s,b,r);
x=xsol(:,1); %projected x values


figure(1)
pbaspect([2 1 1])
plot(t,x,'r-'), hold on
plot(t,x_true,'b-'), hold off
xlabel t; ylabel x;
legend('Perturbation','Unperturbed IC')
title('Perturbed Initial Conditions')

for j = 1:length(t)
    Edat(j,:) = abs((x_true-x));
    errdat(j) = (sum(Edat(j,:))).^2;
end

figure(2)
plot(t,errdat,'b-')
xlabel('Time')
ylabel('SSE')
title('Error Dynamics - Data Assimilation')
grid on
