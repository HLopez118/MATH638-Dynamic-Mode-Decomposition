function dfdt = VdP_eq(t,y,mu)
dx = y(2);
dy = -y(1) + mu*(1-(y(1)^2))*y(2);

dfdt = [dx; dy];
end