function dfdt = Lorenz_eq(t,y,sigma,r,b)
dx = sigma*(y(2)-y(1));
dy = r*y(1) - y(2) - y(1)*y(3);
dz = y(1)*y(2) - b*y(3);

dfdt = [dx; dy; dz];
end