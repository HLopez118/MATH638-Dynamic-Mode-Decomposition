IC = [1; 1];
mu = 3;
tcolm = 10;
trowm = 15;
pmax = tcolm+2*trowm;
pspan = 0:dt:pmax;

[t,Y1] = ode45(@VdP_eq, pspan, IC,[], mu);
Y1 = Y1';