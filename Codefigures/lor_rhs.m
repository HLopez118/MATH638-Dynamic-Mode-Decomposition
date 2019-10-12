function rhs = lor_rhs(t,x,dummy,s,b,r)
    
rhs=[s*(-x(1)+x(2))
     r*(x(1))-x(2)-(x(1)*x(3))
     x(1)*x(2)-(b*x(3))];