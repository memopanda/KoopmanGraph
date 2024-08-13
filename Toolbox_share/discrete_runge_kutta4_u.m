function f_ud= discrete_runge_kutta4_u(f_u, T)
k1 = @(t,x,u) (  f_u(t,x,u) );
k2 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*T/2,u) );
k3 = @(t,x,u) ( f_u(t,x + k2(t,x,u)*T/2,u) );
k4 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*T,u) );
f_ud = @(t,x,u) ( x + (T/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   );
end