function f_d= discrete_runge_kutta4(f, T)
k1_0 = @(t,x) (  f(t,x) );
k2_0 = @(t,x) ( f(t,x + k1_0(t,x)*T/2) );
k3_0 = @(t,x) ( f(t,x + k2_0(t,x)*T/2) );
k4_0 = @(t,x) ( f(t,x + k1_0(t,x)*T) );
f_d = @(t,x) ( x + (T/6) * ( k1_0(t,x) + 2*k2_0(t,x) + 2*k3_0(t,x) + k4_0(t,x)  )   );
end


