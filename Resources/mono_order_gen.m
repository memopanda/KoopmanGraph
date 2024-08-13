function exp_comb = mono_order_gen(n,p,n0)
%  generate order up to p for n variables,e.g. n=2 i.e. {x y}, 
% for p=2, output is [1  0; 2 0; 0 1; 1 1; 0 2]--> [x x^2 y xy y^2]
% n  Number of variables/features
% p   Maximum polynomial degree
% n0 n0=0 include zeros(1,n) otherwise exclude
% YangGuo
vars = cell(n,1);
[vars{1:n}] = ndgrid(0:p);
exp_comb = reshape(cat(n+1, vars{:}), [], n);
exp_comb = exp_comb(sum(exp_comb,2)<=p & sum(exp_comb,2)>0, :);
if exist('n0','var') 
if n0==0
 exp_comb = [zeros(1, n);exp_comb];
end
end
end

