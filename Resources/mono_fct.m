function Y = mono_fct(xx,n,p,yy)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% n: num of variables
% p: max. order
if n==2 && p==3

if (~exist('yy','var') || isempty(yy))
Y= [xx(:,:).^2;xx(:,:).^3;  xx(1,:).*xx(2,:); (xx(1,:).^2).*xx(2,:);  xx(1,:).*(xx(2,:).^2);   ones(1, size(xx,2))]; %Y
else
Y = [2*xx(:,:).*yy(:,:); 3*(xx(:,:).^2).*yy(:,:); xx(1,:).*yy(2,:)+xx(2,:).*yy(1,:); xx(1,:).*xx(1,:).*yy(2,:)+2*xx(1,:).*xx(2,:).*yy(1,:);...
                             2*xx(1,:).*xx(2,:).*yy(2,:)+xx(2,:).*xx(2,:).*yy(1,:);   zeros(1,size(yy,2))];  %dY
end

end




end