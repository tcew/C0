function [a,b,c] = rsttoabc(r,s,t)
  
% function [a,b,c] = rsttoabc(r,s,t)
% Purpose: Transfer from (r,s,t) -> (a,b,c) coordinates in triangle  

Np = length(r);
a = zeros(Np,1); b = zeros(Np,1); c = zeros(Np,1);

tol = 1e-14;
for n=1:Np  
  if(abs(s(n)+t(n))>tol)
    a(n) = 2*(1+r(n))/(-s(n)-t(n))-1;
  else
    a(n) = -1;
  end
  if(abs(t(n)-1)>tol)
    b(n) = 2*(1+s(n))/(1-t(n))-1;
  else
    b(n) = -1;
  end
end
c = t;
return;
