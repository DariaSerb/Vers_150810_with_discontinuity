function alp = AlpShapes
global L;
global S1;
global S2;
global E;
global ro;
global xc;

P = Parameters;

L = P.L;
E = P.E;
ro = P.ro;
S1 = P.S1;
S2 = P.S2;
xc = P.xc;
ModeCnt = P.ModeCnt;

% lambda = eig_val(P);
[lambda, uf] = eig_val(P);
alp = zeros(ModeCnt);

for n = 1:ModeCnt
    for m = 1:ModeCnt
        if n == m
            continue;
        end
       alp(n,m) = integralCalc(lambda(n), lambda(m));
    end
end

end

function alp = integralCalc(lam_i, lam_j)
global L;
global E;
global ro;
global S1;
global S2;
global xc;

N = 200;
dx = L/ N;
X = 0 + dx/2 : dx : (L-dx/2);

[q,qs] = q_calc(X);
[u_i,us_i] = u_calc(X,ro,xc,E,lam_i);
[u_j,us_j] = u_calc(X,ro,xc,E,lam_j);

numer = zeros(N,1);
denom = zeros(N,1);

for n = 1:N
    if X(n) <= xc 
     S = S1;
     else S = S2;
    end
    [numer(n),denom(n)] = integrand_func(S, u_i(n), u_j(n), us_i(n), us_j(n), qs(n), lam_i);
    
end

S_num = ones(1,N)*numer*dx;
S_den = ones(1,N)*denom*dx;

alp = - S_num / S_den;
if lam_i - lam_j ~= 0
   alp = alp / (lam_i - lam_j);   
end
end


function [numer,denom] = integrand_func(S, u_i, u_j, us_i, us_j, qs, lam_i)
%calculation integrand function
global E;
global ro;
denom = ro*S*u_j^2;
numer =(E * S * us_i * us_j + lam_i * ro * S * u_i * u_j)*qs;
end


function [u,us] = u_calc(X, ro, xc, E, lambda)
[N,M] = size(X);
u = zeros(N,M);
us = zeros(N,M);

for n = 1:N
    for m = 1:M
     [u(n,m), us(n,m)] = u_calc_scalar(X(n,m), ro, xc, E, lambda);
    end
end
end


function [u,us] = u_calc_scalar(x, ro, xc, E, lambda)
global S1;
global S2;
%u_calc analytical calculation of derivative of u
B1 = 1;
W = sqrt((ro/E)*lambda);
A2 = B1*sin(W*xc);
B2 = (S1/S2)*B1*cos(W*xc);
 if x <= xc
  u = B1*sin(W*x);    
  else u = A2*cos(W*(x-xc)) + B2*sin(W*(x-xc));
 end

 if x <= xc
  us = B1*W*cos(W*x);    
  else us = -A2*W*sin(W*(x-xc)) + B2*cos(W*(x-xc));
 end
end