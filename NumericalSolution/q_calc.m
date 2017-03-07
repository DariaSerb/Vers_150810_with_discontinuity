function [q, qs] = q_calc(X)
% 04/03/2017
[N,M] = size(X);
q = zeros(N,M);
qs = zeros(N,M);
for n = 1:N
    for m = 1:M
        [q(n,m), qs(n,m)] = q_calc_scalar(X(n,m));
    end
end
end

function [q, qs] = q_calc_scalar(x)

% x0    = 0.98540;
% delta = 0.3;
% eta   = 0.6;
x0    = 0.9000;
delta = 0.10;
eta   = 0.75;
% parameters of graphic q-function
x2 = x0 - delta/2.;
x1 = x2 - eta;
x3 = x0 + delta/2.;
x4 = x3 + eta;
q1 = 0.0;
q2 = 1.0;

if (x < x1)||(x > x4)
    q = q1;
    qs = 0;
    return;
end
if (x < x3)&&(x > x2)
    q = q2;
    qs = 0; 
    return;
end;

if (x < x2)&&(x > x1)
    [f, fs] = func((x-x1)/(x2-x1));
    q = f;
    qs = (1/(x2-x1))*fs; 
    return;
end;

if (x < x4)&&(x > x3)
    [f, fs] = func(1 - (x-x3)/(x4-x3));
    q = f;
    qs = - (1/(x4-x3))*fs; 
    return;
end;

end

function [y, ys] = func(x)
%f_calc_PH analytical calculation f-function and derivative of f
xx1 = 0; xx2 = 1; qq1 = 0; qq2 = 1;
N = 4;
coef = solv_slaeq(xx1,xx2,qq1,qq2,N);
y = 0; ys = 0; 
 for nn=1:N
  y = y + sum(coef(nn)*hermite(nn-1,x));
  ys = ys + sum(coef(nn)*diffH(nn-1,x));
 end
end
  
 function coef = solv_slaeq(xx1,xx2,qq1,qq2,N)
%SOLV_SLAEQ 
% the physical definition of coefficients of linear algebratic equations system 
A = zeros(N,N);
D = zeros(N,1);
for j=1:N
  A(1,j) = hermite(j-1,xx1);
  A(2,j) = hermite(j-1,xx2); 
  A(3,j) = diffH(j-1,xx1); 
  A(4,j) = diffH(j-1,xx2); 
end
  B = [qq1;qq2;0;0];
  D = linsolve(A,B);
  coef = zeros(N,1);
 for n=1:N
   coef(n) = D(n);
 end
 end

function H = hermite(n,x)
%polinome d'Hermite
 if n==0 
    H = 1;
 end
 if n==1
   H = 2*x;
 end
 if n>1
% physical definition
  H = 2*x*hermite(n-1,x)-2*(n-1)*hermite(n-2,x);    
 end

end

function Hs = diffH(n,x)
%diffH calculate the derivative of H - Hs
 if n==0 
    Hs = 0;
 else
% probabilistic physical definition Hs = n*hermite(n-1,x)
Hs = 2*n*hermite(n-1,x);
 end

end
