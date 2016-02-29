function [X,uest] = EigenFunction(dTau)
global P;
P = Parameters;

alp = AlpShapes; 
u = InitialApproximation;
Du = zeros(P.ModeEst,P.PointCnt);
uest = zeros(P.ModeEst,P.PointCnt);

for i = 1:P.ModeEst
  for k = i+1:P.ModeCnt
    for n = 1:P.PointCnt;
         Du(i,n) = Du(i,n) + alp(i,k)*u(k,n);
    end
  end
end
uest = u(1:P.ModeEst,:) + Du*dTau;

X = P.X;
[q,qs] = q_calc(X);
X = X + q*dTau;
end 

function u = InitialApproximation
global P;
% u = zeros(P.ModeCnt, P.PointCnt);
% [lambda, u] = eig_val(P); 
[lambda, u, V, stiffdq, massdq, mui] = eig_val(P); 
end
