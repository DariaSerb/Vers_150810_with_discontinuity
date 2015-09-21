function [X,u] = EigenFunction(dTau)
global P;
P = Parameters;

alp = AlpShapes; 
u = InitialApproximation;
u = u + alp*u*dTau;

X = P.X;
[q,qs] = q_calc(X);
X = X + q*dTau;
end 

function u = InitialApproximation
global P;
% u = zeros(P.ModeCnt, P.PointCnt);
[lambda, u] = eig_val(P); 
end