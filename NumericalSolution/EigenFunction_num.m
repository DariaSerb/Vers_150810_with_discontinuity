function [X,uest] = EigenFunction_num(dTau)
global P;
P = Parameters;
PointCnt = P.PointCnt;

u = InitialApproximation;
% uest = zeros(P.ModeEst,P.PointCnt+2);

[alp_num,DerivLambda,Du] = AlpShapes_num;
Du = Du(1:PointCnt,:);

% u = u + alp*u*dTau;
uest = u(1:P.ModeEst,:) + dTau*Du';

X = P.X;
[q,qs] = q_calc(X);
X = X + q*dTau;
end 

function u = InitialApproximation
global P;
[lambda, u, V, stiffdq, massdq, mui] = eig_val_xc(P); 
end