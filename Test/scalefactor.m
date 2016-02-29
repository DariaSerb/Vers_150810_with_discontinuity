function [coef,vn] = scalefactor(v, w)
%SCALEFACTOR is the coefficient for multiplication in analytical eigenshapes
% v is the analytical eigenshape
% w is the numerical eigenshape
P = Parameters;

 for nn = 1:P.ModeEst
  coef(nn) = sum(w(nn,:).* v(nn,:))/sum(v(nn,:).* v(nn,:));
  vn(nn,:) = v(nn,:).*coef(nn);
 end
 
end

