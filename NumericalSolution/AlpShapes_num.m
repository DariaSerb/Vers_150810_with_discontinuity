function [alp_num,DerivLambda,Du] = AlpShapes_num
P = Parameters;

ModeCnt = P.ModeCnt;
ModeEst = P.ModeEst;

[lambda, uf, V, stiffdq, massdq, mui] = eig_val_xc(P);
alp_num = zeros(ModeCnt,ModeEst);

for n = 1:ModeEst
    var = (stiffdq + lambda(n)*massdq)*V(:,n);
    numer = (V(:,n))'*var;
    DerivLambda(1,n) = - numer/mui(n);
    for m = 1:ModeCnt
        if n == m
            continue;
        end
       numer = (V(:,m))'*var; 
       alp_num(m,n) = numer/(mui(m)*(lambda(m) - lambda(n)));
    end
end

Du = V(:,1:ModeCnt)*alp_num;

end
