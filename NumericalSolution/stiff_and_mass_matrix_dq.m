function [Pdq, Kdq] = stiff_and_mass_matrix_dq(xle,phi1,phi2,A11d,A21d,A12d,A22d)

% definition size of mass matrix and stiffness matrix
xmatrixmassdq = zeros(4,4);
xmatrixkdq = zeros(4,4);
npg=5;
[x,w] = cal_point_gauss(npg);
d = abs(phi1);
% numerical integration (0 - d)
for ig = 1:npg
    s = (d/2)*(x(ig)+1);
    N2 = s/xle;
    N1 = 1-N2;
    [q,qs] = q_calc(s);
    
    F2 = calculation_f2(xle,phi1,phi2,s);
    N = [N1, N2, N1*F2, N2*F2];
    xmatrixmassdq = xmatrixmassdq+A21d*(N'*N)*w(ig)*qs(ig);
    
    dF2 = calculation_df2(xle,phi1,phi2,s);
    
    dN(1) = -1/xle;
    dN(2) = 1/xle;
    dN(3) = dN(1)*F2+N(1)*dF2;
    dN(4) = dN(2)*F2+N(2)*dF2;
    dN = [dN(1), dN(2), dN(3), dN(4)];
    xmatrixkdq = xmatrixkdq+A11d*(dN'*dN)*w(ig)*qs(ig);
    
end
    xmatrixmassdq = xmatrixmassdq*(d/2);
    xmatrixkdq = xmatrixkdq*(d/2);

    K_minusdq = xmatrixkdq;
    P_minusdq = xmatrixmassdq;
    
    xmatrixmassdq = zeros(4,4);
    xmatrixkdq = zeros(4,4);
    
 % numerical integration (d - xle)
for ig = 1:npg
    s =((xle-d)/2)*(x(ig)+1)+d;
    N2 = s/xle;
    N1 = 1-N2;
    [q,qs] = q_calc(s);
    
    F2 = calculation_f2(xle,phi1,phi2,s);
    N = [N1, N2, N1*F2, N2*F2];
    xmatrixmassdq = xmatrixmassdq+A22d*(N'*N)*w(ig)*qs(ig);
    
    dF2 = calculation_df2(xle,phi1,phi2,s);
    
    dN(1) = -1/xle;
    dN(2) = 1/xle;
    dN(3) = dN(1)*F2+N(1)*dF2;
    dN(4) = dN(2)*F2+N(2)*dF2;
    dN = [dN(1), dN(2), dN(3), dN(4)];
    xmatrixkdq = xmatrixkdq+A12d*(dN'*dN)*w(ig)*qs(ig);
    
end

    xmatrixmassdq = xmatrixmassdq*((xle-d)/2);
    xmatrixkdq = xmatrixkdq*((xle-d)/2);   
    
    K_plusdq = xmatrixkdq;
    P_plusdq = xmatrixmassdq;
    
    Kdq = K_minusdq + K_plusdq;
    Pdq = P_minusdq + P_plusdq;
    
end