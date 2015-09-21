function [P, K] = stiffness_and_mass_matrix(xle,phi1,phi2,A11d,A21d,A12d,A22d)

% definition size of mass matrix and stiffness matrix
xmatrixmass = zeros(4,4);
xmatrixk = zeros(4,4);
npg=5;
[x,w] = cal_point_gauss(npg);
d = abs(phi1);
% numerical integration (0 - d)
for ig = 1:npg
    s = (d/2)*(x(ig)+1);
    N2 = s/xle;
    N1 = 1-N2;
    
    F2 = calculation_f2(xle,phi1,phi2,s);
    N = [N1, N2, N1*F2, N2*F2];
    xmatrixmass = xmatrixmass+A21d*(N'*N)*w(ig);
    
    dF2 = calculation_df2(xle,phi1,phi2,s);
    
    dN(1) = -1/xle;
    dN(2) = 1/xle;
    dN(3) = dN(1)*F2+N(1)*dF2;
    dN(4) = dN(2)*F2+N(2)*dF2;
    dN = [dN(1), dN(2), dN(3), dN(4)];
    xmatrixk = xmatrixk+A11d*(dN'*dN)*w(ig);
    
end
    xmatrixmass = xmatrixmass*(d/2);
    xmatrixk = xmatrixk*(d/2);

    K_minus = xmatrixk;
    P_minus = xmatrixmass;
    
    xmatrixmass = zeros(4,4);
    xmatrixk = zeros(4,4);
    
 % numerical integration (d - xle)
for ig = 1:npg
    s =((xle-d)/2)*(x(ig)+1)+d;
    N2 = s/xle;
    N1 = 1-N2;
    
    F2 = calculation_f2(xle,phi1,phi2,s);
    N = [N1, N2, N1*F2, N2*F2];
    xmatrixmass = xmatrixmass+A22d*(N'*N)*w(ig);
    
    dF2 = calculation_df2(xle,phi1,phi2,s);
    
    dN(1) = -1/xle;
    dN(2) = 1/xle;
    dN(3) = dN(1)*F2+N(1)*dF2;
    dN(4) = dN(2)*F2+N(2)*dF2;
    dN = [dN(1), dN(2), dN(3), dN(4)];
    xmatrixk = xmatrixk+A12d*(dN'*dN)*w(ig);
    
end

    xmatrixmass = xmatrixmass*((xle-d)/2);
    xmatrixk = xmatrixk*((xle-d)/2);   
    
    K_plus = xmatrixk;
    P_plus = xmatrixmass;
    
    K = K_minus + K_plus;
    P = P_minus + P_plus;
    
end