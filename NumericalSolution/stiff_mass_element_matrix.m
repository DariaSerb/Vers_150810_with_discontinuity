function [P0, K0] = stiff_mass_element_matrix(xn, xle, A1, A2)

% definition size of mass matrix and stiffness matrix
xmatrixmass = zeros(2,2);
xmatrixk    = zeros(2,2);
npg   = 5;
[x,w] = cal_point_gauss(npg);
d     = xn(2) - xn(1);
% numerical integration (0 - xle)
 for ig = 1:npg
    s   = (d/2)*(x(ig)+1);
    N2  = s/xle;
    N1  = 1 - N2;
    
    N = [N1, N2];
    xmatrixmass = xmatrixmass + A2*(N'*N)*w(ig);
    
    dN(1) = -1/xle;
    dN(2) =  1/xle;
    dN = [dN(1), dN(2)];
    xmatrixk = xmatrixk + A1*(dN'*dN)*w(ig);
 end
    xmatrixmass = xmatrixmass*(d/2);
    xmatrixk    = xmatrixk*(d/2);
    
    P0 = xmatrixmass;
    K0 = xmatrixk;
end
