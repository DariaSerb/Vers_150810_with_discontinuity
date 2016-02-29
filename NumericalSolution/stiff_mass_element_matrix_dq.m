function [P0dq, K0dq] = stiff_mass_element_matrix_dq(xn,xle,A1,A2)

% definition size of mass matrix and stiffness matrix
xmatrixmassdq = zeros(2,2);
xmatrixkdq = zeros(2,2);
npg = 5;

q = zeros(npg,1); 
qs = zeros(npg,1); 

[x,w] = cal_point_gauss(npg);
d = xn(2) - xn(1);
% numerical integration (0 - xle)
 for ig = 1:npg
    s = (d/2)*(x(ig)+1);
    s1 = s + xn(1);
    
    N2 = s/xle;
    N1 = 1-N2;
    [q(ig),qs(ig)] = q_calc(s1);
    
    N = [N1, N2];
    xmatrixmassdq = xmatrixmassdq + A2*(N'*N)*w(ig)*qs(ig);
    
    dN(1) = -1/xle;
    dN(2) = 1/xle;
    [q(ig),qs(ig)] = q_calc(s1);
    
    dN = [dN(1), dN(2)];
    xmatrixkdq = xmatrixkdq + A1*(dN'*dN)*w(ig)*qs(ig);
 end
    xmatrixmassdq = xmatrixmassdq*(d/2);
    xmatrixkdq = xmatrixkdq*(d/2);
    
    P0dq = xmatrixmassdq;
    K0dq = xmatrixkdq;
end
