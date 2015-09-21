function [P0dq, K0dq] = stiff_mass_element_matrix_dq(xle,A1,A2)

% definition size of mass matrix and stiffness matrix
xmatrixmass=zeros(2,2);
xmatrixk=zeros(2,2);
npg=5;
[x,w]=cal_point_gauss(npg);
d=xle;
% numerical integration (0 - xle)
 for ig=1:npg
    s=(d/2)*(x(ig)+1);
    N2=s/xle;
    N1=1-N2;
    [q,qs] = q_calc(s);
    
    N=[N1, N2];
    xmatrixmass=xmatrixmass+A2*(N'*N)*w(ig)*qs(ig);
    
    dN(1)=-1/xle;
    dN(2)=1/xle;
    [q,qs] = q_calc(s);
    
    dN=[dN(1), dN(2)];
    xmatrixk=xmatrixk+A1*(dN'*dN)*w(ig)*qs(ig);
 end
    xmatrixmass=xmatrixmass*(d/2);
    xmatrixk=xmatrixk*(d/2);
    
    P0dq = xmatrixmass;
    K0dq = xmatrixk;
end