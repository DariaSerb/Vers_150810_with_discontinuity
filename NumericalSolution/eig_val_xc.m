function [lambda, eigvec, V, Kdq, Mdq, mui] = eig_val_xc(P)
% function [lambda, eigvec] = eig_val(P)

L  = P.L;
xc = P.xc;
S1 = P.S1;
S2 = P.S2;
E1 = P.E; E2 = P.E;
ro1 = P.ro; ro2 = P.ro; 
PointCnt = P.PointCnt;
ModeCnt = P.ModeCnt;

% numbEl: number of elements
% xle: length of each element
% addterm: two additional terms
addterm         = 2;
numbEl          = PointCnt;
nodeCoordinates = linspace(0,L,numbEl + 1);
xle             = L/numbEl;

p = 0;
% minimal value and number of node enclosed to discontinuity
[val, id] = min(abs(nodeCoordinates - xc));
% correction of discontinuity position
if (val/xle)<(p/100)
xc = nodeCoordinates(id);
end;

% N: number of nodes
N = length(nodeCoordinates);

% stiffnessmatrix: stiffness matrix
% massmatrix: mass matrix
% totalnodes: the sum of number of nodes and additinal nodes
totalnodes = N+addterm;

stiff   = zeros(totalnodes,totalnodes);
mass    = zeros(totalnodes,totalnodes);
stiffdq = zeros(totalnodes,totalnodes);
massdq  = zeros(totalnodes,totalnodes);
Kdq     = zeros(totalnodes-1, totalnodes-1);
Mdq     = zeros(totalnodes-1, totalnodes-1);

% fct_ls: calculation the value of level set function for each node
fct_ls  = nodeCoordinates - xc;

% Ke1:  elementary stiffness matrix for normal element for 1-st domain
% Me1:  elementary mass matrix for normal element for 1-st domain
% A11d: product E1*S1/xle  for 1-st domain
% A21d: product ro1*S1*xle for 1-st domain
A11d  = E1*S1;
A21d  = ro1*S1;
% Ke1 = (A11d/xle)*[1 -1; -1 1];
% Me1 = (A21d*xle)*[1/3 1/6; 1/6 1/3];

% Ke2:  elementary stiffness matrix for normal element for 2-d domain
% Me2:  elementary mass matrix for normal element for 2-d domain
% A12d: product E2*S2/xle  for 2-d domain
% A22d: product ro2*S2*xle for 2-d domain
A12d  = E2*S2;
A22d  = ro2*S2;
% Ke2 = (A12d/xle)*[1 -1; -1 1];
% Me2 = (A22d*xle)*[1/3 1/6; 1/6 1/3];

% calculation of level set's product of two nodes
for in = 1:numbEl
    
    % xn: used in numerical integration of elementary stiffness and mass matrices
    xn = [nodeCoordinates(in); nodeCoordinates(in+1)];
    
    if fct_ls(in)*fct_ls(in+1)>0 
        if fct_ls(in)<0
            [Me1, Ke1] = stiff_mass_element_matrix(xn, xle, A11d, A21d);
            stiff(in:in+1,in:in+1)   = stiff(in:in+1,in:in+1) + Ke1;
            mass(in:in+1,in:in+1)    = mass(in:in+1,in:in+1) + Me1;
            
            [Me1dq, Ke1dq] = stiff_mass_element_matrix_dq(xn, xle, A11d, A21d);
            stiffdq(in:in+1,in:in+1) = stiffdq(in:in+1,in:in+1) + Ke1dq;
            massdq(in:in+1,in:in+1)  = massdq(in:in+1,in:in+1) + Me1dq;
        end
        if fct_ls(in)>0
            [Me2, Ke2] = stiff_mass_element_matrix(xn, xle, A12d, A22d);
            stiff(in:in+1,in:in+1)   = stiff(in:in+1,in:in+1) + Ke2;
            mass(in:in+1,in:in+1)    = mass(in:in+1,in:in+1) + Me2;
            
            [Me2dq, Ke2dq] = stiff_mass_element_matrix_dq(xn, xle, A12d, A22d);
            stiffdq(in:in+1,in:in+1) = stiffdq(in:in+1,in:in+1) + Ke2dq;
            massdq(in:in+1,in:in+1)  = massdq(in:in+1,in:in+1) + Me2dq;
        end;
    end;
    
    if fct_ls(in)*fct_ls(in+1)==0
        if fct_ls(in+1)==0
            [Me1, Ke1] = stiff_mass_element_matrix(xn, xle, A11d, A21d);
            stiff(in:in+1,in:in+1)   = stiff(in:in+1,in:in+1) + Ke1;
            mass(in:in+1,in:in+1)    = mass(in:in+1,in:in+1) + Me1;    
            
            [Me1dq, Ke1dq] = stiff_mass_element_matrix_dq(xn, xle, A11d, A21d);
            stiffdq(in:in+1,in:in+1) = stiffdq(in:in+1,in:in+1) + Ke1dq;
            massdq(in:in+1,in:in+1)  = massdq(in:in+1,in:in+1) + Me1dq;
        end;
        if fct_ls(in)==0
            [Me2, Ke2] = stiff_mass_element_matrix(xn, xle, A12d, A22d);
            stiff(in:in+1,in:in+1)   = stiff(in:in+1,in:in+1) + Ke2;
            mass(in:in+1,in:in+1)    = mass(in:in+1,in:in+1) + Me2;
            
            [Me2dq, Ke2dq] = stiff_mass_element_matrix_dq(xn, xle, A12d, A22d);
            stiffdq(in:in+1,in:in+1) = stiffdq(in:in+1,in:in+1) + Ke2dq;
            massdq(in:in+1,in:in+1)  = massdq(in:in+1,in:in+1) + Me2dq;
        end;
     continue;
    end;
    
    if  fct_ls(in)*fct_ls(in+1)<0
        phi1 = fct_ls(in);
        phi2 = fct_ls(in+1);
        [xmass, xk] = stiffness_and_mass_matrix(xle,phi1,phi2,A11d,A21d,A12d,A22d);
        
        % Assembling global stiffness matrix for bar
        stiff(in:in+1,in:in+1) = stiff(in:in+1,in:in+1) + xk(1:2,1:2);
        stiff(N+1:N+2,N+1:N+2) = stiff(N+1:N+2,N+1:N+2) + xk(3:4,3:4);
        stiff(in:in+1,N+1:N+2) = stiff(in:in+1,N+1:N+2) + xk(1:2,3:4);
        stiff(N+1:N+2,in:in+1) = stiff(N+1:N+2,in:in+1) + xk(3:4,1:2);
        
        % Assembling global mass matrix for each element of bar
        mass(in:in+1,in:in+1) = mass(in:in+1,in:in+1) + xmass(1:2,1:2);
        mass(N+1:N+2,N+1:N+2) = mass(N+1:N+2,N+1:N+2) + xmass(3:4,3:4);
        mass(in:in+1,N+1:N+2) = mass(in:in+1,N+1:N+2) + xmass(1:2,3:4);
        mass(N+1:N+2,in:in+1) = mass(N+1:N+2,in:in+1) + xmass(3:4,1:2);
        
        [xmassdq, xkdq] = stiff_and_mass_matrix_dq(xn,xle,phi1,phi2,A11d,A21d,A12d,A22d);
        
        % Assembling global stiffness matrix*dq for bar
        stiffdq(in:in+1,in:in+1) = stiffdq(in:in+1,in:in+1) + xkdq(1:2,1:2);
        stiffdq(N+1:N+2,N+1:N+2) = stiffdq(N+1:N+2,N+1:N+2) + xkdq(3:4,3:4);
        stiffdq(in:in+1,N+1:N+2) = stiffdq(in:in+1,N+1:N+2) + xkdq(1:2,3:4);
        stiffdq(N+1:N+2,in:in+1) = stiffdq(N+1:N+2,in:in+1) + xkdq(3:4,1:2);
        
        % Assembling global mass matrix*dq for each element of bar
        massdq(in:in+1,in:in+1) = massdq(in:in+1,in:in+1) + xmassdq(1:2,1:2);
        massdq(N+1:N+2,N+1:N+2) = massdq(N+1:N+2,N+1:N+2) + xmassdq(3:4,3:4);
        massdq(in:in+1,N+1:N+2) = massdq(in:in+1,N+1:N+2) + xmassdq(1:2,3:4);
        massdq(N+1:N+2,in:in+1) = massdq(N+1:N+2,in:in+1) + xmassdq(3:4,1:2);
        
    end
end

 Kdq = stiffdq(2:totalnodes, 2:totalnodes);
 Mdq = massdq(2:totalnodes, 2:totalnodes);

% calculation of matrices D and V 
% D:    diagonal matrix containing the eigenvalues
% V:    matrix whose columns are the corresponding right eigenvectors
% freq: natural frequences
% D = eig(stiff(2:totalnodes, 2:totalnodes), mass(2:totalnodes, 2:totalnodes));
    [V, D] = eig(stiff(2:totalnodes, 2:totalnodes), mass(2:totalnodes, 2:totalnodes));
    lambda = zeros(length(D),1);
    eigvec = zeros(length(V),ModeCnt);
% condnumber: condition number of V (matrix of eigenvectors) 
% condnumber=cond(V);
    for in = 1:length(D)
% lambda(in) = sqrt(D(in,in))/(2*pi);
     lambda(in) = D(in,in);  
    end
     lambda = sort(lambda);
     lambda = lambda(1:ModeCnt);
    
    for in=1:ModeCnt
     eigvec(:,in) = V(:,in);
    end
  eigvec = eigvec';
  eigvec = eigvec(:,1:length(V) - 2);
  mui = diag(V'*mass(2:totalnodes, 2:totalnodes)*V);
 end
 
