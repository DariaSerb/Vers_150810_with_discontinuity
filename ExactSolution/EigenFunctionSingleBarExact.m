 function u = EigenFunctionSingleBarExact(X)
% Eigen function for 1D homogeneity bar
% input:
% X - array of discretization points
% output:
% u - array of eigen function values (several modes)
P = Parameters;
N_mode = P.ModeCnt;
N_x = length(X);

u = zeros(N_mode, N_x);

for n_mode = 1:N_mode
    for n_x = 1:N_x
        u(n_mode, n_x) = u_calc(X(n_x), n_mode, P.L);
    end
end
end

function u = u_calc(x, p, L)
%u_exact exact solution of u in scalar form
   W = (2*p-1)*pi/(2*L);
   u = sin(W*x);
end

