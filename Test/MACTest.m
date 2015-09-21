function MACTest
rootpath = startupTest;

P = Parameters;

N_tau = 20;
dTau = 0.5;
% dTaumax = funcdtm(fs_analyt);
% disp(dtmax);
% dtau = linspace(0, dTaumax, N_tau);
dtau = linspace(0, dTau, N_tau);

mac = zeros(P.ModeCnt, N_tau);

for n = 1:N_tau
    [X,u_num] = EigenFunction(dtau(n));
    u_ext = EigenFunctionSectionBarExact(X);
    mac(:,n) = diag(MAC(u_num, u_ext));   
end

figure(1);
plot(dtau, mac);
grid on;
axis([-inf inf -inf 1]);

stopdown(rootpath);
end