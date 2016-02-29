function MACTest
rootpath = startupTest;

P = Parameters;


% dTaumax = funcdtm(fs_analyt);
% disp(dtmax);
% dtau = linspace(0, dTaumax, N_tau);
N_tau = P.Pointdisc;
dTau = P.dTau;
xinf = P.xinf;

% mac = zeros(P.ModeEst, N_tau);
ModeEst2 = P.ModeEst*P.ModeEst;
macfin = zeros(ModeEst2, N_tau);

for n = 1:N_tau
    [X,u_num] = EigenFunction_num(dTau(n));
    u_ext = EigenFunctionSectionBarExact(X,xinf(n));
%   mac(:,n) = diag(MAC(u_num, u_ext));   
    mac = MAC(u_num, u_ext);   
    macfin(:,n) = reshape(mac', 1, ModeEst2);
end

figure(1);
plot(dTau, macfin(1:3,:), 'LineWidth', 2);
legend('Uappr1, Uex1','Uappr1, Uex2','Uappr1, Uex3',3)
grid on;
axis([-inf inf -inf 1]);
title(['Modal Assurance Criterium (Uappr1, Uex)']);
xlabel('dTau');
ylabel('log(MAC)');

figure(2);
plot(dTau, macfin(4:6,:), 'LineWidth', 2);
legend('Uappr2, Uex1','Uappr2, Uex2','Uappr2, Uex3',3)
grid on;
axis([-inf inf -inf 1]);
title(['Modal Assurance Criterium (Uappr2, Uex)']);
xlabel('dTau');
ylabel('log(MAC)');

figure(3);
plot(dTau, macfin(7:9,:), 'LineWidth', 2);
legend('Uappr3, Uex1','Uappr3, Uex2','Uappr3, Uex3',3)
grid on;
axis([-inf inf -inf 1]);
title(['Modal Assurance Criterium (Uappr3, Uex)']);
xlabel('dTau');
ylabel('log(MAC)');

stopdown(rootpath);
end
