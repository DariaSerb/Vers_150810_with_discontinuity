function MACFullMatrixTest
rootpath = startupTest;

dtau = 0.15;

[X,u_num] = EigenFunction(dtau);
u_ext = EigenFunctionSectionBarExact(X);
mac = MAC(u_num, u_ext);   

save 'MATFullMatrix.mat' mac dtau;

stopdown(rootpath);
end