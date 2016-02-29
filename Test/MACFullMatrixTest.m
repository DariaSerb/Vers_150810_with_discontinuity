ffunction MACFullMatrixTest
rootpath = startupTest;

P = Parameters;

dTau = P.dTau;
xinf = P.xinf;
index = 2;

[X,u_num] =  EigenFunction_num(dTau(index));
u_ext = EigenFunctionSectionBarExact(xinf(index));

mac = MAC(u_num, u_ext);   
% mac = MAC(u_num, u_exscale);   

save 'MATFullMatrix.mat' mac dTau(index);

stopdown(rootpath);
end
