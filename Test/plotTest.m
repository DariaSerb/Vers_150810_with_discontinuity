function plotTest()
% test for plot of numerical and exact solutions
rootpath = startupTest;

P = Parameters;

dTau = P.dTau;
index = 100;
[x,u_num] = EigenFunction_num(dTau(index));

X = P.X;
[q,qs] = q_calc(X);
X = X + q*dTau(index);
xinf = P.xinf;
u_exact = EigenFunctionSectionBarExact(X,xinf(index));

% figure(1);
% plot(x, u_num, X, u_exact);
% legend('Uapprox1','Uapprox2','Uapprox3','Uexact1','Uexact2','Uexact3',3)
% grid on;
% title(['plot of numerical and exact solutions \Delta\tau =',num2str(dTau),'(s)']);
% xlabel('x');
% ylabel('u');

[coef,u_exscale] = scalefactor(u_exact, u_num);

% figure(1);
% plot(x, u_num(1,:), X, u_exact(1,:));
% legend('Uapprox1','Uexact1',3)
% grid on;
% title(['plot of numerical and exact solutions \Delta\tau =',num2str(dTau(index)),'(s)']);
% xlabel('x');
% ylabel('u');

figure(1);
plot(x, u_num, X, u_exact);
legend('Uapprox1','Uapprox2','Uapprox3','Uexact1','Uexact2','Uexact3',3)
grid on;
title(['plot of numerical and exact solutions \Delta\tau =',num2str(dTau(index)),'(s)']);
xlabel('x');
ylabel('u');

figure(2);
plot(x, u_num, X, u_exscale);
legend('Uapprox1','Uapprox2','Uapprox3','Uexact1','Uexact2','Uexact3',3)
grid on;
title(['plot of numerical and exact solutions using a scale factor \Delta\tau =',num2str(dTau(index)),'(s)']);
xlabel('x');
ylabel('u');

% figure(2);
% plot(x, u_num(2,:), X, u_exact(2,:));
% legend('Uapprox2','Uexact2',3)
% grid on;
% title(['plot of numerical and exact solutions \Delta\tau =',num2str(dTau(index)),'(s)']);
% xlabel('x');
% ylabel('u');
% 
% figure(3);
% plot(x, u_num(3,:), X, u_exact(3,:));
% legend('Uapprox3','Uexact3',3)
% grid on;
% title(['plot of numerical and exact solutions \Delta\tau =',num2str(dTau(index)),'(s)']);
% xlabel('x');
% ylabel('u');
% 
% for n = 1:3
%    figure(n+3);
%    plot(x, u_num(n,:), X, u_exscale(n,:));
%    legend('Uapprox','Uexact',3)
%    grid on;
%    title(['plot of numerical and exact solutions using a scale factor \Delta\tau =',num2str(dTau(index)),'(s)']);
% %  title(strcat('Mode #', num2str(n)));
%    xlabel('x [m]');
%    ylabel('eigenshape');
% end

% figure(4);
% plot(x, u_num(1,:), X, u_exscale(1,:));
% legend('Uapprox1','Uexact1',3)
% grid on;
% title(['plot of numerical and exact solutions \Delta\tau =',num2str(dTau),'(s)']);
% xlabel('x');
% ylabel('u');

% figure(5);
% plot(x, u_num(2,:), X, u_exscale(2,:));
% legend('Uapprox2','Uexact2',3)
% grid on;
% title(['plot of numerical and exact solutions \Delta\tau =',num2str(dTau),'(s)']);
% xlabel('x');
% ylabel('u');

% figure(6);
% plot(x, u_num(3,:), X, u_exscale(3,:));
% legend('Uapprox3','Uexact3',3)
% grid on;
% title(['plot of numerical and exact solutions \Delta\tau =',num2str(dTau),'(s)']);
% xlabel('x');
% ylabel('u');

% ifig = 0;
% for im = 1 : P.ModeEst
%     ifig = ifig + 1;
%     figure(ifig)
%     plot(u_exact(im,:),u_num(im,:),'LineWidth',2)
%     grid
% end
stopdown(rootpath);
end
