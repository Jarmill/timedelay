options = ddeset('AbsTol', 1e-10, 'RelTol', 1e-7);

T = 1;
xh0 = -1;
% kappa = 0.4;
% K0 = 1;
% K1 = 4;

kappa = 0.2;
K0 = 1;
K1 = 2;

sol = ddesd(@(t, y, z) -K0*y - K1*z,@(t, y) kappa*t,@(t) xh0,[0,T], options);

figure(1)
clf
hold on
plot(sol.x, sol.y)
plot([0, T], [0, 0], ':k')

xlabel('t')
ylabel('x(t)')
title(['$\dot{x}(t) = -', num2str(K0), 'x(t) -', num2str(K1), 'x(', num2str(kappa),'t)$'], 'interpreter', 'latex', 'fontsize', 16)
hold off


% function d = curr_delay(t,y)
% %DDEX1DELAYS  Delays for using with DDEX1DE.
% 
% d = [ 0.2 * t];
% % d = t - 1;
% end

% function y = history(t)
%     y = -1;
% end

% function dydt = ddex(t,y,Z)
% ylag = Z(:,1);
% dydt = -3*y(1) - 3*ylag(1);
% % dydt = -ylag(1);
% end