options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-7, 'Jumps', [0]);

% xh0 = -1;
% T = 1;
% tau = 0.1;
% K0 = 1;
% K1 = 14;
% sol = dde23(@(t,y,z) -K0*y-K1*z, [tau],@(t) xh0,[0,T], options);

% tau = 0.5;
% K0 = 1;
% K1 = 2;

T = 1;
tau = 0.4;
K0 = 1;
K1 = 3;
sol = dde23(@(t,y,z) -K0*y-K1*z, [tau],@(t) xh0,[0,T], options);


figure(5)
clf
hold on
plot([-tau sol.x], [-1 sol.y])
plot([-tau, T], [0, 0], ':k')
xlim([-tau, T])
ylim(abs(xh0)*1.1*[-1, 1])
xlabel('t')
ylabel('x(t)')
title(['$\dot{x}(t) = -', num2str(K0), 'x(t) -', num2str(K1), 'x(t-', num2str(tau),')$'], 'interpreter', 'latex', 'fontsize', 16)
hold off


% function d = curr_delay(t,y)
% 
% d = t - tau;
% end
% 
% function y = history(t)
%     y = -1;
% end
% 
% function dydt = ddex(t,y,Z)
% ylag = Z(:,1);
% % dydt = 0.5*y(1) - ylag(1);
% dydt = -5*ylag(1);
% end