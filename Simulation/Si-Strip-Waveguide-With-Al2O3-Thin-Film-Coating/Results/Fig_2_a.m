clc; clear; close all;
c = 299792458;
lambda_0 = 1.56e-6;
omega_0 = 2 * pi * c / lambda_0;
load('Fig_2_a_Air_cladded.mat')
% for i = 1:51
%     disp(i)
%     % 2nd order polyfit
% %     [p(:, i), ~, mu] = polyfit(2 * pi * f, beta(:, i), 2);
% %     p(:, i) = real(p(:, i));
% %     p(1, i) = p(1, i) / mu(2)^2;
%     % 3rd order polyfit
% %     [p(:, i), ~, mu] = polyfit(2 * pi * f, real(beta(:, i)), 3);
% %     p(2, i) = - 3 * p(1, i) * mu(1) / mu(2)^3 + p(2, i) / mu(2)^2;
% %     p(1, i) = p(1, i) / mu(2)^3;
%     % 4th order polyfit
%     [p(:, i), ~, mu] = polyfit(2 * pi * f, real(beta(:, i)), 4);
%     p(3, i) = 6 * p(1, i) * mu(1)^2 / mu(2)^4 - 3 * p(2, i) * mu(1) / mu(2)^3 + p(3, i) / mu(2)^2;
%     p(2, i) = -4 * p(1, i) * mu(1) / mu(2)^4 + p(2, i) / mu(2)^3;
%     p(1, i) = p(1, i) / mu(2)^4;
% end
% % plot(W * 1e9, 2 * p(1, :) * 1e24, 'mo');
% % plot(W * 1e9, (6 * p(1, :) * omega_0 + 2 * p(2, :)) * 1e24, 'mo');
% plot(W * 1e9, (12 * p(1, :) * omega_0^2 + 6 * p(2, :) * omega_0 + 2 * p(3, :)) * 1e24, 'm-o')
% hold on
% grid on
% xlabel('Width [nm]', 'fontsize', 16)
% ylabel('$\beta_2=\frac{\partial^2\beta}{\partial\omega^2}\,\mathrm{[ps^2 / m]}$', 'interpreter', 'latex', 'fontsize', 16)
% 
% for i = 1:51
%     disp(i)
%     % 2nd order polyfit
% %     [p(:, i), ~, mu] = polyfit(2 * pi * f, vg(:, i).^(-1), 2);
% %     p(2, i) = - 2 * p(1, i) * mu(1) / mu(2)^2 + p(2, i) / mu(2);
% %     p(1, i) = p(1, i) / mu(2)^2;
%     % 3rd order polyfit
% %     [p(:, i), ~, mu] = polyfit(2 * pi * f, real(neff(:, i)), 3);
% %     p(3, i) = 3 * p(1, i) * mu(1)^2 / mu(2)^3 - 2 * p(2, i) * mu(1) / mu(2)^2 + p(3, i) / mu(2);
% %     p(2, i) = - 3 * p(1, i) * mu(1) / mu(2)^3 + p(2, i) / mu(2)^2;
% %     p(1, i) = p(1, i) / mu(2)^3;
%     % 4th order polyfit
%     [p(:, i), ~, mu] = polyfit(2 * pi * f, real(neff(:, i)), 4);
%     p(4, i) = - 4 * p(1, i) * mu(1)^3 / mu(2)^4 + 3 * p(2, i) * mu(1)^2 / mu(2)^3 - 2 * p(3, i) * mu(1) / mu(2)^2 + p(4, i) / mu(2);
%     p(3, i) = 6 * p(1, i) * mu(1)^2 / mu(2)^4 - 3 * p(2, i) * mu(1) / mu(2)^3 + p(3, i) / mu(2)^2;
%     p(2, i) = -4 * p(1, i) * mu(1) / mu(2)^4 + p(2, i) / mu(2)^3;
%     p(1, i) = p(1, i) / mu(2)^4;
% end
% % 1st order derivative
% % dN_domega = 2 * p(1, :) * omega_0 + p(2, :);
% % dN_domega = 3 * p(1, :) * omega_0^2 + 2 * p(2, :) * omega_0 + p(3, :);
% dN_domega = 4 * p(1, :) * omega_0^3 + 3 * p(2, :) * omega_0^2 + 2 * p(3, :) * omega_0 + p(4, :);
% % 2nd order derivative
% % d2N_domega2 = 2 * p(1, :);
% % d2N_domega2 = 6 * p(1, :) * omega_0 + 2 * p(2, :);
% d2N_domega2 = 12 * p(1, :) * omega_0^2 + 6 * p(2, :) * omega_0 + 2 * p(3, :);
% plot(W * 1e9, 2 * pi / c * (2 * dN_domega + omega_0 * d2N_domega2) * 1e24, 'm-o');
% hold on
% grid on
% xlabel('Width [nm]', 'fontsize', 16)
% ylabel('$\beta_2=\frac{\partial^2\beta}{\partial\omega^2}\,\mathrm{[ps^2 / m]}$', 'interpreter', 'latex', 'fontsize', 16)

for i = 1:51
    disp(i)
    % 2nd order polyfit
%     [p(:, i), ~, mu] = polyfit(2 * pi * f, vg(:, i).^(-1), 2);
%     p(2, i) = - 2 * p(1, i) * mu(1) / mu(2)^2 + p(2, i) / mu(2);
%     p(1, i) = p(1, i) / mu(2)^2;
    % 3rd order polyfit
%     [p(:, i), ~, mu] = polyfit(2 * pi * f, vg(:, i).^(-1), 3);
%     p(3, i) = 3 * p(1, i) * mu(1)^2 / mu(2)^3 - 2 * p(2, i) * mu(1) / mu(2)^2 + p(3, i) / mu(2);
%     p(2, i) = - 3 * p(1, i) * mu(1) / mu(2)^3 + p(2, i) / mu(2)^2;
%     p(1, i) = p(1, i) / mu(2)^3;
    % 4th order polyfit
    [p(:, i), ~, mu] = polyfit(2 * pi * f, vg(:, i).^(-1), 4);
    p(4, i) = - 4 * p(1, i) * mu(1)^3 / mu(2)^4 + 3 * p(2, i) * mu(1)^2 / mu(2)^3 - 2 * p(3, i) * mu(1) / mu(2)^2 + p(4, i) / mu(2);
    p(3, i) = 6 * p(1, i) * mu(1)^2 / mu(2)^4 - 3 * p(2, i) * mu(1) / mu(2)^3 + p(3, i) / mu(2)^2;
    p(2, i) = -4 * p(1, i) * mu(1) / mu(2)^4 + p(2, i) / mu(2)^3;
    p(1, i) = p(1, i) / mu(2)^4;
end
% dN_domega = 2 * p(1, :) * omega_0 + p(2, :);
% dN_domega = 3 * p(1, :) * omega_0^2 + 2 * p(2, :) * omega_0 + p(3, :);
dN_domega = 4 * p(1, :) * omega_0^3 + 3 * p(2, :) * omega_0^2 + 2 * p(3, :) * omega_0 + p(4, :);
plot(W * 1e9, dN_domega * 1e24, 'm-o', 'linewidth', 1)
hold on
grid on
xlabel('Width [nm]', 'fontsize', 16)
ylabel('$\beta_2=\frac{\partial^2\beta}{\partial\omega^2}\,\mathrm{\,at\,1560\,nm\,[ps^2 / m]}$', 'interpreter', 'latex', 'fontsize', 16)

load('Fig_2_a_Al2O3_coated.mat')
for i = 1:51
    disp(i)
    % 4th order polyfit
    [p(:, i), ~, mu] = polyfit(2 * pi * f, vg(:, i).^(-1), 4);
    p(4, i) = - 4 * p(1, i) * mu(1)^3 / mu(2)^4 + 3 * p(2, i) * mu(1)^2 / mu(2)^3 - 2 * p(3, i) * mu(1) / mu(2)^2 + p(4, i) / mu(2);
    p(3, i) = 6 * p(1, i) * mu(1)^2 / mu(2)^4 - 3 * p(2, i) * mu(1) / mu(2)^3 + p(3, i) / mu(2)^2;
    p(2, i) = -4 * p(1, i) * mu(1) / mu(2)^4 + p(2, i) / mu(2)^3;
    p(1, i) = p(1, i) / mu(2)^4;
end
dN_domega = 4 * p(1, :) * omega_0^3 + 3 * p(2, :) * omega_0^2 + 2 * p(3, :) * omega_0 + p(4, :);
plot(W * 1e9, dN_domega * 1e24, 'g-v', 'linewidth', 1)

load('Fig_2_a_SiO2_cladded.mat')
for i = 1:51
    disp(i)
    % 4th order polyfit
    [p(:, i), ~, mu] = polyfit(2 * pi * f, vg(:, i).^(-1), 4);
    p(4, i) = - 4 * p(1, i) * mu(1)^3 / mu(2)^4 + 3 * p(2, i) * mu(1)^2 / mu(2)^3 - 2 * p(3, i) * mu(1) / mu(2)^2 + p(4, i) / mu(2);
    p(3, i) = 6 * p(1, i) * mu(1)^2 / mu(2)^4 - 3 * p(2, i) * mu(1) / mu(2)^3 + p(3, i) / mu(2)^2;
    p(2, i) = -4 * p(1, i) * mu(1) / mu(2)^4 + p(2, i) / mu(2)^3;
    p(1, i) = p(1, i) / mu(2)^4;
end
dN_domega = 4 * p(1, :) * omega_0^3 + 3 * p(2, :) * omega_0^2 + 2 * p(3, :) * omega_0 + p(4, :);
plot(W * 1e9, dN_domega * 1e24, 'r-d', 'linewidth', 1)

load('Fig_2_a_proposed.mat')
for i = 1:51
    disp(i)
    % 4th order polyfit
    [p(:, i), ~, mu] = polyfit(2 * pi * f, vg(:, i).^(-1), 4);
    p(4, i) = - 4 * p(1, i) * mu(1)^3 / mu(2)^4 + 3 * p(2, i) * mu(1)^2 / mu(2)^3 - 2 * p(3, i) * mu(1) / mu(2)^2 + p(4, i) / mu(2);
    p(3, i) = 6 * p(1, i) * mu(1)^2 / mu(2)^4 - 3 * p(2, i) * mu(1) / mu(2)^3 + p(3, i) / mu(2)^2;
    p(2, i) = -4 * p(1, i) * mu(1) / mu(2)^4 + p(2, i) / mu(2)^3;
    p(1, i) = p(1, i) / mu(2)^4;
end
dN_domega = 4 * p(1, :) * omega_0^3 + 3 * p(2, :) * omega_0^2 + 2 * p(3, :) * omega_0 + p(4, :);
plot(W * 1e9, dN_domega * 1e24, 'b-*', 'linewidth', 1)

yline(0, 'linewidth', 1);
ylim([-10, 10])
legend('Air-cladded', 'Alumina-coated', 'Silica-cladded', 'Proposed', 'location', 'southeast', 'fontsize', 16)
