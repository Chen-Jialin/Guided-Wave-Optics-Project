clc; clear; close all;
c = 299792458;
lambda_0 = 1.56e-6;
omega_0 = 2 * pi * c / lambda_0;
load('Fig_3_a_proposed-100.mat')

for i = 1:36
    disp(i)
    % 4th order polyfit
    [p(:, i), ~, mu] = polyfit(2 * pi * f, vg(:, i).^(-1), 4);
    p(4, i) = - 4 * p(1, i) * mu(1)^3 / mu(2)^4 + 3 * p(2, i) * mu(1)^2 / mu(2)^3 - 2 * p(3, i) * mu(1) / mu(2)^2 + p(4, i) / mu(2);
    p(3, i) = 6 * p(1, i) * mu(1)^2 / mu(2)^4 - 3 * p(2, i) * mu(1) / mu(2)^3 + p(3, i) / mu(2)^2;
    p(2, i) = -4 * p(1, i) * mu(1) / mu(2)^4 + p(2, i) / mu(2)^3;
    p(1, i) = p(1, i) / mu(2)^4;
end
dN_domega = 4 * p(1, :) * omega_0^3 + 3 * p(2, :) * omega_0^2 + 2 * p(3, :) * omega_0 + p(4, :);
plot(W * 1e9, dN_domega * 1e24, 'm-o', 'linewidth', 1)
hold on
grid on
xlabel('Width [nm]', 'fontsize', 16)
ylabel('$\beta_2=\frac{\partial^2\beta}{\partial\omega^2}\,\mathrm{\,at\,1560\,nm\,[ps^2 / m]}$', 'interpreter', 'latex', 'fontsize', 16)

load('Fig_3_a_proposed-150.mat')
for i = 1:36
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

load('Fig_3_a_proposed-200.mat')
for i = 1:36
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

load('Fig_3_a_proposed-250.mat')
for i = 1:36
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
ylim([-1, 1])
legend('100 nm', '150 nm', '200 nm', '250nm', 'location', 'southeast', 'fontsize', 16)