clc; clear; close all;
c = 3e8;
lambda = 1.45e-6:0.005e-6:1.65e-6;
load('Fig_3_b_proposed_100nm.mat')
omega = 2 * pi * f;
[p, ~, mu] = polyfit(omega, vg.^(-1), 4);
p(4) = - 4 * p(1) * mu(1)^3 / mu(2)^4 + 3 * p(2) * mu(1)^2 / mu(2)^3 - 2 * p(3) * mu(1) / mu(2)^2 + p(4) / mu(2);
p(3) = 6 * p(1) * mu(1)^2 / mu(2)^4 - 3 * p(2) * mu(1) / mu(2)^3 + p(3) / mu(2)^2;
p(2) = -4 * p(1) * mu(1) / mu(2)^4 + p(2) / mu(2)^3;
p(1) = p(1) / mu(2)^4;
dN_domega = 4 * p(1) * omega.^3 + 3 * p(2) * omega.^2 + 2 * p(3) * omega + p(4);
plot(lambda * 1e9, dN_domega * 1e24, 'm-', 'linewidth', 2)
hold on
grid on
xlabel('Wavelenth [nm]', 'fontsize', 16)
ylabel('$\beta_2=\frac{\partial^2\beta}{\partial\omega^2}\,\mathrm{\,at\,1560\,nm\,[ps^2 / m]}$', 'interpreter', 'latex', 'fontsize', 16)

load('Fig_3_b_proposed_150nm.mat')
[p, ~, mu] = polyfit(omega, vg.^(-1), 4);
p(4) = - 4 * p(1) * mu(1)^3 / mu(2)^4 + 3 * p(2) * mu(1)^2 / mu(2)^3 - 2 * p(3) * mu(1) / mu(2)^2 + p(4) / mu(2);
p(3) = 6 * p(1) * mu(1)^2 / mu(2)^4 - 3 * p(2) * mu(1) / mu(2)^3 + p(3) / mu(2)^2;
p(2) = -4 * p(1) * mu(1) / mu(2)^4 + p(2) / mu(2)^3;
p(1) = p(1) / mu(2)^4;
dN_domega = 4 * p(1) * omega.^3 + 3 * p(2) * omega.^2 + 2 * p(3) * omega + p(4);
plot(lambda * 1e9, dN_domega * 1e24, 'g-', 'linewidth', 2)

load('Fig_3_b_proposed_200nm.mat')
[p, ~, mu] = polyfit(omega, vg.^(-1), 4);
p(4) = - 4 * p(1) * mu(1)^3 / mu(2)^4 + 3 * p(2) * mu(1)^2 / mu(2)^3 - 2 * p(3) * mu(1) / mu(2)^2 + p(4) / mu(2);
p(3) = 6 * p(1) * mu(1)^2 / mu(2)^4 - 3 * p(2) * mu(1) / mu(2)^3 + p(3) / mu(2)^2;
p(2) = -4 * p(1) * mu(1) / mu(2)^4 + p(2) / mu(2)^3;
p(1) = p(1) / mu(2)^4;
dN_domega = 4 * p(1) * omega.^3 + 3 * p(2) * omega.^2 + 2 * p(3) * omega + p(4);
plot(lambda * 1e9, dN_domega * 1e24, 'r-', 'linewidth', 2)

load('Fig_3_b_proposed_250nm.mat')
[p, ~, mu] = polyfit(omega, vg.^(-1), 4);
p(4) = - 4 * p(1) * mu(1)^3 / mu(2)^4 + 3 * p(2) * mu(1)^2 / mu(2)^3 - 2 * p(3) * mu(1) / mu(2)^2 + p(4) / mu(2);
p(3) = 6 * p(1) * mu(1)^2 / mu(2)^4 - 3 * p(2) * mu(1) / mu(2)^3 + p(3) / mu(2)^2;
p(2) = -4 * p(1) * mu(1) / mu(2)^4 + p(2) / mu(2)^3;
p(1) = p(1) / mu(2)^4;
dN_domega = 4 * p(1) * omega.^3 + 3 * p(2) * omega.^2 + 2 * p(3) * omega + p(4);
plot(lambda * 1e9, dN_domega * 1e24, 'b-', 'linewidth', 2)

yline(0, 'linewidth', 2);
ylim([-0.6, 0.3])
legend('100 nm', '150 nm', '200 nm', '250 nm', 'location', 'northwest', 'fontsize', 12)