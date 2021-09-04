clc; clear; close all;
c = 299792458;
n_2 = 6e-18;
beta_T = 4.5e-12;
xi_e = 1.25;
tau = 10e-9;
h = 6.63e-34;
hbar = h / 2 / pi;
global alpha_l gamma_e Delta_beta Calpha_fp Calpha_fs Calpha_fi
lambda_p = 1.56e-6;
k_p = 2 * pi / lambda_p;
omega_p = 2 * pi * c / lambda_p;

load('Fig_2_b_proposed.mat')
lambda = 1.45e-6:0.005e-6:1.65e-6;
[p, ~, mu] = polyfit(lambda', loss, 3);
alpha_l = 10^(polyval(p, 1.56e-6, [], mu) / 10);
[p2, ~, mu2] = polyfit(lambda', beta, 3);
beta_p = polyval(p2, lambda_p, [], mu2);

load('Fig_4_a_proposed_50nm.mat')
dS = .75e-8^2;
for i = 1:21
    nc = abs(n(212, 220, i));
    absF2 = abs(Ex(:, :, i)).^2 + abs(Ey(:, :, i)).^2 + abs(Ez(:, :, i)).^2;
    FcdotF = Ex(:, :, i).^2 + Ey(:, :, i).^2 + Ez(:, :, i).^2;
    numerator = (sum(sum(abs(n(:, :, i)).^2 .* absF2)))^2;
%     numerator = (sum(sum((abs(n(:, :, i)).^2 .* absF2)) * (abs(n(:, :, i)) == nc)))^2;
%     denominator = sum(sum(2 * absF2 .* (conj(Ex(:, :, i)) .* Ex(:, :, i) + conj(Ey(:, :, i)) .* Ey(:, :, i) + conj(Ez(:, :, i)) .* Ez(:, :, i)) + FcdotF .* (conj(Ex(:, :, i)).^2 + conj(Ey(:, :, i)).^2 + conj(Ez(:, :, i)).^2)));
    denominator = sum(sum((2 * absF2 .* (conj(Ex(:, :, i)) .* Ex(:, :, i) + conj(Ey(:, :, i)) .* Ey(:, :, i) + conj(Ez(:, :, i)) .* Ez(:, :, i)) + FcdotF .* (conj(Ex(:, :, i)).^2 + conj(Ey(:, :, i)).^2 + conj(Ez(:, :, i)).^2)) .* (abs(n(:, :, i)) == nc)));
    Aeff(i) = real(3 / ng(i)^2 / nc^2 * numerator * dS / denominator);
end
[p, ~, mu] = polyfit(lambda', real(Aeff), 3);
Aeff = polyval(p, lambda_p, [], mu);

eta = [];
for lambda_s = 1.48e-6:.001e-6:1.62e-6
    disp(lambda_s)
    k_s = 2 * pi / lambda_s;
    beta_s = polyval(p2, lambda_s, [], mu2);
    k_i = 2 * k_p - k_s;
    lambda_i = 2 * pi / k_i;
    beta_i = polyval(p2, lambda_i, [], mu2);
    Delta_beta = beta_s + beta_i - 2 * beta_p;
    gamma_e = xi_e * omega_p * n_2 / c / Aeff + 1i * xi_e * beta_T / 2 / Aeff;
    Calpha_fp = 6.04e-10 * lambda_p^2 * xi_e * beta_T * tau / 2 / hbar / omega_p / Aeff^2;
    Calpha_fs = 6.04e-10 * lambda_s^2 * xi_e * beta_T * tau / 2 / hbar / omega_p / Aeff^2;
    Calpha_fi = 6.04e-10 * lambda_i^2 * xi_e * beta_T * tau / 2 / hbar / omega_p / Aeff^2;
    [z, A] = ode45(@CME, [0:0.003e-2:3e-2], [sqrt(10^(15 / 10) * 1e-2), sqrt(10^(-10 / 10) * 1e-2), 0]);
    eta(:, end + 1) = (abs(A(:, 3)).^2 ./ abs(A(:, 2)).^2);
end
lambda_s = 1.48e-6:.002e-6:1.62e-6;
imagesc(lambda_s * 1e9, z * 1e2, log10(eta))
% colormap(jet)
colorbar
axis xy
xlabel('Signal Wavelength [nm]', 'fontsize', 16)
ylabel('Waveguide Length [cm]', 'fontsize', 16)