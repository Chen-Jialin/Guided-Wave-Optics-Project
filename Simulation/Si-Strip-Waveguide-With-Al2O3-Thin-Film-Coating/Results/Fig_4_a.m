%%
clc; clear; close all;
figure(1)
dS = .75e-8^2;
load('Fig_4_a_SiO2_cladded.mat')
for i = 1:21
    nc = abs(n(212, 220, i));
    absF2 = abs(Ex(:, :, i)).^2 + abs(Ey(:, :, i)).^2 + abs(Ez(:, :, i)).^2;
    FcdotF = Ex(:, :, i).^2 + Ey(:, :, i).^2 + Ez(:, :, i).^2;
    numerator = (sum(sum(abs(n(:, :, i)).^2 .* absF2)))^2;
%     numerator = (sum(sum((abs(n(:, :, i)).^2 .* absF2)) * (abs(n(:, :, i)) == nc)))^2;
%     denominator = sum(sum(2 * absF2 .* (conj(Ex(:, :, i)) .* Ex(:, :, i) + conj(Ey(:, :, i)) .* Ey(:, :, i) + conj(Ez(:, :, i)) .* Ez(:, :, i)) + FcdotF .* (conj(Ex(:, :, i)).^2 + conj(Ey(:, :, i)).^2 + conj(Ez(:, :, i)).^2)));
    denominator = sum(sum((2 * absF2 .* (conj(Ex(:, :, i)) .* Ex(:, :, i) + conj(Ey(:, :, i)) .* Ey(:, :, i) + conj(Ez(:, :, i)) .* Ez(:, :, i)) + FcdotF .* (conj(Ex(:, :, i)).^2 + conj(Ey(:, :, i)).^2 + conj(Ez(:, :, i)).^2)) .* (abs(n(:, :, i)) == nc)));
    Aeff(i) = 3 / ng(i)^2 / nc^2 * numerator * dS / denominator;
end
plot(lambda * 1e9, real(Aeff) * 1e12, 'r--', 'linewidth', 2)
% plot(lambda, mode_effective_area, 'r--', 'linewidth', 2)
hold on
grid on
xlabel('Wavelength [nm]', 'fontsize', 16)
ylabel('Effective Mode Area [\mu{}m^2]', 'fontsize', 16)

load('Fig_4_a_proposed_50nm.mat')
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
plot(lambda * 1e9, Aeff * 1e12, 'm-', 'linewidth', 2)
% plot(lambda, mode_effective_area, 'm-', 'linewidth', 2)

load('Fig_4_a_proposed_100nm.mat')
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
plot(lambda * 1e9, Aeff * 1e12, 'k-', 'linewidth', 2)
% plot(lambda, mode_effective_area, 'k-', 'linewidth', 2)

load('Fig_4_a_proposed_150nm.mat')
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
plot(lambda * 1e9, Aeff * 1e12, 'r-', 'linewidth', 2)
% plot(lambda, mode_effective_area, 'r-', 'linewidth', 2)

load('Fig_4_a_proposed_200nm.mat')
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
plot(lambda * 1e9, Aeff * 1e12, 'g-', 'linewidth', 2)
% plot(lambda, mode_effective_area, 'g-', 'linewidth', 2)

load('Fig_4_a_proposed_250nm.mat')
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
plot(lambda * 1e9, Aeff * 1e12, 'b-', 'linewidth', 2)
% plot(lambda, mode_effective_area, 'b-', 'linewidth', 2)

ylim([0.07, 0.11])
legend('Silia-cladded', 'Proposed, T = 50 nm', 'Proposed, T = 100 nm', 'Proposed, T = 150 nm', 'Proposed, T = 200 nm', 'Proposed, T = 250 nm', 'fontsize', 12, 'location', 'southeast')
%%
figure(2)
load('Fig_4_a_proposed_insert.mat')
imagesc(x * 1e6, y * 1e6, abs(mode1_Ex)')
hold on
plot([-1.5, 1.5], [0, 0], 'k-', 'linewidth', 2)
plot([-.23, -.23], [0, .25], 'k-', 'linewidth', 2)
plot([.23, .23], [0, .25], 'k-', 'linewidth', 2)
plot([-.23, .23], [.25, .25], 'k-', 'linewidth', 2)
plot([-1.5, -.28], [0.05, 0.05], 'k-', 'linewidth', 2)
plot([.28, 1.5], [0.05, 0.05], 'k-', 'linewidth', 2)
plot([-.28, -.28], [0.05, 0.3], 'k-', 'linewidth', 2)
plot([.28, .28], [0.05, 0.3], 'k-', 'linewidth', 2)
plot([-.28, .28], [0.3, 0.3], 'k-', 'linewidth', 2)
colormap(jet)
colorbar
axis xy
axis equal
xlim([-1, 1])
ylim([-0.5, 1])
xlabel('x / \mu{}m', 'fontsize', 16)
ylabel('y / \mu{}m', 'fontsize', 16)