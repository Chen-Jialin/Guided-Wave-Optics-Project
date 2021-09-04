function dAdz = CME(z, A)
global alpha_l gamma_e Delta_beta Calpha_fp Calpha_fs Calpha_fi
alpha_fp = Calpha_fp * abs(A(1, 1))^4;
alpha_fs = Calpha_fs * abs(A(2, 1))^4;
alpha_fi = Calpha_fi * abs(A(3, 1))^4;
dAdz(1, 1) = - (alpha_l + alpha_fp) / 2 * A(1, 1) + 1i * gamma_e * abs(A(1, 1))^2 * A(1, 1);
dAdz(2, 1) = - (alpha_l + alpha_fs) / 2 * A(2, 1) + 1i * gamma_e * (2 * abs(A(1, 1))^2 * A(2, 1) + A(1, 1)^2 * conj(A(3, 1)) * exp(- 1i * Delta_beta * z));
dAdz(3, 1) = - (alpha_l + alpha_fi) / 2 * A(3, 1) + 1i * gamma_e * (2 * abs(A(1, 1))^2 * A(3, 1) + A(1, 1)^2 * conj(A(2, 1)) * exp(- 1i * Delta_beta * z));