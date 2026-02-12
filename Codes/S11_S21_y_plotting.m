% -------------------------------
% Données mesurées
% -------------------------------

y = [0, 2, 3, 5, 7, 10];  % Distances (en cm)

% À la fréquence de résonance
S11_dB_fr = [-0.47, -2.21, -3.93, -6.83, -9.96, -22.291];
S21_dB_fr = [-14.79, -3.82, -2.24, -1.08, -0.57, -0.19];

% À la fréquence de pic max
% S11_dB_max = [-18.47, -20, -37.73, -38.86, -35.61, -29.65];
% S21_dB_max = [-1.22, -3.72, -0.29, -0.24, -0.14, -0.15];

% -------------------------------
% Plot séparé pour S11 et S21
% -------------------------------

figure;

% ---- Subplot S11 ----
subplot(2,1,1);
plot(y, S11_dB_fr, 'r-o', 'LineWidth', 2); hold on;
% plot(y, S11_dB_max, 'm-.^', 'LineWidth', 2);
xlabel('Distance y (mm)');
ylabel('S_{11} (dB)');
title('S_{11} vs. Distance');
legend('S_{11} @ f_{res}', 'S_{11} @ f_{max}', 'Location', 'best');
grid on;

% ---- Subplot S21 ----
subplot(2,1,2);
plot(y, S21_dB_fr, 'b--s', 'LineWidth', 2); hold on;
% plot(y, S21_dB_max, 'g:d', 'LineWidth', 2);
xlabel('Distance y (mm)');
ylabel('S_{21} (dB)');
title('S_{21} vs. Distance');
legend('S_{21} @ f_{res}', 'S_{21} @ f_{max}', 'Location', 'best');
grid on;
