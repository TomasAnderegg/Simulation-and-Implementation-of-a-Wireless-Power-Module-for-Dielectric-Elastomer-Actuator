
clc; 
clear; 
close all;


%% Facteur de Qualité/ Frequence de résonance


% On voudrait avoir une capacité entre des [10pF, 100nF]
% ainsi qu'une inductance dans l'ordre des uH (microH)
% En ce qui concerne la capacité parasite de l'inductance, on la voudrait
% la plus petite possible.

f = 13.56e6;
frequency = logspace(log10(1e6), log10(20e9), 1000);
C_tot = logspace(-12, -7, 1000)+ 5e-12; 
R_par = 1;

valid_w = [];
valid_d_out = [];
valid_d_in = [];
valid_N = [];
valid_s = [];
valid_resistance = [];

L_values = zeros(size(C_tot)); 
Q_values = zeros(size(C_tot)); 

for i = 1:length(C_tot)
    [L_values(i), Q_values(i)] = Param(f, C_tot(i), R_par);
end

% C_coude = Coude(Q_values, C_tot);
% C_coude = C_coude + 20e-12;
% [L_coude, Q] = Param(f, (C_coude), R_par);
% [valid_w, valid_d_out,valid_d_in, valid_N, valid_s,valid_resistance] = Planar_Coil(L_coude);
% Affichage(valid_w, valid_d_out,valid_d_in, valid_N, valid_s,valid_resistance, L_coude);

%Design asymetrical coils
% C_coude = 150e-12;
% [L_coude, Q] = Param(f, C_coude, R_par)
% 
% [valid_w, valid_d_out,valid_d_in, valid_N, valid_s,valid_resistance, valid_longueur_coil] = Planar_Coil(L_coude);
% Affichage(valid_w, valid_d_out,valid_d_in, valid_N, valid_s,valid_resistance, L_coude, valid_longueur_coil);

% % Pour une capacité de 75 pF
C_coude = 75e-12;
[L_digikey, Q] = Param(f, C_coude, R_par);
[valid_w, valid_d_out,valid_d_in, valid_N, valid_s,valid_resistance, valid_longueur_coil] = Planar_Coil(L_digikey);
Affichage(valid_w, valid_d_out,valid_d_in, valid_N, valid_s,valid_resistance, L_digikey, valid_longueur_coil);
L_coude = L_digikey;



% Check the validity of the found components
if (C_coude < 20e-12 || C_coude > 1e-9)
    fprintf('Attention la valeur de C est trop faible\n')
end

if L_coude < 1e-7
    fprintf('Attention la valeur de L est trop faible, quality factor trop faible\n')
end

% Display of the found component 
fprintf('Inductance L = %.3e H\n', L_coude);
fprintf('Quality Factor Q = %.2f\n', Q);
fprintf('Value of C = %.3e F\n', C_coude-5e-12);


% Plotting of the evolution Q in function of C
figure;
yyaxis left
semilogx(C_tot, Q_values, 'b', 'LineWidth', 2);
ylabel('Quality Factor Q');
xlabel('Capacitance C (F)');

grid on;

yyaxis right
semilogx(C_tot, L_values, 'r', 'LineWidth', 2);
ylabel('Inductance L (H)');

title('Evolution of the inductance and Quality Factor in function of C');
legend('Quality Factor Q', 'Inductor L', 'Location', 'best');



% m = metal('copper');
% inductor = spiralInductor(SpiralShape="Circle", InnerDiameter=valid_d_in(2), Width=valid_w(2), Spacing=valid_s(2),NumTurns=valid_N(2),GroundPlaneLength=0.049,GroundPlaneWidth=0.049, Conductor=m);
% ind = inductance(inductor,13.56e6);

% Efficiency(C_coude, L_coude, R_par, frequency,1);

% Functions section

function Q = Quality(L, R_par, frequency)
    Q = (2 * pi * frequency * L) ./ R_par; % pour un circuit en série

    figure;
    semilogx(frequency, Q, 'b', 'LineWidth', 2);
    ylabel('Quality Factor Q');
    xlabel('Frequency (Hz)');
    title('Evolution of the Quality Factor in function of the frequency');
    grid on;
end

function [L, Q] = Param(f, C, R_par)
    
    L = (1 ./ ((2 * pi .* f) .^ 2 .* C));
    Q = (2 * pi .* f .* L) ./ R_par;

end

function C_coude = Coude(Q_values,C) %L-method

    P1 = [log10(C(1)), Q_values(1)];
    P2 = [log10(C(end)), Q_values(end)];

    distances = zeros(size(C));
    for i = 1:length(C)
        P = [log10(C(i)), Q_values(i)];
        distances(i) = abs(det([P2-P1; P-P1])) / norm(P2-P1);
    end
    
    [~, idx_coude] = max(distances);
    C_coude = C(idx_coude);
    
    figure;
    semilogx(C, Q_values, 'b', 'LineWidth', 2);
    hold on;
    plot(C_coude, Q_values(idx_coude), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    fprintf('Elbow capacitor (C) : %.2e F\n', C_coude);
    fprintf('Max Quality Factor(Q) : %.2f\n', Q_values(idx_coude));
    xlabel('Capacitancer C (F)');
    ylabel('Quality Factor Q');
    title('Plotting of the Quality Factor in function of the Capacitance');
    legend('Quality Factor Q(C)', 'Elbow detected');
    grid on;

end

% function Total_Q(Q_values,K)
% 
%     Q_total = zeros(size(Q_values));
% 
%     for i = 1:length(Q_values)
%        Q_total(i) = (Q_values(i)*Q_values(i))/K(i);
%     end
% 
%     % Affichage des résultats avec deux axes Y
%     figure;
%     semilogx(K, Q_total, 'b', 'LineWidth', 2);
%     ylabel('Facteur de Qualité Total Q');
%     ylim([0 max(Q_values)*1.1]); % Ajustement de l'échelle
%     xlabel('Valeur de K');
%     grid on;
% 
% end

% function Efficiency(C, L, R_par, frequency,k)
%     Q_values = (2 * pi * frequency .* L) ./ R_par; 
%     efficiency = k.^2*Q_values.^2 ./ (k.^2*Q_values.^2 + 1); 
% 
%     figure;
%     semilogx(frequency, efficiency, 'r', 'LineWidth', 2);
% 
%     xlabel('Fréquence (Hz)');
%     ylabel('Efficacité η');
%     title(sprintf('l''efficacité en fct de la fréquence pour C = %.3e F', C));
%     grid on;
% end


%pour Rx: [d_min,d_max]=[3cm,5cm] -> [r_min,r_max] = [1.5cm,2.5cm] 
%pour Tx: [d_min,d_max]=[5cm,8cm] -> [r_min,r_max] = [5cm,8cm]
% check standard h avec le manufact.

function [valid_w, valid_d_out,valid_d_in, valid_N, valid_s,valid_resistance, valid_long_coil] = Planar_Coil(L_target)
% Based on the value we give to the function + the different constraints we
% have for the geometry, we calculate different geometry values to match
% the inductance. A table gives differents possibilities. Also possible to
% see it visually 


% s: Spacing of turns (m)
% N: Number of turns
% d_in: Inner diameter (m)
% d_out: Outer diameter (m)
% w: Wire width, paramètre determiné par PCBWays ? selon PCBWays pour du copper 
% h =[35, 455]um, on prend une valeur fixe de 35um

    % Constraint values for the different parameters
    w_min = 0.1e-3;
    % w_min = 0.5e-3;   
    w_max = 1e-3;   
    d_out_min = 2e-2; 
    d_out_max = 7e-2; 
    s_min = 50e-6;    
    s_max = 500e-6;
    % s_min = 0.2e-3;%0.5e-3;    
    % s_max = 1e-3;

    % w_min = 35e-6;    
    % w_max = 455e-6;   
    % d_out_min = 3e-2; 
    % d_out_max = 5e-2; 
    % s_min = 50e-6;    
    % s_max = 500e-6;  
    % % N_min = 1;        
    % % N_max = 15; 

    N_min = 3;        
    N_max = 20;    
    
    valid_w = [];
    valid_d_out = [];
    valid_d_in = [];
    valid_N = [];
    valid_s = [];
    valid_resistance = [];
    valid_long_coil = [];
 
    for w = linspace(w_min, w_max, 20)  
        for d_out = linspace(d_out_min, d_out_max, 100) 
            for s = linspace(s_min, s_max, 20) 
                for N = N_min:N_max  

                    C_1 = 1.0;
                    C_2 = 2.46;
                    C_3 = 0.0;
                    C_4 = 0.2;
                    


                    % d_in = d_out - 2 * N * (w + s);
                    % if d_in <= 0
                    %     continue;
                    % end
                    d_in = d_out - 2*N*w - 2*(N-1)*s;
                    if d_out > d_in
                        d_avg = (d_out + d_in)/2;
                        kho = (d_out - d_in)/(d_out + d_in);
                        mu = 1.256*10^(-6); % H/m 
    
                        % Make sure the spire don't go beyond d_out
                        % if N * (w + s) >= d_out
                        %     continue; % On ignore cette configuration impossible
                        % end
                        % 
                        % Inductance calculation
                        % L = (N^2 * (d_out - N * (w + s))^2) / (16 * d_out + 28 * N * (w + s)) * 39.38e-6;
    
                        
                            L = mu * N^2 * d_avg *C_1 *(log(C_2/kho) + C_3*kho + C_4*kho^2)/2;
                    end
                    
                    % Check if the inductance is closed to the target
                    if abs(L - L_target) < 0.1e-9  % Tolérance de 0.1 µH
                        valid_w = [valid_w, w];
                        valid_d_out = [valid_d_out, d_out];
                        valid_d_in = [valid_d_in, d_in];
                        valid_N = [valid_N, N];
                        valid_s = [valid_s, s];

                        [R_parasite,long_coil] = Res_para(w, N, d_out, d_in);
                        valid_resistance = [valid_resistance, R_parasite];
                        valid_long_coil = [valid_long_coil, long_coil];

                    end
                end
            end
        end
    end
    
    % Plotting
    figure;
  
    color_resistance = (valid_resistance - min(valid_resistance)) / (max(valid_resistance) - min(valid_resistance));

    scatter3(valid_w * 1e6, valid_d_out * 100, valid_N, 50, color_resistance, 'filled');
    xlabel('Coil width w (\mum)');
    ylabel('Outter diameter d_{out} (cm)');
    zlabel('Nomber of spire N');
    cb = colorbar;
    ylabel(cb, 'Parasitic resistance (Ohms)');
    title(['Valide configurations for L = ', num2str(L_target * 1e6), ' \muH']);
    grid on;
end


function  [resistance_para, l_est] = Res_para(w, N, d_out, d_in)
    kho = 1.68e-8;
    h = 35e-6;  %Copper thickness in the PCB
    section = w*h;
    l_est = pi*N*(d_out+d_in)/2;

    resistance_para = kho*l_est/section;
end

function Affichage(valid_w, valid_d_out, valid_d_in, valid_N, valid_s, valid_resistance, L_coude, longueur_coil)

    fprintf('\n----------------------\n');
    fprintf('Table of the different geometrical parameters of the ''inductor :\n');
    fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------\n');
    fprintf('| %15s | %10s | %15s | %15s | %10s | %15s | %15s | %15s |\n', 'L(µH)','w (µm)', 'd_out (cm)', 'd_in (cm)', 'N', 's (mm)', 'R (Ω)', 'Longueur Coil (cm)');
    fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------\n');

    for i = 1:length(valid_w)
         fprintf('| %15.2f | %10.2f | %15.2f | %15.2f | %10d | %15.2f | %15.6f | %15.6f |\n', ...
            L_coude* 1e6,valid_w(i) * 1e6, valid_d_out(i) * 100, valid_d_in(i) * 100, valid_N(i), valid_s(i) * 1e3, valid_resistance(i), longueur_coil(i)* 100);
    end

fprintf('----------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
end


% function Affichage(valid_w, valid_d_out, valid_d_in, valid_N, valid_s, valid_resistance, L_coude)
% 
%     fprintf('\n----------------------\n');
%     fprintf('Table of the different geometrical parameters of the ''inductor :\n');
%     fprintf('--------------------------------------------------------------------------------------------------------------------------\n');
%     fprintf('| %15s | %10s | %15s | %15s | %10s | %15s | %15s |\n', 'L(µH)','w (µm)', 'd_out (cm)', 'd_in (cm)', 'N', 's (µm)', 'R (Ω)');
%     fprintf('--------------------------------------------------------------------------------------------------------------------------\n');
% 
%     for i = 1:length(valid_w)
%          fprintf('| %15.2f | %10.2f | %15.2f | %15.2f | %10d | %15.2f | %15.6f |\n', ...
%             L_coude* 1e6,valid_w(i) * 1e6, valid_d_out(i) * 100, valid_d_in(i) * 100, valid_N(i), valid_s(i) * 1e6, valid_resistance(i));
%     end
% 
% fprintf('--------------------------------------------------------------------------------------------------------------------------\n');
% end
