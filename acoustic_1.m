clear; clc; close all;

% Parameters for the first plot (Sound Transmission Loss)
freq1 = linspace(200, 500, 1000); % Frequency range
omega = 2 * pi * freq1;
S = 0.0079; % Surface area (m^2)
rho_cl = 1240; % Density of cantilever beam material (kg/m^3)
rho_w = 1000; % Density of surrounding medium
d_cl = 0.0009; % Thickness of cantilever beam (m)
m_cl = 0.0014; % Moving mass of the resonator (kg)
f_cl = 310; % Mechanical resonance frequency (Hz)
omega_cl = 2 * pi * f_cl;
xi_cl = 2.6/100; % Damping ratio for cantilever
j = sqrt(-1); % Imaginary unit

% STL (Sound Transmission Loss) calculation
m_avg = (rho_cl + (m_cl/(S*d_cl)) * ((2*j*xi_cl*omega/omega_cl + 1) ./ ...
    (1 + 2*j*xi_cl*omega/omega_cl - (omega/omega_cl).^2))) / rho_w;

% Create a figure with 3 subplots
figure;

% First subplot: Sound Transmission Loss
subplot(3, 1, 1);
plot(freq1, real(m_avg), 'b-', 'LineWidth', 2); % Plot in blue
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('\rho_{eff}/\rho_{w}', 'FontSize', 12);
grid on;
yline(0, 'k--', 'LineWidth', 1); % Add horizontal line at y = 0

% Parameters for the second plot (Effective Density Ratio)
freq2 = linspace(250, 500); % Frequency range for the second plot
S = 0.0079;
rho_cl = 1240;
d_cl = 0.0009;
m_cl = 0.0014;
f_cl = 310;
xi_cl = 0.26;
d_hr = 0.35;
phi_hr = 0.034;
phi_w = 0.048;
f_hr = 414;
xi_hr = 0.0064;

% STL (Effective Density Ratio) calculation
rho_eff_hr_by_rho_o = (2 + phi_hr) / (2 * (1 - phi_hr));

% Second subplot: Effective Density Ratio
subplot(3, 1, 2);
y2 = ones(size(freq2)) * rho_eff_hr_by_rho_o; % Create a constant vector for plotting
plot(freq2, y2, 'r-', 'LineWidth', 2); % Plot in red
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('\rho_{eff/hr}/\rho_{0}', 'FontSize', 12);
grid on;
yline(0, 'k--', 'LineWidth', 1); % Add horizontal line at y = 0

% Parameters for the third plot (Effective Bulk Modulus)
freq3 = linspace(250, 500); % Frequency range for the third plot
phi_w = 0.048; % Re-defining as required
omega_hr = freq3 / f_hr; % Calculate omega_hr

% STL (Effective Bulk Modulus) calculation
K_eff_hr_by_K_o = 1 ./ (1 - phi_w + (phi_hr ./ (1 + 2*j*xi_hr*omega_hr - omega_hr.^2)));

% Third subplot: Absolute value of K_eff_hr_by_K_o
subplot(3, 1, 3);
y3 = real(K_eff_hr_by_K_o); % Calculate the absolute value
plot(freq3, y3, 'g-', 'LineWidth', 2); % Plot in green
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('(K_{eff/hr}/K_{0})', 'FontSize', 12);
grid on;
yline(0, 'k--', 'LineWidth', 1); % Add horizontal line at y = 0

% Adjust layout
sgtitle('Acoustic Properties vs Frequency', 'FontSize', 16); % Super title for the figure
