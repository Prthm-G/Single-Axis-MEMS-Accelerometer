% Material Properties
E_silicon = 170e9; % Young's modulus of silicon in N/m^2 (170 GPa)
t = 25e-6; % Thickness of the beam in meters (25 µm)
n_eff = 18.5e-6;
finger_dim = [850e-6, 7e-6]; % meters
comb_l = finger_dim(1);
comb_w = finger_dim(2);
N=70;


% Adjustable Beam Dimensions
w_s = 15e-6; % Width of the beam in meters
L_s1 = 750e-6;
L_s2 = 832.5e-6;% Length of the beam in meters

% Spring Constant Calculation using Euler–Bernoulli beam theory
mems_k = 2*((pi^4 / 6) * (E_silicon * t * w_s^3 / ((2*L_s1)^3+(2*L_s1)^3)));

% Display the Spring Constant
fprintf('Beam Width (w_s): %.2e m\n', w_s);
fprintf('Beam Length (L_s): %.2e m\n', L_s1);
fprintf('Calculated Spring Constant (k): %.2e N/m\n', mems_k);

% Proof Mass Calculation
length_m = 3.5*1e-3; % length in m
width = 1.9*1e-3; % width in mm
thickness = 25; % thickness in um
density_si = 2330; % density of silicon in kg/m^3

mems_mass = (calculate_proof_mass(length_m, width, thickness, density_si, N, comb_l,comb_w));
%mass = 4.305e-7;
fprintf('The proof mass is: %.2e kg\n', mems_mass);

% Natural Frequency Calculation
fn = (1/(2*pi)) * sqrt(mems_k/mems_mass);
fprintf('The natural frequency is: %f Hz\n', fn);

% Constants for x_min Calculation
g = 9.81; % acceleration due to gravity in m/s^2

% Minimum Acceleration (0.1% of full scale 80g)
a_min = 0.08 * g;

% Minimum Displacement Calculation
xmin = a_min / (2*pi*fn)^2;

% Convert x_min to micrometers (µm)
xmin_um = xmin * 1e6;

fprintf('The minimum displacement is: %f micrometers\n', xmin_um);

% Constants for Capacitive Sensing
e = 8.854e-12; % permittivity of air in F/m
l = 810e-6; % overlap length of combs in meters (example value)
t = 25e-6; % thickness of the beam in meters (25 µm)
d1 = 13e-6; % smaller gap between electrodes in meters (example value)
d2 = 23e-6; % larger gap between electrodes in meters (example value)
wc = 7e-6; % width of the comb in meters (example value)
    
% Minimum Length of the Proofmass Calculation
%lm = (((d2 + d1 + 2 * wc) * N)-(d2+d1+wc))*1e3;

% Differential Capacitance Calculation
DC = ((2 * e * l * t * N * ((d2^2 - d1^2)*xmin))/((d1^2-xmin^2)*(d2^2-xmin^2)))*1e12;

% Display Capacitive Sensing Results
fprintf('Differential Capacitance (DC): %f pF\n', DC);
%fprintf('Number of Combs (N): %f\n', N);
%fprintf('Minimum Length of the Proofmass (l_m): %f mm\n', lm);
A = l * t;

% Calculate Capacitances for d1 and d2
C1 = e * A*N / d1; % Capacitance with smaller gap (d1)
C2 = e * A*N / d2; % Capacitance with larger gap (d2)

% Convert C1 and C2 to picofarads (pF)
C1_pF = C1 * 1e12;
C2_pF = C2 * 1e12;

% Display Capacitance Values
fprintf('Capacitance with smaller gap (C1): %f pF\n', C1_pF);
fprintf('Capacitance with larger gap (C2): %f pF\n', C2_pF);

b1 = 7.2 * n_eff * t * (l / d1)^3;
b2 = 7.2 * n_eff * t * (l / d2)^3;

mems_b = b1+b2;

zeta = (mems_b / (2 * sqrt(mems_mass * mems_k)));

Q=1/(2*zeta);

fprintf('Damping Coefficient: %.2e Ns/m\n', mems_b);
fprintf('Damping Ratio: %.5f\n', zeta);
fprintf('Quality Factor: %.5f\n', Q);

max_acc = 9.81;
force_at_40g = mems_mass*max_acc;
max_displacement = (force_at_40g / mems_k);

fprintf('Force: %.2e N\n', force_at_40g);


fprintf('Max Displacement: %.2e m\n', max_displacement);

% Given values for FEA
F = 4.14e-6; % Force in Newtons
d_mm = 4.359e-5; % Displacement in millimeters

% Convert displacement from mm to m
d_m = d_mm * 0.001; % 1 mm = 0.001 m

% Calculate spring constant k
k_1 = F / d_m;

% Display the result
disp(['Spring constant for FEA k is: ', num2str(k_1), ' N/m'])

% Function Definitions
function proof_mass = calculate_proof_mass(length_m, width, thickness, density_si,N,comb_l,comb_w)
    % Convert length, width, and thickness to meters
    thickness_m = thickness * 1e-6;

    % Calculate the volume
    volume = ((length_m * width)+(2*N*comb_l*comb_w)) * thickness_m;

    % Calculate the proof mass
    proof_mass = density_si * volume;
end
