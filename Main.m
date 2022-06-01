%% AS5213 UAV Design
% Author: Team SARAS
% Last Modified: 25/ 5/ 2022


%% INTIAL SETTINGS
clc; clear; close all;
PS = PLOT_STANDARDS();


%% GENERAL CONSTANTS AND MISSION SPECIFICATIONS

% Earth gravitational acceleration
global g
g = 9.81;

% Design quantities
global V_cruise V_takeoff V_descent Max_range Max_altitude Max_cruise_height
V_cruise = 20;
V_takeoff = 5;
V_descent = 2.5;
% Range for IITM is 42 km
Max_range = 110 * 1000 * (1 + 5/100);
Max_altitude = 4200;
Max_cruise_height = 45.81;

% Number of motors
global VTOL_motor_count FixedWing_motor_count
VTOL_motor_count = 4;
FixedWing_motor_count = 1;


%% FIRST WEIGHT ESTIMATE

%==================================================
% PAYLOAD WEIGHT ESTIMATE

global payload_redundancy_ratio W_payload W_payload_dropped
payload_redundancy_ratio = 5 / 100;
W_payload = 1.309 * (1 + payload_redundancy_ratio);
W_payload_dropped = 0;
% W_payload = 1.109 * (1 + payload_redundancy_ratio);
% W_payload_dropped = 0.495;


%==================================================
% IMPORT HISTORICAL DATA FROM GOOGLE SHEETS

ID = '1g2JynuGUuHT1ay0MLLTOBL6tM6G6r_BdMQI8w2SrnpY';
sheet_name = 'HistoricData';
url_name = sprintf('https://docs.google.com/spreadsheets/d/%s/gviz/tq?tqx=out:csv&sheet=%s', ID, sheet_name);
HistoricData = webread(url_name);


%==================================================
% STORE THE DATA IN VARIABLES

range_historic = HistoricData.Range;
endurance_historic = HistoricData.Endurance;
cruise_speed_historic = HistoricData.CruiseSpeed;
MTOW_historic = HistoricData.MTOW;
payload_weight_historic = HistoricData.PayloadWeight;
empty_weight_historic = HistoricData.EmptyWeight;
battery_weight_historic = HistoricData.BatteryWeight;


%==================================================
% GET THE EMPTY WEIGHT FRACTIONS AND FIT FUNCTIONAL RELATION WITH MTOW

empty_weight_fraction_historic = empty_weight_historic ./ MTOW_historic;
global battery_redundancy_ratio
battery_redundancy_ratio = 5 / 100;

% Fit function to the data, fit power law relation
% Fitting We/W0 = A*(W0)^c
% Or log(We/W0) = c*log(W0) + log(A)
x = MTOW_historic;
y = empty_weight_fraction_historic;
f = fit(x, y, 'power1');

A = f.a;
c = f.b;

% Plot the Original Data and the Curve Fit
figure(1);
fig1_comps.fig = gcf;
% Plot the data
hold on
fig1_comps.p1 = plot(linspace(min(x), max(x), 100), f(linspace(min(x), max(x), 100)), 'DisplayName', sprintf("$$%.4fW_{0}^{%.4f}$$", A, c), 'LineWidth', 1.25, 'Color', PS.Green2);
fig1_comps.p2 = plot(x, y, 'DisplayName', 'Historic Data', 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', PS.Blue2, 'MarkerEdgeColor', PS.DBlue3);

% Set properties
xlabel('$$W_{0}$$');
ylabel('$$\frac{W_{e}}{W_{0}}$$');
legend();
legendX = 0.74; legendY = 0.8; legendWidth = .1; legendHeight = .1;
fig1_comps.legendPosition = [legendX, legendY, legendWidth, legendHeight];


%==================================================
% ESTIMATE BATTERY WEIGHT FRACTION

battery_weight_fraction_historic = battery_weight_historic ./ MTOW_historic;
average_battery_weight_fraction_historic = mean(battery_weight_fraction_historic);

% Add redundancy in battery weight fraction
battery_redundancy_ratio = 5/100;
battery_weight_fraction = average_battery_weight_fraction_historic;

figure(2);
fig2_comps.fig = gcf;
fig2_comps.p1 = scatter(MTOW_historic, battery_weight_fraction_historic, 'Marker', 'o', 'SizeData', 30, 'MarkerFaceColor', PS.Purple1, 'MarkerEdgeColor', PS.Purple2);

% Set properties
xlabel('$$W_{0}$$');
ylabel('$$\frac{W_{b}}{W_{0}}$$');
legendX = 0.65; legendY = 0.85; legendWidth = .1; legendHeight = .1;
fig2_comps.legendPosition = [legendX, legendY, legendWidth, legendHeight];


%==================================================
% SOLVE TO GET TOTAL WEIGHT

syms W0

eqn = W0 == (W_payload) / (1 - (battery_weight_fraction) - (A*(W0^(c))) );

W0_FirstEstimate = double(vpasolve(eqn, W0));

W_First_EmptyWeight = W0_FirstEstimate * f(W0_FirstEstimate);
W_First_BatteryWeight = W0_FirstEstimate * battery_weight_fraction;
%Can we use the iteration method, so that I can put it in the report?

%========================================================
% PRINT RESULTS OF THIS SECTION

disp('FIRST WEIGHT ESTIMATE');
disp('--------------------');

fprintf('Payload Weight = %.2f kg\n', W_payload);
fprintf('Empty Weight Estimate = %.2f kg\n', W_First_EmptyWeight);
fprintf('Battery Weight Estimate = %.2f kg\n', W_First_BatteryWeight);
fprintf('First Weight Estimate = %.2f kg\n', W0_FirstEstimate);

fprintf('\n\n');


%========================================================
% SAVE FIGURE AS AN IMAGE

figure_location = 'Figures';
SAVE_MY_FIGURE(fig1_comps, sprintf('%s\\EmptyWeightFraction_vs_MTOW.png', figure_location), 'small');
SAVE_MY_FIGURE(fig2_comps, sprintf('%s\\BatteryWeightFraction_vs_MTOW.png', figure_location), 'small');

%========================================================
% CLOSE ALL FIGURES
close all;


%% SECOND WEIGHT ESTIMATE

W_payload;
W_payload_dropped;
W_EmptyWeight = W_First_EmptyWeight;
W_BatteryWeight = W_First_BatteryWeight;
W0_FirstEstimate;

% Empty weight fraction - A * W0^c
A;
c;


%==================================================
% BATTERY SPECIFICATIONS

% We choose Battery with 22.2V
% Commonly available Li-ion battery energy density 300W-hr/kg (without fireproofing)
global Battery_specific_energy battery_SOH battery_SOC battery_PIF battery_discharge_efficiency battery_Voltage
Battery_specific_energy = 300 * 3600;
battery_SOH = 0.9;
battery_SOC = 0.85;
battery_PIF = 0.65;
battery_discharge_efficiency = 0.9;
battery_Voltage = 22.2;

%==================================================
% DESIGN, GEOMETRIC AND ATMOSPHERIC QUANTITIES OF WING AND PROPELLER

% Geometric quantities - taken from similar aircraft
Wing_span = 2.53;
Root_chord = 0.3044;
Tip_chord = 0.2150;
% Tip_chord = 0.1993;
Mean_chord = 0.2623;
% Mean_chord = 0.259701;
% Mean_chord = 0.255495;
% P330_MTOW = 14;
% P330_WingArea = 0.63716;
% P330_WingLoading = (P330_MTOW) / P330_WingArea; 
% Wing_area = W0_FirstEstimate / P330_WingLoading;
Wing_area = 0.63716;
AR = 10.046;
Wing_area = 0.8;
AR = 7.2;
Taper_ratio = 0.7064;
%Taper_ratio = 0.654768

% Propeller - 16 inch diameter - taken from similar aircraft
prop_diameter = (16 * 2.54) / 100;
A_prop = pi * (prop_diameter^2) / 4;

% Atmospheric quantities
std_atm = standardAtmosphere();
Density_sea_level = std_atm.density(0 / 1000);
Density_max_altitude = std_atm.density(Max_altitude / 1000);


%==================================================
% EFFICIENCY FACTORS

global prop_efficiency motor_efficiency ESC_efficiency oswald_efficiency figure_of_merit_hoverpower
prop_efficiency = 0.85;
motor_efficiency = 0.9;
ESC_efficiency = 0.85;
oswald_efficiency = 0.75;
figure_of_merit_hoverpower = 0.8;


%==================================================
% TIME REQUIRED IN DIFFERENT MISSION SEGMENTS

global TakeOff_time Cruise_time Descent_time
TakeOff_time = Max_cruise_height / V_takeoff;
Cruise_time = Max_range / V_cruise;
Descent_time = Max_cruise_height / V_descent;


%==================================================
% ITERATIVELY GET BETTER ESTIMATES OF WEIGHT AT DIFFERENT ALTITUDES

N_h = 100;
h_min = 0;
h_max = Max_altitude;
h_list = linspace(h_min, h_max, N_h);

% Different weights at all altitudes
W_empty_weight_list = cell(N_h, 1);
W_battery_weight_list = cell(N_h, 1);
W_total_weight_list = cell(N_h, 1);

W_Second_EmptyWeight_list = zeros(N_h, 1);
W_Second_BatteryWeight_list = zeros(N_h, 1);
W0_SecondEstimate_list = zeros(N_h, 1);

% Set first element as First Weight Estimates
for i_h = 1:N_h
    W_empty_weight_list{i_h} = [W_EmptyWeight];
    W_battery_weight_list{i_h} = [W_BatteryWeight];
    W_total_weight_list{i_h} = [W0_FirstEstimate];
end

% Find Second Weight Estimate for all the altitudes
for i_h = 1:N_h
    
    iteration_counter = 1;
    % tolerance of 1e-3 corresponds to error of 1 gm
    tolerance = 1e-3;
    error_weight_estimate = 10;
    
    while error_weight_estimate > tolerance
        
        h = h_list(i_h);
        rho_h = std_atm.density(h / 1000);
        
        %==================================================
        % AERODYNAMIC QUANTITIES

        % Design lift coefficient estimation
        % In cruise W = L = (1/2) * rho * V^2 * S * Cl

        CL_cruise = (g * W_total_weight_list{i_h}(iteration_counter)) / (0.5 * rho_h * V_cruise * V_cruise * Wing_area);
        CD0 = CD0_func(Wing_area);
        CD_cruise = CD0 + (1 / (pi * AR * oswald_efficiency)) * (CL_cruise^2);
        

        %==================================================
        % POWER REQUIRED CALCULATIONS

        % Cruise
        T_cruise = 0.5 * rho_h * V_cruise * V_cruise * Wing_area * CD_cruise;
        P_cruise = V_cruise * T_cruise;

        % Vertical take off 1
        K_T = 1.2;
        T_takeoff_1 = K_T * (g * W_total_weight_list{i_h}(iteration_counter) / VTOL_motor_count);
        P_takeoff_1 = VTOL_motor_count * (V_takeoff * T_takeoff_1 / 2) * (1 + sqrt(1 + ((2 * T_takeoff_1) / (rho_h * (V_takeoff^2) * A_prop))));

        % Vertical descent 1
        T_hover_1 = g * W_total_weight_list{i_h}(iteration_counter) / VTOL_motor_count;
        V_propeller_hover = sqrt(T_hover_1 / (2 * rho_h * A_prop));
        x = -1 * (V_descent / V_propeller_hover);
        K = 1.2;
        V_propeller_induced = (K - 1.125*x - 1.372*(x^2) - 1.718*(x^3) - 0.655*(x^4)) * V_propeller_hover;
        P_descent_1 = VTOL_motor_count * K * (g * W_total_weight_list{i_h}(iteration_counter) / VTOL_motor_count) * (V_propeller_induced - V_descent);

        % Vertical take off 2
        K_T = 1.2;
        T_takeoff_2 = K_T * (g * (W_total_weight_list{i_h}(iteration_counter) - W_payload_dropped) / VTOL_motor_count);
        P_takeoff_2 = VTOL_motor_count * (V_takeoff * T_takeoff_2 / 2) * (1 + sqrt(1 + ((2 * T_takeoff_2) / (rho_h * (V_takeoff^2) * A_prop))));

        % Vertical descent 2
        T_hover_2 = g * (W_total_weight_list{i_h}(iteration_counter) - W_payload_dropped) / VTOL_motor_count;
        V_propeller_hover = sqrt(T_hover_2 / (2 * rho_h * A_prop));
        x = -1 * (V_descent / V_propeller_hover);
        K = 1.2;
        V_propeller_induced = (K - 1.125*x - 1.372*(x^2) - 1.718*(x^3) - 0.655*(x^4)) * V_propeller_hover;
        P_descent_2 = VTOL_motor_count * K * (g * (W_total_weight_list{i_h}(iteration_counter) - W_payload_dropped) / VTOL_motor_count) * (V_propeller_induced - V_descent);


        %==================================================
        % POWER REQUIRED ACCOUNTING FOR INEFFICIENCIES

        P_cruise = P_cruise / (prop_efficiency * motor_efficiency * ESC_efficiency);
        P_takeoff_1 = P_takeoff_1 / (figure_of_merit_hoverpower * motor_efficiency * ESC_efficiency);
        P_descent_1 = P_descent_1 / (figure_of_merit_hoverpower * motor_efficiency * ESC_efficiency);
        P_takeoff_2 = P_takeoff_2 / (figure_of_merit_hoverpower * motor_efficiency * ESC_efficiency);
        P_descent_2 = P_descent_2 / (figure_of_merit_hoverpower * motor_efficiency * ESC_efficiency);
        
        
        %==================================================
        % ENERGY REQUIRED CALCULATIONS

        Energy_cruise = P_cruise * Cruise_time;
        Energy_takeoff_1 = P_takeoff_1 * TakeOff_time;
        Energy_descent_1 = P_descent_1 * Descent_time;
        Energy_takeoff_2 = P_takeoff_2 * TakeOff_time;
        Energy_descent_2 = P_descent_2 * Descent_time;
        

        %==================================================
        % BATTERY WEIGHT REQUIRED CALCULATIONS

        Total_energy_required = Energy_cruise + Energy_takeoff_1 + Energy_descent_1 + Energy_takeoff_2 + Energy_descent_2;

        % Accounting for inefficiencies
        Total_energy_required = Total_energy_required / (battery_SOH * battery_SOC * battery_discharge_efficiency);
        Battery_weight = (Total_energy_required / Battery_specific_energy) / battery_PIF;
        Battery_weight = Battery_weight * (1 + battery_redundancy_ratio);
        
        
        %==================================================
        % WEIGHT ESTIMATE

        syms W0

        eqn = W0 == (W_payload + Battery_weight) / (1 - (A*(W0^(c))) );

        % Get new weight estimate
        iteration_counter = iteration_counter + 1;
        W_total_weight_list{i_h}(iteration_counter) = double(vpasolve(eqn, W0));
        % W_total_weight_list{i_h}(iteration_counter) = double(solve(eqn, W0));

        error_weight_estimate = abs(W_total_weight_list{i_h}(iteration_counter) - W_total_weight_list{i_h}(iteration_counter - 1));

        W_empty_weight_list{i_h}(iteration_counter) = A*(W_total_weight_list{i_h}(iteration_counter)^c) * W_total_weight_list{i_h}(iteration_counter);
        W_battery_weight_list{i_h}(iteration_counter) = Battery_weight;

        % iteration_counter


    end


    %==================================================
    % ACCOUNT FOR REDUNDANCIES

    W_Second_EmptyWeight_list(i_h) = W_empty_weight_list{i_h}(iteration_counter);
    W_Second_BatteryWeight_list(i_h) = W_battery_weight_list{i_h}(iteration_counter);
    W0_SecondEstimate_list(i_h) = W_payload + W_Second_EmptyWeight_list(i_h) + W_Second_BatteryWeight_list(i_h);

end


%==================================================
% FIND THE MAXIMUM SECOND WEIGHT ESTIMATE AND PARAMETER CALCULATIONS

[W0_SecondEstimate_max, idx_W0_SecondEstimate_max] = max(W0_SecondEstimate_list);

W_Second_EmptyWeight = W_Second_EmptyWeight_list(idx_W0_SecondEstimate_max);
W_Second_BatteryWeight = W_Second_BatteryWeight_list(idx_W0_SecondEstimate_max);
W0_SecondEstimate = W0_SecondEstimate_list(idx_W0_SecondEstimate_max);

h_W0_SecondEstimate_max = h_list(idx_W0_SecondEstimate_max);
rho_h = std_atm.density(h_W0_SecondEstimate_max / 1000);
CL_cruise = (g * W0_SecondEstimate) / (0.5 * rho_h * V_cruise * V_cruise * Wing_area);
CD0 = CD0_func(Wing_area);
CD_cruise = CD0 + (1 / (pi * AR * oswald_efficiency)) * (CL_cruise^2);

[P_cruise, P_takeoff_1, P_descent_1, P_takeoff_2, P_descent_2, Total_energy_required, Battery_weight] = Power_Calculation_func(h_W0_SecondEstimate_max, W0_SecondEstimate, Wing_area, AR, A_prop);
Energy_cruise = P_cruise * Cruise_time;
Energy_takeoff_1 = P_takeoff_1 * TakeOff_time;
Energy_descent_1 = P_descent_1 * Descent_time;
Energy_takeoff_2 = P_takeoff_2 * TakeOff_time;
Energy_descent_2 = P_descent_2 * Descent_time;


%==================================================
% PLOT FIGURES FOR WEIGHT ESTIMATES THROUGH THE ITERATIONS

% Show the convergence of the total weight estimate
figure(1);
fig1_comps.fig = gcf;
% Plot the data
hold on
fig1_comps.p1 = plot(1:length(W_total_weight_list{idx_W0_SecondEstimate_max}), W_total_weight_list{idx_W0_SecondEstimate_max}, 'LineWidth', 1.25, 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', PS.Blue1, 'MarkerEdgeColor', PS.DBlue2);

% Set properties
xlabel('$$Iteration Count$$');
ylabel('$$W_{0}$$');
% legend();
legendX = 0.78; legendY = 0.77; legendWidth = .1; legendHeight = .1;
fig1_comps.legendPosition = [legendX, legendY, legendWidth, legendHeight];


%==================================================
% PRINT RESULTS OF THIS SECTION

disp('SECOND WEIGHT ESTIMATE');
disp('--------------------');

fprintf('TakeOff_time = %.2f s or %.2f min\n', TakeOff_time, TakeOff_time / 60);
fprintf('Cruise_time = %.2f s or %.2f min\n', Cruise_time, Cruise_time / 60);

fprintf('Height Second Weight Estimate max = %.2f\n', h_W0_SecondEstimate_max);

fprintf('CL_cruise = %.2f\n', CL_cruise);
fprintf('CD_cruise = %.2f\n', CD_cruise);

fprintf('P_cruise = %.2f Watts\n', P_cruise);
fprintf('P_takeoff_1 = %.2f Watts\n', P_takeoff_1);
fprintf('P_descent_1 = %.2f Watts\n', P_descent_1);
fprintf('P_takeoff_2 = %.2f Watts\n', P_takeoff_2);
fprintf('P_descent_2 = %.2f Watts\n', P_descent_2);

fprintf('Energy_cruise = %.2f (%.4e) J\n', Energy_cruise, Energy_cruise);
fprintf('Energy_takeoff_1 = %.2f (%.4e) J\n', Energy_takeoff_1, Energy_takeoff_1);
fprintf('Energy_descent_1 = %.2f (%.4e) J\n', Energy_descent_1, Energy_descent_1);
fprintf('Energy_takeoff_2 = %.2f (%.4e) J\n', Energy_takeoff_2, Energy_takeoff_2);
fprintf('Energy_descent_2 = %.2f (%.4e) J\n', Energy_descent_2, Energy_descent_2);

fprintf('Total_energy_required = %.2f J\n', Total_energy_required);
fprintf('Battery_weight = %.2f kg\n', Battery_weight);

fprintf('Payload Weight = %.2f kg\n', W_payload);
fprintf('Empty Weight Estimate = %.2f kg\n', W_Second_EmptyWeight);
fprintf('Battery Weight Estimate = %.2f kg\n', W_Second_BatteryWeight);
fprintf('Final Weight Estimate = %.2f kg\n', W0_SecondEstimate);

fprintf('\n\n');


%========================================================
% SAVE FIGURE AS AN IMAGE

figure_location = 'Figures';
SAVE_MY_FIGURE(fig1_comps, sprintf('%s\\SecondWeightEstimate_ConvergencePlot.png', figure_location), 'small');


%========================================================
% CLOSE ALL FIGURES
close all;


%% WING LOADING AND POWER LOADING

%==================================================
% FACTORS GOVERNING WING LOADING
% 1) Prescribed Flight Speed
% 2) Absolute Ceiling
% 3) Range

% POWER LOADING = (P/W)
% GOAL: Minimize Power Loading or Minimize Battery Capacity or Fix Cruise Speed at Absolute Ceiling
% NOTE: Take into account the efficiency factors


%==================================================
% SET REQUIRED DATA VALUES

MTOW = W0_SecondEstimate;

N_h = 100;
h_min = 0;
h_max = Max_altitude;
h_list = linspace(h_min, h_max, N_h);
% Get std_atm object
std_atm = standardAtmosphere();
% Get area list
N_S = 1000;
S_min = 0.25;
S_max = 2;
S_list = linspace(S_min, S_max, N_S);
WL_list = (MTOW * g) ./ S_list;
% Expressions for CD0 and K taken from https://nptel.ac.in/courses/101106035 lecture 8
CD0 = CD0_func(S_list);
K = 1.333 / (pi * AR);

global Variation_percentage
Variation_percentage = 5 / 100;


%==================================================
% PRESCRIBED FLIGHT SPEED
% Idea: Fix V_cruise, Go to each height find the variation of PL vs WL. Find min at each altitude.
% Find max of all the minimas. Take 5% variation of that particular WL. That is our I1.

PL_list = zeros(N_h, N_S);
for i = 1:N_h
    h = h_list(i);
    rho_h = std_atm.density(h / 1000);
    % Get Power Loading
    q = (1/2) * rho_h * V_cruise * V_cruise;
    CL = WL_list / q;
    CD = CD0 + K * (CL.^2);
    D = q .* S_list .* CD;
    T = D;
    PL_list(i, :) = (T * V_cruise) / (MTOW * g);
    % PL_list(i, :) = (V_cruise * q .* CD) ./ (WL_list);
end

% Plot variation of Power Loading with Wing Loading at different altitudes
figure(1);
fig1_comps.fig = gcf;
% Plot the data
hold on
i_list = 0: (N_h / 5): N_h;
i_list(1) = 1;
for i = i_list
    fig1_comps.p1(i) = plot(WL_list, PL_list(i, :), 'DisplayName', sprintf('h = %.2f', h_list(i)), 'LineWidth', 1.25);
end

% Set properties
xlabel('$$Wing Loading \; \frac{W}{S}$$');
ylabel('$$Power Loading \; \frac{P}{W}$$');
legend();
legendX = 0.75; legendY = 0.77; legendWidth = .1; legendHeight = .1;
fig1_comps.legendPosition = [legendX, legendY, legendWidth, legendHeight];

STANDARDIZE_FIGURE(fig1_comps);
figure_location = 'Figures';
SAVE_MY_FIGURE(fig1_comps, sprintf('%s\\PowerLoading_vs_WingLoading_PrescribedVelocity_at_different_Altitudes.png', figure_location), 'small');


% Find Maximum of the Minimum Power Loading at each Each Altitude. Take 5% variation around that as the Interval for Wing Loading.
PL_min_list = zeros(N_h, 1);
for i = 1:N_h
    PL_min_list(i) = min(PL_list(i, :));
end
[max_PL_min_list, idx_max_PL_min_list] = max(PL_min_list);

% Get Power Loading at that Altitude and the Interval of Wing Loading for 5% Variation of Power Loading
h_prescribedV = h_list(idx_max_PL_min_list);
rho_h_prescribedV = std_atm.density(h_prescribedV / 1000);
% Get Power Loading
q = (1/2) * rho_h_prescribedV * V_cruise * V_cruise;
CL = WL_list / q;
CD = CD0 + K * (CL.^2);
D = q .* S_list .* CD;
T = D;
PL = (T * V_cruise) / (MTOW * g);

PL_min = min(PL);
PL_percent_variation_idx = find(abs(PL - PL_min)/PL_min < Variation_percentage);
WL_interval_prescribedV = WL_list(PL_percent_variation_idx);
PL_interval_prescribedV = PL(PL_percent_variation_idx);

% Plot variation of Power Loading with Wing Loading and 5% Variation
figure(2);
fig2_comps.fig = gcf;
% Plot the data
fig2_comps.p1 = plot(WL_list, PL, 'LineWidth', 1.25);
hold on
fig2_comps.p2 = xline(min(WL_interval_prescribedV), '--', 'LineWidth', 1);
fig2_comps.p3 = xline(max(WL_interval_prescribedV), '--', 'LineWidth', 1);
fig2_comps.p4 = yline(min(PL_interval_prescribedV), '--', 'LineWidth', 1);
fig2_comps.p5 = yline(max(PL_interval_prescribedV), '--', 'LineWidth', 1);

% Set properties
xlabel('$$Wing Loading \; \frac{W}{S}$$');
ylabel('$$Power Loading \; \frac{P}{W}$$');

STANDARDIZE_FIGURE(fig2_comps);
figure_location = 'Figures';
SAVE_MY_FIGURE(fig2_comps, sprintf('%s\\PowerLoading_vs_WingLoading_PrescribedVelocity_5PercentVariation.png', figure_location), 'small');


%==================================================
% ABSOLUTE CEILING

rho_ceil = std_atm.density(h_max / 1000);

% The CL value which minimizes Power Required is used
CL = sqrt((3 * CD0) / K);
CD = CD0 + K * CL.^(2);
PL_1 = sqrt((2 * MTOW * g) ./ (rho_ceil * S_list)) .* (CD ./ (CL.^(1.5)));

% IMPORTANT: What was wrong with the code below?
% q = (1/2) * rho_ceil * V_cruise * V_cruise;
% CL = WL_list / q;
% CD = CD0 + K * (CL.^2);
% D = q .* S_list .* CD;
% PL_2 = (D * V_cruise) / (MTOW * g);

q = (1/2) * rho_ceil * V_cruise * V_cruise;
CL = sqrt((3 * CD0) / K);
CD = CD0 + K * CL.^(2);
D = q .* S_list .* CD;
PL_2 = (D * V_cruise) / (MTOW * g);

% Plot figures for variation with Altitude
figure(3);
fig3_comps.fig = gcf;
% Plot the data
fig3_comps.p1 = plot(WL_list, PL_1, 'DisplayName', 'Minimum Power', 'LineWidth', 1.25);
hold on
fig3_comps.p2 = plot(WL_list, PL_2, 'DisplayName', 'Cruise Velocity 20', 'LineWidth', 1.25);
% Set properties
xlabel('$$Wing Loading \; \frac{W}{S}$$');
legend();
legendX = 0.72; legendY = 0.8; legendWidth = .1; legendHeight = .1;
fig3_comps.legendPosition = [legendX, legendY, legendWidth, legendHeight];

STANDARDIZE_FIGURE(fig3_comps);
figure_location = 'Figures';
SAVE_MY_FIGURE(fig3_comps, sprintf('%s\\PowerLoading_vs_WingLoading_at_AbsoluteCeiling.png', figure_location), 'small');

% 5percent Variation
% Power Loading at the point of nearest approach
PL_intersection = Point_of_Closest_Approach(WL_list, PL_1, PL_2);
PL_1_percent_variation_idx = find(abs(PL_1 - PL_intersection)/PL_intersection < Variation_percentage);
PL_2_percent_variation_idx = find(abs(PL_2 - PL_intersection)/PL_intersection < Variation_percentage);
WL_1_interval_absceil = WL_list(PL_1_percent_variation_idx);
WL_2_interval_absceil = WL_list(PL_2_percent_variation_idx);
PL_1_interval_absceil = PL_1(PL_1_percent_variation_idx);
PL_2_interval_absceil = PL_2(PL_2_percent_variation_idx);

% Plot variation of Power Loading with Wing Loading and 5% Variation
figure(4);
fig4_comps.fig = gcf;
% Plot the data
fig4_comps.p1 = plot(WL_list, PL_1, 'DisplayName', 'PL_1', 'LineWidth', 1.25);
hold on
fig4_comps.p2 = plot(WL_list, PL_2, 'DisplayName', 'PL_2', 'LineWidth', 1.25);
fig4_comps.p3(1) = xline(min(WL_1_interval_absceil), '--r', 'LineWidth', 1, 'HandleVisibility', 'off');
fig4_comps.p3(2) = xline(max(WL_1_interval_absceil), '--r', 'LineWidth', 1, 'HandleVisibility', 'off');
fig4_comps.p4(1) = xline(min(WL_2_interval_absceil), '--b', 'LineWidth', 1, 'HandleVisibility', 'off');
fig4_comps.p4(2) = xline(max(WL_2_interval_absceil), '--b', 'LineWidth', 1, 'HandleVisibility', 'off');
fig4_comps.p5(1) = yline(min(PL_1_interval_absceil), '--r', 'LineWidth', 1, 'HandleVisibility', 'off');
fig4_comps.p5(2) = yline(max(PL_1_interval_absceil), '--r', 'LineWidth', 1, 'HandleVisibility', 'off');
fig4_comps.p6(1) = yline(min(PL_2_interval_absceil), '--b', 'LineWidth', 1, 'HandleVisibility', 'off');
fig4_comps.p6(2) = yline(max(PL_2_interval_absceil), '--b', 'LineWidth', 1, 'HandleVisibility', 'off');

% Set properties
xlabel('$$Wing Loading \; \frac{W}{S}$$');
ylabel('$$Power Loading \; \frac{P}{W}$$');
legend();
legendX = 0.78; legendY = 0.77; legendWidth = .1; legendHeight = .1;
fig3_comps.legendPosition = [legendX, legendY, legendWidth, legendHeight];

STANDARDIZE_FIGURE(fig4_comps);
figure_location = 'Figures';
SAVE_MY_FIGURE(fig4_comps, sprintf('%s\\PowerLoading_vs_WingLoading_at_AbsoluteCeiling_5PercentVariation.png', figure_location), 'small');


%==================================================
% BATTERY WEIGHT

C_list = zeros(N_h, N_S);
for i = 1:N_h
    h = h_list(i);
    rho_h = std_atm.density(h / 1000);
    % Get Power Loading
    q = (1/2) * rho_h * V_cruise * V_cruise;
    [~, ~, ~, ~, ~, Total_energy_required, Battery_weight] = Power_Calculation_func(h, MTOW, S_list, AR, A_prop);
    Battery_weight_list(i, :) = Battery_weight;
    % C_list(i, :) = Max_range * (q .* S_list .* CD0 + ((MTOW*g)^2 * K)./(q .* S_list)) / (prop_efficiency * battery_Voltage);
end

% Plot figures for variation with Altitude
figure(5);
fig5_comps.fig = gcf;
% Plot the data
hold on
i_list = 0: (N_h / 5): N_h;
i_list(1) = 1;
for i = i_list
    fig5_comps.p(i) = plot(WL_list, Battery_weight_list(i, :), 'DisplayName', sprintf('h = %.2f', h_list(i)), 'LineWidth', 1.25);
end

% Set properties
xlabel('$$Wing Loading \; \frac{W}{S}$$');
ylabel('$$Battery Weight$$');
legend();
legendX = 0.74; legendY = 0.77; legendWidth = .1; legendHeight = .1;
fig5_comps.legendPosition = [legendX, legendY, legendWidth, legendHeight];

STANDARDIZE_FIGURE(fig5_comps);
figure_location = 'Figures';
SAVE_MY_FIGURE(fig5_comps, sprintf('%s\\BatteryWeight_vs_WingLoading_Range_at_different_Altitudes.png', figure_location), 'small');

% Find Maximum of the Minimum Battery Capacity at each Each Altitude. Take 5% variation around that as the Interval for Wing Loading.
Battery_weight_min_list = zeros(N_h, 1);
for i = 1:N_h
    Battery_weight_min_list(i) = min(Battery_weight_list(i, :));
end
[max_Battery_weight_min_list, idx_max_Battery_weight_min_list] = max(Battery_weight_min_list);

% Get Power Loading at that Altitude and the Interval of Wing Loading for 5% Variation of Power Loading
h_rangeV = h_list(idx_max_Battery_weight_min_list);
rho_h_rangeV = std_atm.density(h_rangeV / 1000);

% Get Battery weight
[~, ~, ~, ~, ~, Total_energy_required, Battery_weight] = Power_Calculation_func(h_rangeV, MTOW, S_list, AR, A_prop);

Battery_weight_min = min(Battery_weight);
Battery_weight_percent_variation_idx = find(abs(Battery_weight - Battery_weight_min)/Battery_weight_min < Variation_percentage);
WL_interval_rangeV = WL_list(Battery_weight_percent_variation_idx);
Battery_weight_interval_rangeV = Battery_weight(Battery_weight_percent_variation_idx);

% Plot variation of Power Loading with Wing Loading and 5% Variation
figure(6);
fig6_comps.fig = gcf;
% Plot the data
fig6_comps.p1 = plot(WL_list, Battery_weight, 'LineWidth', 1.25);
hold on
fig6_comps.p2 = xline(min(WL_interval_rangeV), '--', 'LineWidth', 1);
fig6_comps.p3 = xline(max(WL_interval_rangeV), '--', 'LineWidth', 1);
fig6_comps.p4 = yline(min(Battery_weight_interval_rangeV), '--', 'LineWidth', 1);
fig6_comps.p5 = yline(max(Battery_weight_interval_rangeV), '--', 'LineWidth', 1);

% Set properties
xlabel('$$Wing Loading \; \frac{W}{S}$$');
ylabel('$$Battery Weight$$');

STANDARDIZE_FIGURE(fig6_comps);
figure_location = 'Figures';
SAVE_MY_FIGURE(fig6_comps, sprintf('%s\\BatteryWeight_vs_WingLoading_Range_5PercentVariation.png', figure_location), 'small');


%==================================================
% COMMON INTERVAL FOR WING LOADING THAT SATISFIES ALL THE CONSTRAINTS
WL_interval_min_list = [min(WL_interval_prescribedV), min(WL_1_interval_absceil), min(WL_2_interval_absceil), min(WL_interval_rangeV)];
WL_interval_max_list = [max(WL_interval_prescribedV), max(WL_1_interval_absceil), max(WL_2_interval_absceil), max(WL_interval_rangeV)];
WL_interval_common_min = max(WL_interval_min_list);
WL_interval_common_max = min(WL_interval_max_list);
WL_interval_common = WL_list((WL_list >= WL_interval_common_min) & (WL_list <= WL_interval_common_max));
S_interval_common_min = (MTOW * g) / WL_interval_common_max;
S_interval_common_max = (MTOW * g) / WL_interval_common_min;

figure(7);
fig7_comps.fig = gcf;
% Plot the data
fig7_comps.p1 = plot(WL_interval_prescribedV, 6 * ones(length(WL_interval_prescribedV), 1), 'DisplayName', 'Prescribed Velocity', 'LineWidth', 1.25);
hold on
fig7_comps.p2 = plot(WL_1_interval_absceil, 5 * ones(length(WL_1_interval_absceil), 1), 'DisplayName', 'Absolute Ceiling I1', 'LineWidth', 1.25);
fig7_comps.p3 = plot(WL_2_interval_absceil, 4 * ones(length(WL_2_interval_absceil), 1), 'DisplayName', 'Absolute Ceiling I2', 'LineWidth', 1.25);
fig7_comps.p4 = plot(WL_interval_rangeV, 3 * ones(length(WL_interval_rangeV), 1), 'DisplayName', 'Battery Weight', 'LineWidth', 1.25);
fig7_comps.p5 = plot(WL_interval_common, 1 * ones(length(WL_interval_common), 1), 'DisplayName', 'Common Interval', 'LineWidth', 1.75);

fig7_comps.p6 = xline(min(WL_interval_common), '--', 'LineWidth', 1, 'HandleVisibility', 'off');
fig7_comps.p7 = xline(max(WL_interval_common), '--', 'LineWidth', 1, 'HandleVisibility', 'off');

x = [WL_interval_common_min, WL_interval_common_max, WL_interval_common_max, WL_interval_common_min];
y = [0, 0, 8, 8];
patch(x, y, PS.Blue1, 'HandleVisibility', 'off', 'FaceAlpha', .3);

ylim([0, 8]);
set(gca, 'YTick', []);
legend('Location', 'northwest');

% Set properties
xlabel('$$Wing Loading \; \frac{W}{S}$$');

STANDARDIZE_FIGURE(fig7_comps);
figure_location = 'Figures';
SAVE_MY_FIGURE(fig7_comps, sprintf('%s\\WingLoading_CommonInterval.png', figure_location), 'small');


%==================================================
% PRINT RESULTS OF THIS SECTION

disp('WING LOADING AND POWER LOADING');
disp('--------------------');

fprintf('Wing Loading Min = %.2f\n', WL_interval_common_min);
fprintf('Wing Loading Max = %.2f\n', WL_interval_common_max);

fprintf('MTOW = %.2f\n', MTOW);

fprintf('Wing Area Min = %.2f\n', S_interval_common_min);
fprintf('Wing Area Max = %.2f\n', S_interval_common_max);

fprintf('\n\n');


%========================================================
% CLOSE ALL FIGURES
close all;


%% DISK LOADING AND POWER LOADING

%==================================================
% FACTORS GOVERNING DISK LOADING
% 1) Take-off
% 2) Descent

% DISK LOADING = (MTOW / VTOL_motor_count) / (A_prop)
% POWER LOADING = (P_prop/W) - defined per propeller
% GOAL: Minimize Power Loading
% NOTE: Take into account the efficiency factors


%==================================================
% SET REQUIRED DATA VALUES

MTOW = W0_SecondEstimate;

N_h = 100;
h_min = 0;
h_max = Max_altitude;
h_list = linspace(h_min, h_max, N_h);
% Get std_atm object
std_atm = standardAtmosphere();
% Get area list
N_Aprop = 1000;
% A_prop from past data is 0.1297
Aprop_min = 0.05;
Aprop_max = .3;
Aprop_list = linspace(S_min, S_max, N_S);
DL_list = (MTOW * g / VTOL_motor_count) ./ Aprop_list;
Weight_per_prop = (MTOW * g) / VTOL_motor_count;

% Find out the power required at different disk loading at varying altitudes
PL_takeoff_list = zeros(N_h, N_Aprop);
PL_descent_list = zeros(N_h, N_Aprop);
for i = 1:N_h
    h = h_list(i);
    rho_h = std_atm.density(h / 1000);
    % Get Power required
    % Dummy Wing area - this will not affect power in hover
    Wing_area = 0.6;
    % Only power required in takeoff and descent are of interest
    [~, P_takeoff_1, P_descent_1, ~, ~, ~, ~] = Power_Calculation_func(h, MTOW, Wing_area, AR, Aprop_list);
    PL_takeoff_list(i, :) = (P_takeoff_1 / VTOL_motor_count) / (Weight_per_prop);
    PL_descent_list(i, :) = (P_descent_1 / VTOL_motor_count) / (Weight_per_prop);
end

% Plot variation of Power required takeoff with Disk Loading at different altitudes
figure(1);
fig1_comps.fig = gcf;
% Plot the data
hold on
i_list = 0: (N_h / 5): N_h;
i_list(1) = 1;
for i = i_list
    fig1_comps.p1(i) = plot(DL_list, PL_takeoff_list(i, :), 'DisplayName', sprintf('h = %.2f', h_list(i)), 'LineWidth', 1.25);
end

% Set properties
xlabel('$$Disk Loading \; \frac{W}{S}$$');
ylabel('$$Power Loading (Takeoff)\; \frac{P}{W}$$');
legend();
legendX = 0.75; legendY = 0.77; legendWidth = .1; legendHeight = .1;
fig1_comps.legendPosition = [legendX, legendY, legendWidth, legendHeight];

STANDARDIZE_FIGURE(fig1_comps);
figure_location = 'Figures';
SAVE_MY_FIGURE(fig1_comps, sprintf('%s\\PowerLoadingTakeoff_vs_DiskLoading_at_different_Altitudes.png', figure_location), 'small');

% Plot variation of Power required Descent with Disk Loading at different altitudes
figure(2);
fig2_comps.fig = gcf;
% Plot the data
hold on
i_list = 0: (N_h / 5): N_h;
i_list(1) = 1;
for i = i_list
    fig2_comps.p1(i) = plot(DL_list, PL_descent_list(i, :), 'DisplayName', sprintf('h = %.2f', h_list(i)), 'LineWidth', 1.25);
end

% Set properties
xlabel('$$Disk Loading \; \frac{W}{S}$$');
ylabel('$$Power Loading (Takeoff)\; \frac{P}{W}$$');
legend();
legendX = 0.75; legendY = 0.77; legendWidth = .1; legendHeight = .1;
fig2_comps.legendPosition = [legendX, legendY, legendWidth, legendHeight];

STANDARDIZE_FIGURE(fig2_comps);
figure_location = 'Figures';
SAVE_MY_FIGURE(fig2_comps, sprintf('%s\\PowerLoadingDescent_vs_DiskLoading_at_different_Altitudes.png', figure_location), 'small');

% Find Maximum of the Minimum Power Loading at each Each Altitude.
PL_takeoff_min_list = zeros(N_h, 1);
PL_descent_min_list = zeros(N_h, 1);
for i = 1:N_h
    PL_takeoff_min_list(i) = min(PL_takeoff_list(i, :));
    PL_descent_min_list(i) = min(PL_descent_list(i, :));
end
[max_PL_takeoff_min_list, idx_max_PL_takeoff_min_list] = max(PL_takeoff_min_list);
[max_PL_descent_min_list, idx_max_PL_descent_min_list] = max(PL_descent_min_list);

% Plot the variation of power loading with disk loading for take off and descent at this altitude
figure(3);
fig3_comps.fig = gcf;
% Plot the data
fig3_comps.p1 = plot(DL_list, PL_takeoff_list(idx_max_PL_takeoff_min_list, :), 'DisplayName', 'Takeoff', 'LineWidth', 1.25);
hold on
fig3_comps.p2 = plot(DL_list, PL_descent_list(idx_max_PL_descent_min_list, :), 'DisplayName', 'Descent', 'LineWidth', 1.25);

% Set properties
xlabel('$$Disk Loading \; \frac{W}{S}$$');
ylabel('$$Power Loading \; \frac{P}{W}$$');
legend();
legendX = 0.75; legendY = 0.77; legendWidth = .1; legendHeight = .1;
fig2_comps.legendPosition = [legendX, legendY, legendWidth, legendHeight];

STANDARDIZE_FIGURE(fig3_comps);
figure_location = 'Figures';
SAVE_MY_FIGURE(fig3_comps, sprintf('%s\\PowerLoading_vs_DiskLoading_MaxValues.png', figure_location), 'small');


%==================================================
% PRINT RESULTS OF THIS SECTION

disp('DISK LOADING AND POWER LOADING');
disp('--------------------');

fprintf('PL max in takeoff at h = %.2f\n', h_list(idx_max_PL_takeoff_min_list));
fprintf('PL max in descent at h = %.2f\n', h_list(idx_max_PL_descent_min_list));

fprintf('\n\n');


%========================================================
% CLOSE ALL FIGURES
close all;


%% MACH NUMBER AND REYNOLDS NUMBER VARIATION WITH ALTITUDE

%==================================================
% VARIATION OF M AND Re WITH ALTITUDE

std_atm = standardAtmosphere();
N_h = 100;
h_min = 0;
h_max = Max_altitude;
h_list = linspace(h_min, h_max, N_h);

% M = V_cruise / c_h
% Re = ((rho_h * V_cruise * L) / mu_h) - where L is a characteristic length. We choose it to be chord length

rho_h = std_atm.density(h_list / 1000);
T_h = std_atm.temperature(h_list / 1000);
gamma = 1.4;
R = 287.058;
c_h = sqrt(gamma * R * T_h);

mu_h = std_atm.viscosity(h_list / 1000);

% Characteristic Length
L_characteristic = .25;

% Mach Number
M_h = V_cruise ./ c_h;

% Reynolds Number
Re_h = (rho_h * V_cruise * L_characteristic) ./ mu_h;

% Minimum and Maximum Mach No. and Reynolds No.
[M_min, idx_M_min] = min(M_h);
h_M_min = h_list(idx_M_min);
[M_max, idx_M_max] = max(M_h);
h_M_max = h_list(idx_M_max);

[Re_min, idx_Re_min] = min(Re_h);
h_Re_min = h_list(idx_Re_min);
[Re_max, idx_Re_max] = max(Re_h);
h_Re_max = h_list(idx_Re_max);


%==================================================
% PLOT THE VARIATION OF MACH NUMBER AND REYNOLDS NUMBER WITH HEIGHT

% Mach Number
figure(1);
fig1_comps.fig = gcf;
% Plot the data
fig1_comps.p1 = plot(h_list, M_h, 'LineWidth', 1.25);

% Set properties
title('Mach Number variation with Altitude');
xlabel('Altitude h');
ylabel('Mach Number M');

figure_location = 'Figures';
SAVE_MY_FIGURE(fig1_comps, sprintf('%s\\MachNumber_vs_Altitude.png', figure_location), 'small');

% Reynolds Number
figure(2);
fig2_comps.fig = gcf;
% Plot the data
fig2_comps.p1 = plot(h_list, Re_h, 'LineWidth', 1.25);

% Set properties
title('Reynolds Number variation with Altitude');
xlabel('Altitude h');
ylabel('Reynolds Number Re');

figure_location = 'Figures';
SAVE_MY_FIGURE(fig2_comps, sprintf('%s\\ReynoldsNumber_vs_Altitude.png', figure_location), 'small');


%==================================================
% PRINT RESULTS OF THIS SECTION

disp('MACH NUMBER AND REYNOLDS NUMBER VARIATION WITH ALTITUDE');
disp('--------------------');

fprintf('Minimum Mach No. (at h = %.2f) = %.4f\n', h_M_min, M_min);
fprintf('Maximum Mach No. (at h = %.2f) = %.4f\n', h_M_max, M_max);
fprintf('Minimum Reynolds No. (at h = %.2f) = (%.2f)%.4e\n', h_Re_min, Re_min, Re_min);
fprintf('Maximum Reynolds No. (at h = %.2f) = (%.2f)%.4e\n', h_Re_max, Re_max, Re_max);

fprintf('\n\n');


%========================================================
% CLOSE ALL FIGURES
close all;


%% VARIATION OF CL WITH WING LOADING

%==================================================
% OPERATING RANGE FOR MACH NUMBER AND REYNOLDS NUMBER

%==================================================
% CL RANGE BASED ON THE RANGE OF WING LOADING OBTAINED EARLIER

WL_min = WL_interval_common_min;
WL_max = WL_interval_common_max;


%==================================================
% SET REQUIRED DATA VALUES

MTOW = W0_SecondEstimate;

N_h = 100;
h_min = 0;
h_max = Max_altitude;
h_list = linspace(h_min, h_max, N_h);
% Get std_atm object
std_atm = standardAtmosphere();
% Get area list
N_S = 1000;
S_min = (MTOW * g) / WL_max;
S_max = (MTOW * g) / WL_min;
S_list = linspace(S_min, S_max, N_S);
WL_list = (MTOW * g) ./ S_list;


%==================================================
% FIND THE RANGE OF Cl REQUIRED FOR THE RANGE OF WING LOADING AT DIFFERENT ALTITUDE

CL_list = zeros(N_h, N_S);
for i = 1:N_h
    h = h_list(i);
    rho_h = std_atm.density(h / 1000);
    q = (1/2) * rho_h * V_cruise * V_cruise;
    CL_list(i, :) = WL_list / q;
end

% Plot variation of CL with Wing Loading at different altitudes
figure(1);
fig1_comps.fig = gcf;
% Plot the data
hold on
i_list = 0: (N_h / 5): N_h;
i_list(1) = 1;
for i = i_list
    fig1_comps.p1(i) = plot(WL_list, CL_list(i, :), 'DisplayName', sprintf('h = %.2f', h_list(i)), 'LineWidth', 1.25);
end

% Set properties
xlabel('$$Wing Loading \; \frac{W}{S}$$');
ylabel('$$Lift Coefficient \; C_{l}$$');
legend();
legendX = 0.75; legendY = 0.77; legendWidth = .1; legendHeight = .1;
fig1_comps.legendPosition = [legendX, legendY, legendWidth, legendHeight];

STANDARDIZE_FIGURE(fig1_comps);
figure_location = 'Figures';
SAVE_MY_FIGURE(fig1_comps, sprintf('%s\\CL_variation_with_WingLoading_at_different_altitudes.png', figure_location), 'small');

% Find Maximum range of CL possible. The following variables are for the maximum of the lower value of each interval of CL at an altitude.
CL_max_lower_idx = 1;
CL_max_lower = 0;
for i = 1:N_h
    if min(CL_list(i, :)) > CL_max_lower
        CL_max_lower = min(CL_list(i, :));
        CL_max_lower_idx = i;
    end
end

% Get Power Loading at that Altitude and the Interval of Wing Loading for 5% Variation of Power Loading
h_maxCL = h_list(CL_max_lower_idx);
rho_h_maxCL = std_atm.density(h_maxCL / 1000);
q = (1/2) * rho_h_maxCL * V_cruise * V_cruise;
CL_max_lower = WL_min / q;
CL_max_upper = WL_max / q;

CL_min_all = min(CL_list(:));
CL_max_all = max(CL_list(:));


%==================================================
% PRINT RESULTS OF THIS SECTION

disp('VARIATION OF CL WITH WING LOADING');
disp('--------------------');

fprintf('CL maximum interval lower value = %.2f\n', CL_max_lower);
fprintf('CL maximum interval upper value = %.2f\n', CL_max_upper);

fprintf('CL minimum across altitudes = %.2f\n', CL_min_all);
fprintf('CL maximum across altitudes = %.2f\n', CL_max_all);

fprintf('\n\n');


%========================================================
% CLOSE ALL FIGURES
close all;


%% SELECT PARTICULAR WING LOADING

%========================================================
% SELECT WING LOADING

% We select that wing loading that corresponds to the least Cl_cruise at max altitude
WingLoading = WL_min;
Wing_area = (W0_SecondEstimate * g) / WingLoading;

% Corresponding Cl
q = (1/2) * rho_ceil * V_cruise * V_cruise;
CL_cruise_ceil = WingLoading / q;

% Get Battery Weight Required
h = 0;
[~, ~, ~, ~, ~, ~, Battery_weight] = Power_Calculation_func(h, W0_SecondEstimate, Wing_area, AR, A_prop);

% Check that the predicted Battery_weight is lower than predicted by second weight estimate
if Battery_weight > W_Second_BatteryWeight
    disp('ALERT!!   BATTERY WEIGHT MORE THAN SECOND WEIGHT ESTIMATE   ALERT!!');
    return
end


%==================================================
% PRINT RESULTS OF THIS SECTION

disp('SELECT PARTICULAR WING LOADING');
disp('--------------------');

fprintf('WingLoading = %.2f\n', WingLoading);
fprintf('CL_cruise_ceil = %.2f\n', CL_cruise_ceil);

fprintf('Battery_weight = %.2f\n', Battery_weight);

fprintf('\n\n');


%========================================================
% CLOSE ALL FIGURES
close all;


%% SELECT PARTICULAR DISK LOADING

%========================================================
% SELECT DISK LOADING

% We select that wing loading that corresponds to the least Cl_cruise at max altitude
WingLoading = WL_min;
Wing_area = (W0_SecondEstimate * g) / WingLoading;

% Corresponding Cl
q = (1/2) * rho_ceil * V_cruise * V_cruise;
CL_cruise_ceil = WingLoading / q;

% Get Battery Weight Required
h = 0;
[~, ~, ~, ~, ~, ~, Battery_weight] = Power_Calculation_func(h, W0_SecondEstimate, Wing_area, AR, A_prop);

% Check that the predicted Battery_weight is lower than predicted by second weight estimate
if Battery_weight > W_Second_BatteryWeight
    disp('ALERT!!   BATTERY WEIGHT MORE THAN SECOND WEIGHT ESTIMATE   ALERT!!');
    return
end


%==================================================
% PRINT RESULTS OF THIS SECTION

disp('SELECT PARTICULAR WING LOADING');
disp('--------------------');

fprintf('WingLoading = %.2f\n', WingLoading);
fprintf('CL_cruise_ceil = %.2f\n', CL_cruise_ceil);

fprintf('\n\n');


%========================================================
% CLOSE ALL FIGURES
close all;


%% SELECT A PARTICULAR AIRFOIL

%==================================================
% SELECT A PARTICULAR AIRFOIL

% Based on the Cl requirements as given above we pick this airfoil - GOE 802
airfoil_name = 'GOE802';

% Airfoil characteristics
file_location = 'DataFiles/GOE802_Data_Airfoiltools.csv';
[Cl_alpha, Cl_0_alpha] = AirfoilData_func(file_location);


%==================================================
% PLOTS OF THE AIRFOIL CHARACTERISTICS


%==================================================
% PRINT RESULTS OF THIS SECTION

% Print airfoil characteristics and make plots
disp('AIRFOIL');
disp('--------------------');

fprintf('Airfoil = %s\n', airfoil_name);

fprintf('Cl_alpha = %.4f\n', Cl_alpha);
fprintf('Cl_0_alpha = %.4frad (%.4f deg)\n', Cl_0_alpha, Cl_0_alpha * (180 / pi));

fprintf('\n\n');


%========================================================
% CLOSE ALL FIGURES
close all;


%% WING DESIGN

%==================================================
% EFFECTS OF AR, SWEEP AND TAPER RATIO

% Aspect ratio
AR;
% Sweep
sweep = 0;
% Taper ratio
taper_ratio = 0.45;

M_avg = (M_min + M_max) / 2;
beta = sqrt(1 - M_avg^2);
eta = Cl_alpha / (2*pi);

% Lift curve slope
CL_alpha = (2 * pi * AR) ./ (2 + sqrt(4 + ((AR.^2 .* beta.^2) ./ eta.^2)*(1 + (tan(sweep).^2 ./ beta.^2))));
% For untwisted wing zero lift angle of attack is same for the wing and airfoil
CL_0_alpha = Cl_0_alpha;


%==================================================
% OTHER WING PARAMETERS

% Wing area
Wing_area;
Wing_span = sqrt(AR * Wing_area);
% Chord lengths
Root_chord = (2 * Wing_area) / (Wing_span * (1 + taper_ratio));
Tip_chord = taper_ratio * Root_chord;
Mean_chord = (2/3) * Root_chord * ((1 + taper_ratio + taper_ratio^2) / (1 + taper_ratio));
% Wing twist
Wing_twist = 0;

% Wing setting angle determination
std_atm = standardAtmosphere();
h_min = 0;
h_max = Max_altitude;
h_design = h_min * (65/100) + h_max * (35/100);
% h_design = 0;
% h_design = Max_altitude;
rho_design = std_atm.density(h_design / 1000);
CL_design = WingLoading / ((1/2) * rho_design * V_cruise * V_cruise);
Wing_setting = (CL_design / CL_alpha) + CL_0_alpha;

% Dihedral angle (in radians)
Wing_dihedral = 6 * (pi / 180);

% Wing vertical location
Wing_vertical_location = 'Mid-Wing Configuration';


%==================================================
% PRINT RESULTS OF THIS SECTION

% Print airfoil characteristics and make plots
disp('WING DESIGN');
disp('--------------------');

fprintf('AR = %.4f\n', AR);
fprintf('Sweep = %.4f rad (%.4f deg)\n', sweep, sweep * (180 / pi));
fprintf('Taper Ratio = %.4f\n', taper_ratio);
fprintf('CL_alpha = %.4f\n', CL_alpha);
fprintf('Wing Area = %.4f m sq.\n', Wing_area);
fprintf('Wing Span = %.4f m\n', Wing_span);
fprintf('Root Chord = %.4f m\n', Root_chord);
fprintf('Tip Chord = %.4f m\n', Tip_chord);
fprintf('Mean Chord = %.4f m\n', Mean_chord); %Shouldn't it be mean aerodynamic chord, because mean chord is just an arithmetic mean of chord lengths
fprintf('Wing Twist = %.4f rad (%.4f deg)\n', Wing_twist, Wing_twist * (180 / pi)); % Yes it is mean aerodynamic chord. The name would be too long so I used mean chord. Also it fits with the other names, same structure.
fprintf('Wing Setting = %.4f rad (%.4f deg)\n', Wing_setting, Wing_setting * (180 / pi));
fprintf('Dihedral = %.4frad (%.4f deg)\n', Wing_dihedral, Wing_dihedral * (180 / pi));
fprintf('Wing Vertical Location = %s\n', Wing_vertical_location);

fprintf('\n\n');


%========================================================
% CLOSE ALL FIGURES
close all;


%% TAIL DESIGN

%=======================================================
% HORIZONTAL TAIL
V_h_bar = 0.5;
d_f_aft = 40/1000;
l_opt = sqrt(2*Mean_chord*Wing_area*V_h_bar/(pi*d_f_aft));
S_h = V_h_bar*Mean_chord*Wing_area/l_opt;
C_m_af = -0.1136;
C_m0_wf = C_m_af*AR/(AR+2);
eta_h = 0.9;
delta_h_cg_ac = 0.2-0.25;
CL_h = (C_m0_wf+(CL_cruise_ceil*delta_h_cg_ac))/(eta_h*V_h_bar);
Cl_alpha_h = AirfoilData_func('DataFiles/xf-n0009sm-il-200000.csv'); %NACA0009
Sweep_h = sweep;
Dihedral_h = Wing_dihedral;
AR_h = 2*AR/3;
Taper_ratio_h = 0.7;
CL_alpha_h = Cl_alpha_h/(1+(Cl_alpha_h/(pi*AR_h)));
alpha_h = CL_h/CL_alpha_h;
epsilon0 = 2*CL_cruise_ceil/(pi*AR);
epsilon_alpha = 2*CL_alpha/(pi*AR);
i_h = alpha_h+epsilon0+(epsilon_alpha*Wing_setting);
b_h = sqrt(AR_h*S_h);
c_r_h = 2*S_h/(b_h*(1+Taper_ratio_h));
c_t_h = c_r_h*Taper_ratio_h;
c_h_bar = (2/3)*c_r_h*(1+Taper_ratio_h+(Taper_ratio_h^2))/(1+Taper_ratio_h);

%==================================================
% VERTICAL TAIL
V_v_bar = 0.04;
S_v = Wing_span*Wing_area*V_v_bar/l_opt;
AR_v = 1.5;
b_v = sqrt(AR_v*S_v);
Taper_ratio_v = 1;
c_r_v = 2*S_v/(b_v*(1+Taper_ratio_v));
c_t_v = c_r_v*Taper_ratio_v;
c_v_bar = (2/3)*c_r_v*(1+Taper_ratio_v+(Taper_ratio_v^2))/(1+Taper_ratio_v);
Sweep_v = atand((c_r_v-c_t_v)/b_v);

%==============================================
% CLOSE ALL FIGURES
close all;

%% DRAG ESTIMATION
%Parasitic Drag
%Skin Friction Coefficient
mu = 1.825e-5;
[~,a,~,rho] = atmosisa(0);
M = V_cruise/a;
l_f = 0.61047+l_opt;
w_f_fore = 0.136;
h_f_fore = 0.07405;

Re_f = Reynolds_number(rho,V_cruise,l_f,mu);
Re_W = Reynolds_number(rho,V_cruise,Mean_chord,mu);
Re_h = Reynolds_number(rho,V_cruise,c_h_bar,mu);
Re_v = Reynolds_number(rho,V_cruise,c_v_bar,mu);

Cf_f = Skin_Friction(Re_f,M);
Cf_W = Skin_Friction(Re_W,M);
Cf_h = Skin_Friction(Re_h,M);
Cf_v = Skin_Friction(Re_v,M);

%Form Factors
FF_W = Wing_Form_Factor(0.3,0.098,M);
FF_h = Wing_Form_Factor(0.309,0.09,M);
FF_v = Wing_Form_Factor(0.309,0.09,M);
d_f = sqrt(4*w_f_fore*h_f_fore/pi);
f = l_f/d_f;
FF_f = 0.9+(5/(f^1.5)+(f/400));

%Interference
Q_f = 1;
Q_W = 1;
Q_h = 1.04;
Q_v = 1.04;

%Wetted area
S_wet_W = Wing_Wetted_Area(Wing_Exposed_Area(Root_chord,Tip_chord,Wing_span,w_f_fore),0.098);
S_wet_h = Wing_Wetted_Area(Wing_Exposed_Area(c_r_h,c_t_h,b_h,d_f_aft),0.09);
S_wet_v = Wing_Wetted_Area(S_v,0.09);
S_wet_f = pi*d_f*l_f*((1-(2/f))^(2/3))*(1+(f^(-2)));

%Wetted Area Ratio
S_wet_S_ref_W = S_wet_W/Wing_area;
S_wet_S_ref_h = S_wet_h/Wing_area;
S_wet_S_ref_v = S_wet_v/Wing_area;
S_wet_S_ref_f = S_wet_f/Wing_area;

CD0_f = Cf_f*FF_f*Q_f*S_wet_S_ref_f;
CD0_W = Cf_W*FF_W*Q_W*S_wet_S_ref_W;
CD0_h = Cf_h*FF_h*Q_h*S_wet_S_ref_h;
CD0_v = Cf_v*FF_v*Q_v*S_wet_S_ref_v;
CD0_final = (CD0_f+CD0_W+CD0_h+CD0_v)*1.1;

%Induced Drag
e_final = (1.78*(1-(0.045*(AR^0.68))))-0.64;
K_final = 1/(pi*e_final*AR);
function Re = Reynolds_number(rho,V,l,mu)
    Re = rho*V*l/mu;
end

function Cf = Skin_Friction(Re,M)
    Cf = 0.455/(((log10(Re))^2.58)*((1+(0.144*(M^2)))^0.65));
end

function FF = Wing_Form_Factor(x_by_c,t_by_c,M)
    FF = (1+(0.6*t_by_c/x_by_c)+(100*(t_by_c^4)))*(1.34*(M^0.18));
end

function S_wet = Wing_Wetted_Area(S_exp,t_by_c_root)
    S_wet = 2*S_exp*(1+(0.25*t_by_c_root));
end

function S_exp = Wing_Exposed_Area(c_r,c_t,b,w_f)
    S_exp = ((2*c_t)+((b-w_f)*(c_r-c_t)/b))*(b-w_f)/2;
end