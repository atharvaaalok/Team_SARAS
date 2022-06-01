function [P_cruise, P_takeoff_1, P_descent_1, P_takeoff_2, P_descent_2, Total_energy_required, Battery_weight] = Power_Calculation_func(h, MTOW, Wing_area, AR, A_prop)
    
    global g V_cruise oswald_efficiency VTOL_motor_count V_takeoff V_descent W_payload_dropped Cruise_time TakeOff_time Descent_time prop_efficiency motor_efficiency ESC_efficiency figure_of_merit_hoverpower battery_SOH battery_SOC battery_discharge_efficiency Battery_specific_energy battery_PIF battery_redundancy_ratio 
    
    std_atm = standardAtmosphere();
    rho_h = std_atm.density(h / 1000);

    %==================================================
    % AERODYNAMIC QUANTITIES

    % Design lift coefficient estimation
    % In cruise W = L = (1/2) * rho * V^2 * S * Cl

    CL_cruise = (g * MTOW) ./ (0.5 * rho_h * V_cruise * V_cruise * Wing_area);
    CD0 = CD0_func(Wing_area);
    CD_cruise = CD0 + (1 / (pi * AR * oswald_efficiency)) * (CL_cruise.^2);


    %==================================================
    % POWER REQUIRED CALCULATIONS

    % Cruise
    T_cruise = 0.5 * rho_h * V_cruise * V_cruise .* Wing_area .* CD_cruise;
    P_cruise = V_cruise * T_cruise;

    % Vertical take off 1
    K_T = 1.2;
    T_takeoff_1 = K_T * (g * MTOW / VTOL_motor_count);
    P_takeoff_1 = VTOL_motor_count * (V_takeoff * T_takeoff_1 / 2) * (1 + sqrt(1 + ((2 * T_takeoff_1) ./ (rho_h * (V_takeoff^2) * A_prop))));

    % Vertical descent 1
    T_hover_1 = g * MTOW / VTOL_motor_count;
    V_propeller_hover = sqrt(T_hover_1 ./ (2 * rho_h * A_prop));
    x = -1 * (V_descent ./ V_propeller_hover);
    K = 1.2;
    V_propeller_induced = (K - 1.125*x - 1.372*(x.^2) - 1.718*(x.^3) - 0.655*(x.^4)) .* V_propeller_hover;
    P_descent_1 = VTOL_motor_count * K * (g * MTOW / VTOL_motor_count) * (V_propeller_induced - V_descent);

    % Vertical take off 2
    K_T = 1.2;
    T_takeoff_2 = K_T * (g * (MTOW - W_payload_dropped) / VTOL_motor_count);
    P_takeoff_2 = VTOL_motor_count * (V_takeoff * T_takeoff_2 / 2) * (1 + sqrt(1 + ((2 * T_takeoff_2) ./ (rho_h * (V_takeoff^2) * A_prop))));

    % Vertical descent 2
    T_hover_2 = g * (MTOW - W_payload_dropped) / VTOL_motor_count;
    V_propeller_hover = sqrt(T_hover_2 ./ (2 * rho_h * A_prop));
    x = -1 * (V_descent ./ V_propeller_hover);
    K = 1.2;
    V_propeller_induced = (K - 1.125*x - 1.372*(x.^2) - 1.718*(x.^3) - 0.655*(x.^4)) .* V_propeller_hover;
    P_descent_2 = VTOL_motor_count * K * (g * (MTOW - W_payload_dropped) / VTOL_motor_count) * (V_propeller_induced - V_descent);

    
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


end