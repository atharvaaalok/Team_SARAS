function [Cl_alpha, Cl_0_alpha] = AirfoilData_func(file_location)
    
    airfoil_data = readtable(file_location, 'NumHeaderLines', 1);
    airfoil_data = airfoil_data{:, :};

    alpha = airfoil_data(:, 1);
    Cl = airfoil_data(:, 2);
    Cd = airfoil_data(:, 3);
    Cdp = airfoil_data(:, 4);
    Cm = airfoil_data(:, 5);
    
    figure(1);
    plot(alpha, Cl);
    grid on
    
    tolerance = 0.06;
    alpha_min = 0;
    alpha_min_idx = find(abs(alpha-alpha_min) < tolerance);
    alpha_min_idx = alpha_min_idx(1);
    alpha_max = 5;
    alpha_max_idx = find(abs(alpha-alpha_max) < tolerance);
    alpha_max_idx = alpha_max_idx(1);
    
    Cl_alpha_min = Cl(alpha_min_idx);
    Cl_alpha_max = Cl(alpha_max_idx);
    
    Cl_alpha = (Cl_alpha_max - Cl_alpha_min) / (alpha_max - alpha_min);
    Cl_alpha = Cl_alpha * (180 / pi);
    
    CL_vs_alpha_poly_fit = polyfit(alpha, Cl, 5);
    alpha_finer = linspace(min(alpha), max(alpha), 10000);
    Cl_finer = polyval(CL_vs_alpha_poly_fit, alpha_finer);
    tolerance = 0.01;
    Cl_0_idx = find(abs(Cl_finer - 0) < tolerance);
    Cl_0_idx = Cl_0_idx(1);
    
    Cl_0_alpha = alpha_finer(Cl_0_idx) * (pi / 180);
    

end