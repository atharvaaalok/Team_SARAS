clc; clear; close all;

%Parasitic Drag
%Skin Friction Coefficient
mu = 1.825e-5;
[~,a,~,rho] = atmosisa(0);
V = 20;
M = V/a;
l_f = 0.6+0.6252;
l_W = 0.2830;
l_h = 0.159;
l_v = 0.2164;
w_F = 0.126;
h_F = 0.06405;
S_ref = 0.5250;

Re_f = Reynolds_number(rho,V,l_f,mu);
Re_W = Reynolds_number(rho,V,l_W,mu);
Re_h = Reynolds_number(rho,V,l_h,mu);
Re_v = Reynolds_number(rho,V,l_v,mu);

Cf_f = Skin_Friction(Re_f,M);
Cf_W = Skin_Friction(Re_W,M);
Cf_h = Skin_Friction(Re_h,M);
Cf_v = Skin_Friction(Re_v,M);

%Form Factors
FF_W = Wing_Form_Factor(0.3,0.098,M);
FF_h = Wing_Form_Factor(0.309,0.09,M);
FF_v = Wing_Form_Factor(0.309,0.09,M);
f = l_f/sqrt(4*w_F*h_F/pi);
FF_f = 0.9+(5/(f^1.5)+(f/400));

%Interference
Q_f = 1;
Q_W = 1;
Q_h = 1.04;
Q_v = 1.04;

%Wetted area
S_wet_W = Wing_Wetted_Area(Wing_Exposed_Area(0.3725,0.1676,1.9443,w_F),0.098);
S_wet_h = Wing_Wetted_Area(Wing_Exposed_Area(0.1851,0.1296,0.7552,w_F),0.09);
S_wet_v = Wing_Wetted_Area(0.0653,0.09);
d_f = 2*(w_F+h_F)/pi;
lambda_f = l_f/d_f;
S_wet_f = pi*d_f*l_f*((1-(2/lambda_f))^(2/3))*(1+(lambda_f^(-2)));

%Wetted Area Ratio
S_wet_S_ref_W = S_wet_W/S_ref;
S_wet_S_ref_h = S_wet_h/S_ref;
S_wet_S_ref_v = S_wet_v/S_ref;
S_wet_S_ref_f = S_wet_f/S_ref;

CD0_f = Cf_f*FF_f*Q_f*S_wet_S_ref_f;
CD0_W = Cf_W*FF_W*Q_W*S_wet_S_ref_W;
CD0_h = Cf_h*FF_h*Q_h*S_wet_S_ref_h;
CD0_v = Cf_v*FF_v*Q_v*S_wet_S_ref_v;
CD0 = CD0_f+CD0_W+CD0_h+CD0_v;

%Induced Drag
AR = 7.2;
e = (1.78*(1-(0.045*(AR^0.68))))-0.64;
K = 1/(pi*e*AR);
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
    
    