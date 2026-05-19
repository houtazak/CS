%---------------------------------------------------
%  NAME:      main_nonlinear_CS_model_electrical_machine_GUO2020.m
%  WHAT:      Nonlinear current sheet of an electrical machine (https://doi.org/10.1109/TMAG.2019.2950614)
%  REQUIRED:  CSmodel toolbox 20200321
%  AUTHOR:    20200321, J. Guo, L. Quéval (loic.queval@gmail.com), F. Trillaud, B. Roucaries
%  MODIFICATIONS: 2026, Zakaria HOUTA
%   - Nonlinear current sheet model of an electrical machine with computational improvements,
%     and calculation of stator electromagnetic power using stator flux linkages and analytical torque.
%   - Loop vectorisation, improved nonlinear convergence using an Anderson accelerator,
%     and RMS-based evaluation of iron permeability dependence on magnetic flux density.
%   - General numerical optimisations to speed up execution. 
%  COPYRIGHT: 2020, Loïc Quéval, BSD License (http://opensource.org/licenses/BSD-3-Clause)
%----------------------------------------------------
clear all, close all, clc

addpath('..\CS Core') %add CSmodel toolbox functions to the path

%% Constants
mu0 = 4*pi*1e-7;

%% Variables
isnonlin = 1; %0 for linear case, 1 for nonlinear case
h_min = 1; %min harmonic number to consider (must be odd)
h_max = 13; %max harmonic number to consider



%% Define machine   
% Machine geometry
P = 6;             %number pole pair []
r1 = 1.320;        %rotor inner radius [m]
r2 = 1.470;        %rotor outer radius [m]
r3 = 1.546;        %radius of field coil [m]
rTe = 1.619;       %radius for estimation of Te [m]
r4 = 1.683;        %radius of amature coil [m]
r5 = 1.750;        %stator inner radius [m]
r6 = 2.000;        %Stator outer radius [m]
L_eff = 1.540;     %Effective length of the machine [m]

% Rotor winding
Nf = 100;           %field coil turns []
hf = 0.057;         %field coil height [m]
theta_1f = 0.163;   %field coil width angle [elec. rad]
theta_2f = 2.703;   %field coil aperture angle [elec. rad]
wf = 0.042;         %field coil width [m]
alfa_m_deg_ini = -15;	%initial rotor angle [mec. deg]  -15
alfa_m_rad_ini = alfa_m_deg_ini*pi/180;   %rotor angle [mec. rad]
i_f = 5.03e3;      %field winding current [A]

% Stator winding
Na= 120;            %amature coil turns []
ha = 0.041;         %armature coil height [m]
theta_1a = 0.692;   %armature coil width angle [elec. rad]
theta_2a = 0.664;   %armature coil aperture angle [elec. rad]
wa = 0.194;         %armature coil width [m]
Im = 1e3;     %amplitude of the stator currents

% Nominal speed of this synchronous machine : ws=(60*f/P) [tr/min] = (2*pi*f/P) [rad/s]
freq= 50; % frequancy [Hz]
ws = 2*pi*freq/P; % nominal speed [rad/s]

% iron
fct_mur_B([]); %Define iron (and plot)

% build machine geom
r_l = [r1,r2,r3,r4,r5,r6]; %vector of current sheet radius
is_currentsheet = [0,0,1,1,0,0]; %0 for nothing, 1 for current sheet
mur_l = [1,1200,1,1,1,1200,1]; %vector of domain (initial) relative permeability (air, rotor iron, airgap, airgap, airgap, stator iron, air)
is_nonlinear = [0,1,0,0,0,1,0]*isnonlin; %0 for linear, 1 for nonlinear domain


% plot machine geometry
% fig1 = figure();
% box on, grid on, hold on
% fct_plot_geom_generic(r_l,is_currentsheet)
% xlabel('x [m]'), ylabel('y [m]')


% Disable warnings regarding the singularity of the matrix
w = warning('off', 'MATLAB:nearlySingularMatrix')



%% Calculation of electromagnetic torque, magnetic flux and electromagnetic power 
Nt = 300; %Number of simulation steps
Tf = 2*pi/ws; %Final moment to define (Here is a full rotation of the machine)
tsim = linspace(0, Tf, Nt); %Discretisation of time points

phivalue = 0; %Phase shift of the stator currents [rad]


% Storage of the instantaneous electromagnetic torque (calculated using the Maxwell stress tensor)
Torque_Maxwell = zeros(size(tsim));

% Initialisation for storing the instantaneous magnetic fluxes associated with the A, B and C phases of the stator
lambdaA = zeros(size(tsim));
lambdaB = zeros(size(tsim));
lambdaC = zeros(size(tsim));

% Initialisation for storing the instantaneous stator currents
ia_vec = zeros(size(tsim));    %Current in phase A
ib_vec = zeros(size(tsim));    %Current in phase B
ic_vec = zeros(size(tsim));    %Current in phase C

% Initialize current sheets
K_lh_s = zeros(length(r_l),h_max);
K_lh_c = zeros(length(r_l),h_max);

% Initialisation for storing the number of iterations of the non-linear solver
nbre_iter_NL_history = zeros(size(tsim));

tic;
for k = 1:length(tsim) 
    k
    t_loop = tsim(k); %time [s]
    theta_m = ws*t_loop + alfa_m_rad_ini; %rotor angle [mec. rad]
    theta_e = P*theta_m; %rotor angle [elec. rad]

    % Balanced stator currents (as a function of the electrical angle)
    ia = Im * sin(theta_e + phivalue); ia_vec(k) = ia;
    ib = Im * sin(theta_e + phivalue - 2*pi/3); ib_vec(k) = ib;
    ic = Im * sin(theta_e + phivalue - 4*pi/3); ic_vec(k) = ic;

    % CS at the rotor : K3  
    [K_lh_c(3,:), K_lh_s(3,:)] = fct_get_K3(h_min,h_max,Nf,i_f,theta_1f,theta_2f,theta_m,P,wf);

    % CS at the stator : K4
    [K_lh_c(4,:), K_lh_s(4,:)] = fct_get_K4(Na, ia, ib, ic, h_min, h_max, theta_1a, theta_2a, wa);
    
    % Resolution
    [a_lh,b_lh,c_lh,d_lh,mur_l, nbre_iter_NL] = fct_get_XX_nonlin(h_min,h_max,P,r_l,mur_l,mu0,K_lh_s,K_lh_c,is_nonlinear);
    nbre_iter_NL_history(k) = nbre_iter_NL; %Storing the number of iterations of the non-linear solver

    % Instantaneous electromagnetic torque 
    Torque_Maxwell(k) = Get_AnalyticalTorque(L_eff, P, mu0, a_lh, b_lh, c_lh, d_lh, 4, h_min, h_max);
    
    % Instantaneous stator fluxes
    lambdaA(k) = Get_flux_CS_in_stator_phase( a_lh, b_lh, c_lh, d_lh, h_min, h_max, P, r4, ha, theta_1a, theta_2a, Na, L_eff, 0); %At Phase A
    lambdaB(k) = Get_flux_CS_in_stator_phase( a_lh, b_lh, c_lh, d_lh, h_min, h_max, P, r4, ha, theta_1a, theta_2a, Na, L_eff, 2*pi/3); %At Phase B
    lambdaC(k) = Get_flux_CS_in_stator_phase( a_lh, b_lh, c_lh, d_lh, h_min, h_max, P, r4, ha, theta_1a, theta_2a, Na, L_eff, 4*pi/3); %At Phase C
 
end
Tfinal=toc


% Average electromagnetic torque [N.m]
C_mean = trapz(tsim, Torque_Maxwell) / (tsim(end)-tsim(1));
% Instantaneous electromagnetic power calculated from the electromagnetic torque [W]
P_em_torque_inst = Torque_Maxwell * ws;
% Electromagnetic power calculated from the electromagnetic torque [W]   
P_em_torque = C_mean * ws; 

% Calculation of induced phase EMFs from the magnetic fluxes in phases A, B and C [V]
Ea = - gradient(lambdaA, tsim);
Eb = - gradient(lambdaB, tsim);
Ec = - gradient(lambdaC, tsim);
% Instantaneous electromagnetic power calculated from the electromagnetic fluxes [W]
P_em_flux_inst = Ea .* ia_vec + Eb .* ib_vec + Ec .* ic_vec;
% Electromagnetic power calculated from the electromagnetic fluxes [W]
P_em_flux = trapz(tsim, P_em_flux_inst) / (tsim(end)-tsim(1));
    
% Displaying results
fprintf('\n=== ENERGY AUDIT ===\n');
fprintf('Average electromagnetic torque: %.2f N.m\n', C_mean);
fprintf('Electromagnetic power calculated from the electromagnetic torque: %.2f MW\n', P_em_torque/1e6);
fprintf('Electromagnetic power calculated from the electromagnetic fluxes: %.2f MW\n', P_em_flux/1e6);
fprintf('Relative difference in the electromagnetic power calculated using the two methods: %.4f %%\n', abs(P_em_torque - P_em_flux)/abs(P_em_torque)*100);


% Visualisation
figure('Position',[100 100 1000 900]);
%%% =========================================================
% Instantaneous stator fluxes
%%% =========================================================
subplot(3,1,1);
plot(tsim, lambdaA, 'r', tsim, lambdaB, 'g', tsim, lambdaC, 'b', 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('Fluxes [Wb]');
title('Instantaneous stator fluxes');
legend('\lambda_a','\lambda_b','\lambda_c','Location','best');
grid on;

%%% =========================================================
% Induced phase EMFs
%%% =========================================================
subplot(3,1,2);
plot(tsim, Ea, 'r', tsim, Eb, 'g', tsim, Ec, 'b', 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('EMFs [V]');
title('Induced phase EMFs');
legend('e_a','e_b','e_c','Location','best');
grid on;

%%% =========================================================
% Electromagnetic power [MW]   
%%% =========================================================
subplot(3,1,3);
plot(tsim, P_em_flux_inst/1e6, 'k', 'LineWidth', 1.5); hold on;
plot(tsim, P_em_torque_inst/1e6,  'r', 'LineWidth', 1.5);

% Moyennes
yline(P_em_flux/1e6, 'k--', 'LineWidth', 1.5);
yline(P_em_torque/1e6,  'r--', 'LineWidth', 1.5);

xlabel('Time [s]');
ylabel('Electromagnetic power [MW]');
title('Instantaneous and average electromagnetic power');

legend({'P_{em, flux}(t)', ...
        'P_{em, torque}(t)', ...
        sprintf('<P>_{em, flux} = %.2f MW', P_em_flux/1e6), ...
        sprintf('<P>_{em, torque} = %.2f MW', P_em_torque/1e6)}, ...
        'Location','best');

grid on;







  