function lambda = Get_flux_CS_in_stator_phase( a_lh, b_lh, c_lh, d_lh, h_min, h_max, P, r_coil, h_coil, theta_1a, theta_2a, Nturns, L_eff, theta_shift)
% %---------------------------------------------------
%  NAME:      Get_flux_CS_in_stator_phase.m
%  WHAT:      Calculation of the total flux lambda in a stator phase.
%  AUTHORS:  Zakaria HOUTA
%
%  USE:
%    [lambda] =
%    Get_flux_CS_in_stator_phase(a_lh, b_lh, c_lh, d_lh, h_min, h_max, P, r_coil, h_coil, theta_1a, theta_2a, Nturns, L_eff, theta_shift);
%
%  INPUTS:
%    a_lh,b_lh,c_lh,d_lh   = coefficients of the x-matrix (Eq.6)
%    h_min,h_max           = min and max harmonic number to be calculated
%    P          = nb of pairs of pole
%    r_coil     = radius of amature coil
%    h_coil     = armature coil height
%    theta_1a   = armature coil width angle [elec. rad]
%    theta_2a   = armature coil aperture angle [elec. rad]
%    Nturns     = amature coil turns
%    L_eff      = effective length of the machine
%    theta_shift = phase shift angle
%
%  OUTPUTS:
%    lambda         = The total flux lambda in a stator phase.
%----------------------------------------------------
    
    % index of the stator windings sub-domain
    l_stator_core = 5;
    
    % Extract the coefficients from stator windings sub-domain
    a_h = a_lh(l_stator_core, h_min:h_max);
    b_h = b_lh(l_stator_core, h_min:h_max);
    c_h = c_lh(l_stator_core, h_min:h_max);
    d_h = d_lh(l_stator_core, h_min:h_max);
    
    % Radii limiting the winding
    r_int = r_coil - h_coil/2;
    r_ext = r_coil + h_coil/2;
    
    % Mechanical angles
    theta1_plus = (theta_2a/2 + theta_shift) / P;
    theta2_plus = (theta_2a/2 + theta_1a + theta_shift) / P;
    theta1_minus = (-theta_2a/2 - theta_1a + theta_shift) / P;
    theta2_minus = (-theta_2a/2 + theta_shift) / P;
    
    % Calculation of integrals and areas
    [I_plus, A_plus] = integrate_Az_over_section(a_h, b_h, c_h, d_h, r_int, r_ext, theta1_plus, theta2_plus, P, h_min, h_max);
    [I_minus, A_minus] = integrate_Az_over_section(a_h, b_h, c_h, d_h, r_int, r_ext, theta1_minus, theta2_minus, P, h_min, h_max);
    
    % Calculation of averages
    mean_Az_plus = I_plus / A_plus;
    mean_Az_minus = I_minus / A_minus;
    
    % Total Flux 
    lambda = -P * L_eff * Nturns * (mean_Az_plus - mean_Az_minus);
end