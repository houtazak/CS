function T = Get_AnalyticalTorque(L_eff, P, mu0, a_lh, b_lh, c_lh, d_lh, l_airgap, h_min, h_max)
% %---------------------------------------------------
%  NAME:      Get_AnalyticalTorque.m
%  WHAT:      Calculates AZ at field points.
%  AUTHORS:  Zakaria HOUTA
%
%  USE:
%    T =
%    Get_AnalyticalTorque(L_eff, P, mu0, a_lh, b_lh, c_lh, d_lh, l_airgap, h_min, h_max);
%
%  INPUTS:
%    L_eff      = effective length of the machine
%    P          = nb of pairs of pole
%    mu0        = magnetic permeability of air
%    a_lh,b_lh,c_lh,d_lh   = coefficients of the x-matrix (Eq.6)
%    l_airgap              = the index of the domain corresponding to the machine's air gap
%    h_min,h_max           = min and max harmonic number to be calculated
%
%  OUTPUTS:
%    T         = Analytical electromagnetic torque
%----------------------------------------------------
    
    h_vec = h_min:h_max;
    % Extract the Fourier coefficients for the air gap subdomain 
    a_h = a_lh(l_airgap, h_vec);
    b_h = b_lh(l_airgap, h_vec);
    c_h = c_lh(l_airgap, h_vec);
    d_h = d_lh(l_airgap, h_vec);
    
    % Calculation of the sum
    sum_h = sum( (h_vec * P).^2 .* (a_h .* d_h - b_h .* c_h) );
    
    % Analytical electromagnetic torque
    T = (L_eff / mu0) * 2 * pi * sum_h;
end