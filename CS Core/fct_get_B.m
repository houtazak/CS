function [RHO,THETA,BR,BTH,NORMB] = fct_get_B(a_lh,b_lh,c_lh,d_lh,h_min,h_max,P,r_l,RHO,THETA)
%---------------------------------------------------
%  NAME:      fct_get_B.m
%  WHAT:      Calculates B at field points.
%  REQUIRED:  CSmodel toolbox 20200321
%  AUTHORS:   20200321, J. Guo, L. Quéval (loic.queval@gmail.com)
%  MODIFICATIONS:
%  2026, Zakaria Houta
%  - Vectorisation instead of using lots of loops.
%  COPYRIGHT: 2020, Loïc Quéval, BSD License (http://opensource.org/licenses/BSD-3-Clause)
%
%  USE:
%    [RHO,THETA,BR,BTH,NORMB] = fct_get_B(a_lh,b_lh,c_lh,d_lh,h_min,h_max,P,r_l,RHO,THETA);
%
%  INPUTS:
%    a_lh,b_lh,c_lh,d_lh   = coefficients of the x-matrix (Eq.6)
%    h_min,h_max           = min and max harmonic number to be calculated
%    P          = nb of pairs of pole
%    r_l        = geometry of the machine
%    RHO        = Field points rho-coordinate vector or matrix
%    THETA      = Field points theta-coordinate vector or matrix
%
%  OUTPUTS:
%    RHO        = Field points rho-coordinate vector or matrix
%    THETA      = Field points theta-coordinate vector or matrix
%    BR         = Field points B rho-component vector or matrix
%    BTH        = Field points B theta-component vector or matrix
%    NORMB      = Field points B norm vector or matrix
%---------------------------------------------------
    % Get dimensions
    [M, Q] = size(RHO);
    
    % Convert RHO and THETA into column vectors for processing
    rho_vec = RHO(:);
    theta_vec = THETA(:);
    
    % Initialize BR et BTH
    BR = zeros(M,Q);
    BTH = zeros(M,Q);
    
    % Determine the sectors for all points at once
    sectors = zeros(size(rho_vec));
    for i = 1:length(r_l)
        sectors = sectors + (rho_vec > r_l(i));
    end
    sectors = sectors + 1;
    
    % For each harmonic
    for h_idx = h_min:1:h_max    
        hP = h_idx * P;
    
        % Calculate the powers of rho for all points
        rho_pow_pos = rho_vec.^(hP - 1);
        rho_pow_neg = rho_vec.^(-hP - 1);
    
        % Calculate trigonometric functions
        hP_theta = hP * theta_vec;
        cos_hP_theta = cos(hP_theta);
        sin_hP_theta = sin(hP_theta);
    
        % Extract the coefficients for this harmonic and the corresponding sectors
        % Note: a_lh, b_lh, c_lh, d_lh are of size [number_of_sectors × number_of_harmonics]
        a_vals = a_lh(sectors, h_idx);
        b_vals = b_lh(sectors, h_idx);
        c_vals = c_lh(sectors, h_idx);
        d_vals = d_lh(sectors, h_idx);
    
        % Calculate the terms in the expression for magnetic flux density for BR and BTH
        term1_br = (a_vals .* rho_pow_pos + b_vals .* rho_pow_neg) * hP;
        term2_br = (c_vals .* rho_pow_pos + d_vals .* rho_pow_neg) * hP;
    
        term1_bth = (-a_vals .* rho_pow_pos + b_vals .* rho_pow_neg) * hP;
        term2_bth = (-c_vals .* rho_pow_pos + d_vals .* rho_pow_neg) * hP;
    
        % Calculate the contributions for this harmonic
        br_contrib = term1_br .* cos_hP_theta - term2_br .* sin_hP_theta;
        bth_contrib = term1_bth .* sin_hP_theta + term2_bth .* cos_hP_theta;
    
        % Add the contributions (restoring them to their original form)
        BR = BR + reshape(br_contrib, M, Q);
        BTH = BTH + reshape(bth_contrib, M, Q);
    end
    
    NORMB = sqrt(BR.^2+BTH.^2);






