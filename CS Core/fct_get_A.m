function [RHO,THETA,AZ] = fct_get_A(a_lh,b_lh,c_lh,d_lh,h_min,h_max,P,r_l,RHO,THETA)
% %---------------------------------------------------
%  NAME:      fct_get_A.m
%  WHAT:      Calculates AZ at field points.
%  AUTHORS:   Zakaria HOUTA
%
%  USE:
%    [RHO,THETA,AZ] =
%    fct_get_A(a_lh,b_lh,c_lh,d_lh,h_min,h_max,P,r_l,RHO,THETA);
%
%  INPUTS:
%    a_lh,b_lh,c_lh,d_lh   = coefficients of the x-matrix (Eq.6)
%    h_min,h_max           = min and max harmonic number to be calculated
%    P          = nb of pairs of pole
%    r_l        = geometry of the machine
%    RHO        = Field points rho-coordinate vector or matrix
%    TH         = Field points theta-coordinate vector or matrix
%
%  OUTPUTS:
%    RHO        = Field points rho-coordinate vector or matrix
%    THETA      = Field points theta-coordinate vector or matrix
%    AZ         = Field points Az vector or matrix
%----------------------------------------------------
    AZ = zeros(size(RHO));
    secteurs = discretize(RHO, [-Inf, r_l(:)', Inf]);
    orig_size = size(RHO);

    % Pre-calculation of recurrence factors
    rho_pow_P = RHO .^ P;          
    rho_neg_pow_P = RHO .^ (-P);   

    % Initialisation for the first h
    p = h_min * P;
    rho_pos = RHO .^ p;
    rho_neg = RHO .^ (-p);

    for h = h_min:h_max
        sin_theta = sin(p * THETA);
        cos_theta = cos(p * THETA);

        a = reshape(a_lh(secteurs, h), orig_size);
        b = reshape(b_lh(secteurs, h), orig_size);
        c = reshape(c_lh(secteurs, h), orig_size);
        d = reshape(d_lh(secteurs, h), orig_size);

        AZ = AZ + (a .* rho_pos + b .* rho_neg) .* sin_theta ...
                + (c .* rho_pos + d .* rho_neg) .* cos_theta;

        % Update for the next harmonic (unless it is the last one)
        if h < h_max
            p = p + P;
            rho_pos = rho_pos .* rho_pow_P;
            rho_neg = rho_neg .* rho_neg_pow_P;
        end
    end
end

