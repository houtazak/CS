function [I, Area] = integrate_Az_over_section(a_h, b_h, c_h, d_h, r_in, r_out, theta1, theta2, P, h_min, h_max)
% %--------------------------------------------------------------------------------------
%  NAME:      integrate_Az_over_section.m
%  WHAT:      Calculate the integral of A_z over the cross-section, as well as its area.
%  AUTHORS:  Zakaria HOUTA
%
%  USE:
%    [RHO,THETA,AZ] =
%    integrate_Az_over_section(a_lh,b_lh,c_lh,d_lh,h_min,h_max,P,r_l,RHO,THETA);
%
%  INPUTS:
%    a_h, b_h, c_h, d_h   = coefficients of the x-matrix for the relevant part only
%    r_in, r_out          = internal and external radii that limit the winding
%    theta1, theta2       = the angles that limit the winding
%    P                    = nb of pairs of pole
%    h_min,h_max          = min and max harmonic number to be calculated
%
%  OUTPUTS:
%    I        = the integral of A_z over the cross-section (=∫∫ A_z r dr dθ)
%    Area         = the winding area ( = ∫∫ r dr dθ )
%----------------------------------------------------------------------------------------
    
    h = (h_min:h_max)';
    hP = h * P;
    
    % Radial integrals
    I_r_plus = (r_out.^(hP+2) - r_in.^(hP+2)) ./ (hP + 2);
    I_r_minus = zeros(size(hP));
    mask = (abs(hP - 2) < 1e-12);
    if any(mask)
        I_r_minus(mask) = log(r_out) - log(r_in);
    end
    I_r_minus(~mask) = (r_out.^(-hP(~mask)+2) - r_in.^(-hP(~mask)+2)) ./ (-hP(~mask) + 2);
    
    % Angular integrals
    I_theta_sin = (cos(hP*theta1) - cos(hP*theta2)) ./ hP;
    I_theta_cos = (sin(hP*theta2) - sin(hP*theta1)) ./ hP;
    
    % Formatting of Fourier coefficients (row vectors)
    if size(a_h,1) > 1 && size(a_h,2) == 1
        a_h = a_h'; b_h = b_h'; c_h = c_h'; d_h = d_h';
    end
    
    % Sum over harmonics
    I = sum( (a_h .* I_r_plus' + b_h .* I_r_minus') .* I_theta_sin' + ...
             (c_h .* I_r_plus' + d_h .* I_r_minus') .* I_theta_cos' );
    
    % Area of the cross-section (surface element r dr dθ)
    Area = (theta2 - theta1) * (r_out^2 - r_in^2) / 2;
end