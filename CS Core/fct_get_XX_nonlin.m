function [a_lh,b_lh,c_lh,d_lh,mur_l_new,nbre_iter_NL_history] = fct_get_XX_nonlin(h_min,h_max,P,r_l,mur_l,mu0,K_lh_s,K_lh_c,isnonlin)
%---------------------------------------------------
%  NAME:      fct_get_XX_nonlin.m
%  WHAT:      iterative procedure for the nonlinear case
%  AUTHORS:   20200321, J. Guo, L. Quéval (loic.queval@gmail.com), F. Trillaud, B. Roucaries
%  2026, Zakaria Houta
%   - Nonlinear current sheet model with improved computational efficiency
%     and better physical consistency in the nonlinear permeability update.
%   - Introduction of a Gauss-Legendre quadrature in the radial direction
%     to replace point evaluation of B. This enables a volumetric RMS
%     estimation of B^2 and defines an effective field Beff.
%   - Replacement of max(B) criterion by RMS-based averaging.
%   - Anderson acceleration is used to improve nonlinear convergence.
%   - Loop vectorisation and numerical optimisations for faster harmonic evaluation.
%  COPYRIGHT: 2020, Loïc Quéval, BSD License (http://opensource.org/licenses/BSD-3-Clause)
%
%  USE:
%    [a_lh,b_lh,c_lh,d_lh,mur_l_new,nbre_iter_NL_history] = fct_get_XX_nonlin(h_min,h_max,P,r_l,mur_l,mu0,K_lh_s,K_lh_c,isnonlin)
%
%  INPUTS:
%    h_min,h_max         = min and max harmonic number to be calculated
%    P                   = nb of pairs of pole
%    r_l                 = vector of current sheet radius
%    mur_l               = vector of (initial) domain relative permeability
%    mu0                 = permeability of vacuum
%    K_lh_s,K_lh_c       = Fourier coefficients of the current sheets
%    isnonlin            = 0 for linear, 1 for nonlinear case
%        
%  OUTPUTS:
%    a_lh,b_lh,c_lh,d_lh  = coefficients of the x-matrix (Eq.6)
%    mur_l_new            = updated material of the machine
%    nbre_iter_NL_history = Number of iterations taken by the non-linear solver
%----------------------------------------------------

% Paramètres 
Number_iter = 300;              % number of iterations
mmax = 2;                       % Anderson memory
n_radial = 5;                   % number of radial points radial averaging via Gauss-Legendre quadrature

N = length(mur_l);

% Initialize loop
if sum(isnonlin) > 0
    loop_counter = 1;
else
    loop_counter = Number_iter - 1;
end

mur_l_old = mur_l.^2;
mur_l_new = mur_l;
e_mur = 1e-6;

% History for Anderson acceleration
Xhist = [];
Fhist = [];

while ( max(abs(mur_l_old - mur_l_new)) >= e_mur && loop_counter < Number_iter )

    xk = mur_l_new;
    mu_l = xk * mu0;

    % Linear solution with current mu
    [a_lh,b_lh,c_lh,d_lh] = fct_get_XX_lin(h_min,h_max,P,r_l,mu_l,K_lh_s,K_lh_c);

    mur_l_old = xk;

    % --------------------------------------------------------------
    % Fixed-point candidate with radial averaging (Gauss-Legendre)
    % --------------------------------------------------------------
    mur_candidate = xk;

    for l = 1:N
        if isnonlin(l) > 0
            r_min = r_l(l-1);
            r_max = r_l(l);

            % Get Gauss-Legendre points and weights on [r_min, r_max]
            [r_quad, w_quad] = gauss_legendre(n_radial, r_min, r_max);

            sum_B2_weighted = 0;
            total_weight = 0;

            for ir = 1:n_radial
                r = r_quad(ir);
                % Angular grid (full 0.5° resolution)
                [THETA, RHO] = ndgrid((0:0.5:360)*pi/180, r);% 0.5
                [~,~,~,~,NORMXXX_j] = fct_get_B(a_lh,b_lh,c_lh,d_lh,...
                                                  h_min,h_max,P,r_l,RHO,THETA);

                % Volume element: r dr dθ
                % Gauss weight already accounts for dr, so we multiply by r
                weight = w_quad(ir) * r;
                sum_B2_weighted = sum_B2_weighted + weight * sum(NORMXXX_j(:).^2);
                total_weight = total_weight + weight * numel(NORMXXX_j);
            end

            Beff = sqrt(sum_B2_weighted / total_weight); % An RMS is used to determine the dependence of permeability on magnetic flux density
            mur_candidate(l) = fct_mur_B(Beff);
        end
    end

    % --------------------------------------------------------------
    % Anderson acceleration
    % --------------------------------------------------------------
    rk = (mur_candidate - xk).';

    Xhist(:, end+1) = xk;
    Fhist(:, end+1) = mur_candidate;
    if size(Xhist, 2) > mmax + 1
        Xhist(:, 1) = [];
        Fhist(:, 1) = [];
    end

    mur_next = mur_candidate;

    if size(Xhist, 2) >= 2
        m = min(mmax, size(Xhist, 2) - 1);
        idx0 = size(Xhist, 2) - m;

        X = Xhist(:, idx0:end);
        F = Fhist(:, idx0:end);

        G = F - X;
        dG = diff(G, 1, 2);
        dF = diff(F, 1, 2);

        if ~isempty(dG) && all(isfinite(dG(:)))
            gamma = pinv(dG) * rk;
            mur_accel = mur_candidate - (dF * gamma).';
            mur_accel = max(1, mur_accel);
            if all(isfinite(mur_accel)) && norm(mur_accel - xk, inf) <= norm(rk, inf)
                mur_next = mur_accel;
            end
        end
    end

    mur_l_new = mur_next;
    % fprintf('== Counter = %d ==\n', loop_counter);
    loop_counter = loop_counter + 1;
end

nbre_iter_NL_history = loop_counter;

end





% ----------------------------------------------------------------------
% Local function: Gauss-Legendre quadrature on [a,b]
% Returns points x and weights w such that ∫_a^b f(x) dx ≈ Σ w_i f(x_i)
% ----------------------------------------------------------------------
function [x, w] = gauss_legendre(n, a, b)
    % Based on the standard Golub-Welsch algorithm
    % (taken from many public domain implementations)
    if n < 1
        error('n must be >= 1');
    end

    % Legendre polynomial coefficients
    beta = 0.5 ./ sqrt(1 - (2*(1:n-1)).^(-2));
    T = diag(beta, 1) + diag(beta, -1);
    [V, D] = eig(T);
    x = diag(D);                     % eigenvalues (roots)
    w = 2 * V(1,:).^2;               % weights

    % Sort points and weights
    [x, idx] = sort(x);
    w = w(idx);

    % Map from [-1,1] to [a,b]
    x = (b - a)/2 * x + (a + b)/2;
    w = (b - a)/2 * w;
end