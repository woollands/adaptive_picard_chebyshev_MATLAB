% Robyn Woollands 2017
% Texas A&M University - Department of Aerospace Engineering
% File name     : grav_var_fid.m
% Description   : Computes variable fidelity force model
% Date Written  : April 4, 2014
% Date Modified : April 29, 2017
% References
%
% 1. Macomber, B., Probe, A., Woollands, R., Read, J. and Junkins, J.,
%    "Enhancements of Modified Chebyshev Picard Iteration Efficiency for
%    Perturbed Orbit Propagation", Computational Modelling in Engineering &
%    Sciences, Vol. 111, 2016.
%
% 2. Macomber, B., "Enhancements of Chebyshev-Picard Iteration Efficiency
%    for Generally Perturbed Orbits and Constrained Dynamics Systems", PhD
%    Dissertation, Texas A&M University, College Station, TX, 2015.
%
% Input:  X     -- State: position (km) and velocity (km/s)
%         errX  -- Picard integration error
%
% Output: G     -- Gravity (km/s^2)
%         model -- which model (1, 2 or 3)
%================================================================

function [G,model] = grav_var_fid(X,errX)

global itr tol hot grav_vec
persistent del_G ITR1 ITR2 ITR3 ITR4

% Initialization
if itr == 0 && hot == 0
    ITR1 = 0;
    ITR2 = 0;
    ITR3 = 0;
    ITR4 = 0;
end

% Initialization with Hot Start
if itr == 0 && hot == 1
    ITR1 = -1;
    ITR2 = -1;
    ITR3 = -1;
    ITR4 = -1;
end

if errX > 1e-1
%     disp('J2-J6')
    G       = grav_approx(X);
    ITR1    = itr;
    ITR2    = itr;
    ITR3    = itr;
    ITR4    = itr;
    model   = 1;
elseif errX <= 1e-1 && errX > 1e-4 && ITR1 == itr - 1 % 1e-1 1e-4
%     disp('FULL Gravity 1')
    G       = grav_full(X);
    del_G   = G - grav_approx(X);
    model   = 3;
elseif errX <= 1e-1 && errX > 1e-4                     % 1e-1 1e-4
%     disp('Approx Gravity 1')
    G       = grav_approx(X) + del_G;
    ITR2    = itr;
    ITR3    = itr;
    ITR4    = itr;
    model   = 2;
elseif errX <= 1e-4 && errX > 1e-7 && ITR2 == itr - 1  % 1e-4 1e-7
%     disp('FULL Gravity 2')
    G       = grav_full(X);
    del_G   = G - grav_approx(X);
    model   = 3;
elseif errX <= 1e-4 && errX > 1e-7                     % 1e-4 1e-7
%     disp('Approx Gravity 2')
    G       = grav_approx(X) + del_G;
    ITR3    = itr;
    ITR4    = itr;
    model   = 2;
elseif errX <= 1e-7 && errX > 1e-10 && ITR3 == itr - 1 % 1e-7 1e-10
%     disp('FULL Gravity 3')
    G       = grav_full(X);
    del_G   = G - grav_approx(X);
    model   = 3;
elseif errX <= 1e-7 && errX > 1e-10                    % 1e-7 1e-10
%     disp('Approx Gravity 3')
    G       = grav_approx(X) + del_G;
    ITR4    = itr;
    model   = 2;
elseif errX <= 1e-10 && errX > 1e-12 && ITR4 == itr - 1 % 1e-10 1e-12
%     disp('FULL Gravity 4')
    G       = grav_full(X);
    del_G   = G - grav_approx(X);
    model   = 3;
elseif errX > tol
%     disp('Approx Gravity 4')
    G       = grav_approx(X) + del_G;
    model   = 2;
end
itr = itr + 1;

return

