% Texas A&M University - Department of Aerospace Engineering
% Robyn Woollands 2017
% Texas A&M University - Department of Aerospace Engineering
% File name     : grav_full.m
% Description   : Computes spherical harmonic gravity model
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
% 3. Junkins, J., Bani Younes, A., Woollands, R., and Bai, X., “Efficient 
%    and Adaptive Orthogonal Finite Element Representation of the Geopotential”,
%    Journal of the Astronautical Sciences, January, 2017.
%
% Input:  X  -- States: position (m) and velocity (m/s)
%
% Output: G  -- Gravity (m/s^2)
%================================================================

function G = grav_full(X)

global grav_vec Feval Deg

% Compute Spherical Harmonic Gravity
for i = 1:length(X(:,1))
    
    G(i,:) = EGMGravMex([X(i,1), X(i,2), X(i,3)],grav_vec(i));
    
    % Count Full Force Evaluations
    Feval = Feval + grav_vec(i)^2 / Deg^2;
end

return
