% Robyn Woollands 2017
% Texas A&M University - Department of Aerospace Engineering
% File name     : perturbed_gravity.m
% Description   : Computes a variable fidelity force model
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
% Input:  X  -- States: position (m) and velocity (m/s)
%
% Output: G  -- Gravity (m/s^2)
%================================================================

function G = grav_approx(X)

global Jeval

G   = grav_j2j6(X(:,1),X(:,2),X(:,3));

Jeval = Jeval + 1;

return
