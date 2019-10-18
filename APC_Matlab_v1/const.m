% Robyn Woollands 2017
% Texas A&M University - Department of Aerospace Engineering
% File name     : const.m
% Description   : Defines constants
% Date Written  : May, 2017
% Date Modified : May, 2017
%================================================================

global mu muCan Re omega DU TU VU AU

mu          = 398600.4418;          % Gravitational Parameter
muCan       = 1;                    % Non-dimensional Gravitational Parameter
Re          = 6378.137;             % Earth Radius
omega       = 0;%7.2921151e-5;         % Earth Rotation Rate

% Canonical Units
DU          = 6378.137;             % Distance Unit
TU          = sqrt(DU^3 / mu);      % Time Unit
VU          = DU/TU;                % Velocity Unit
AU          = DU/TU/TU;             % Acceleration Unit