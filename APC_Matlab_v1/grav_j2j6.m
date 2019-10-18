% Texas A&M University - Department of Aerospace Engineering
% Author        : Robyn Woollands
% File name     : grav_j2j6.m
% Description   : Computes acceleration for two-body + J2 to J6 zonals
% Date Written  : March 20, 2014
% Date Modified : April 29, 2017
% References
%
% 1. Schaub, H., and Junkins, J., "Analytical Mechanics of Space Systems",
%    2nd ed., 2009.
%
%================================================================

% Robyn Woollands 2017
% Texas A&M University - Department of Aerospace Engineering
% File name     : grav_j2j6.m
% Description   : Computes acceleration for two-body + J2 to J6 zonals
% Date Written  : March 20, 2014
% Date Modified : April 29, 2017
% References
%
% 1. Schaub, H., and Junkins, J., "Analytical Mechanics of Space Systems",
%    2nd ed., 2009.
%
% Input:  x     -- Cartesian state x (m)
%         y     -- Cartesian state y (m)
%         z     -- Cartesian state z (m)
%
% Output: aJ2J6 -- Acceleration (first 6 zonals) (m/s^2)
%================================================================


function [aJ2J6] = grav_j2j6(x,y,z)

GM      = 3.986004418e14;   % Gravitational Parameter (m^3 / s)    
Req     = 6378137;          % Earth Radius (m)

r       = sqrt(x.^2 + y.^2 + z.^2);

% Definitions
J2      = 1082.63e-6;
J3      = -2.52e-6;
J4      = -1.61e-6;
J5      = -0.15e-6;
J6      = 0.57e-6;

x_r_1   = x./r;
y_r_1   = y./r;
z_r_1   = z./r;
z_r_2   = (z./r).^2;
z_r_3   = (z./r).^3;
z_r_4   = (z./r).^4;
z_r_5   = (z./r).^5;
z_r_6   = (z./r).^6;

% Two-body Acceleration
aTB    = [-(GM./r.^3).*x, -(GM./r.^3).*y, -(GM./r.^3).*z];

% Zonal Acceleration (J2-J6)
aJ2     = [(-3./2.*J2.*(GM./r.^2).*(Req./r).^2).*[(1 - 5.*z_r_2).*x_r_1],...
    (-3./2.*J2.*(GM./r.^2).*(Req./r).^2).*[(1 - 5.*z_r_2).*y_r_1],...
    (-3./2.*J2.*(GM./r.^2).*(Req./r).^2).*[(3 - 5.*z_r_2).*z_r_1]];

aJ3     = [(1./2.*J3.*(GM./r.^2).*(Req./r).^3).*[5.*(7.*z_r_3 - 3.*z_r_1).*x_r_1],...
    (1./2.*J3.*(GM./r.^2).*(Req./r).^3).*[5.*(7.*z_r_3 - 3.*z_r_1).*y_r_1],...
    (1./2.*J3.*(GM./r.^2).*(Req./r).^3).*[3.*(1 - 10.*z_r_2 + 35./3.*z_r_4)]];

aJ4     = [(5./8.*J4.*(GM./r.^2).*(Req./r).^4).*[(3 - 42.*z_r_2 + 63.*z_r_4).*x_r_1],...
    (5./8.*J4.*(GM./r.^2).*(Req./r).^4).*[(3 - 42.*z_r_2 + 63.*z_r_4).*y_r_1],...
    (5./8.*J4.*(GM./r.^2).*(Req./r).^4).*[(15 - 70.*z_r_2 + 63*z_r_4).*z_r_1]];

aJ5     = [(J5./8.*(GM./r.^2).*(Req./r).^5).*[3*(35.*z_r_1 - 210.*z_r_3 + 231.*z_r_5).*x_r_1],...
    (J5./8.*(GM./r.^2).*(Req./r).^5).*[3*(35.*z_r_1 - 210.*z_r_3 + 231.*z_r_5).*y_r_1],...
    (J5./8.*(GM./r.^2).*(Req./r).^5).*[693.*z_r_6 - 945.*z_r_4 + 315.*z_r_2 - 15]];

aJ6     = [(-J6./16.*(GM./r.^2).*(Req./r).^6).*[(35 - 945.*z_r_2 + 3465.*z_r_4 - 3003.*z_r_6).*x_r_1],...
    (-J6./16.*(GM./r.^2).*(Req./r).^6).*[(35 - 945.*z_r_2 + 3465.*z_r_4 - 3003.*z_r_6).*y_r_1],...
    (-J6./16.*(GM./r.^2).*(Req./r).^6).*[(245 - 2205.*z_r_2 + 4851.*z_r_4 - 3003.*z_r_6).*z_r_1]];

% Two-body + Zonals
aJ2J6   = aTB + aJ2 + aJ3 + aJ4 + aJ5 + aJ6;

return