% Robyn Woollands 2016
% Texas A&M University - Department of Aerospace Engineering
% File name     : ecef2eci.m
% Description   : Converts acceleration from body to inertial frame
% Date Written  : June 15, 2016
% Date Modified : June 15, 2016
%
% Input:  t     -- time (s)
%         aB    -- ECEF Position (km/s^2)
%         omega -- Earth rotational velocity (rad/s)
%
% Output: acc   -- ECI Acceleration (km/s^2)  
%================================================================

function acc = ecef2eci(t,aB,omega)

th           = t*omega;
cos_th       = cos(th);
sin_th       = sin(th);

% Convert to Acceleration to Inertial Frame
acc(:,1)     = cos_th.*aB(:,1) - sin_th.*aB(:,2);
acc(:,2)     = sin_th.*aB(:,1) + cos_th.*aB(:,2);
acc(:,3)     = aB(:,3);

return