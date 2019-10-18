% Robyn Woollands 2016
% Texas A&M University - Department of Aerospace Engineering
% File name     : orb_feasibility.m
% Description   : Convert Cartesian to Keplerian elements
% Date Written  : May, 2017
% Date Modified : May, 2017
%
% Inputs: r   -- Cartesian Position (km)
%         v   -- Cartesian Velocity (km/s)
%         mu  -- Gravitational Parameter (km^3 / s^2)


function res = orb_feasibility(r,v)

global mu Re

res = 1;

% Check for non-elliptic orbits
if norm(v) > 11.2
    disp(['### ERROR: Orbit not Elliptic (v = ',num2str(norm(v)),' km/s > escape speed)'])
    res = 0;
end

% Check for Earth Collisions
elm = rv2elm(r,v,mu,1e-10);
rp = elm(2)*(1 - elm(3));
if rp < 1.05*Re
    disp(['### ERROR: Earth Collision (rp = ',num2str(rp),' km)'])
    res = 0;
end

return