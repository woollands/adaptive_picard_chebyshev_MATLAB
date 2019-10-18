% Texas A&M University - Department of Aerospace Engineering
% Author        : Robyn Woollands
% File name     : exp_drag.m
% Description   : Computes atmospheric drag
% Date Written  : April 4, 2014
% Date Modified : February 18, 2017
% References    : D. Vallado, "Fundamentals of Astrodynamics and Applications", 4th ed, 2013.
%
% Input:  xB   -- ECEF position (km)
%         vB   -- ECEF velocity (km/s)
%
% Output: drag -- Atmospheric Drag (m/s^2)
%================================================================

function drag = exp_drag(xB,vB,F,Re)

% Table 8.4: Vallado
% Altitude (km) rho0 (kg/m^3) Scale Height (km)
data = [0	1.225	    7.249;
    25	3.90E-02	6.349;
    30	1.77E-02	6.682;
    40	3.97E-03	7.554;
    50	1.06E-03	8.382;
    60	3.21E-04	7.714;
    70	8.77E-05	6.549;
    80	1.91E-05	5.799;
    90	3.40E-06	5.382;
    100	5.30E-07	5.877;
    110	9.66E-08	7.263;
    120	2.44E-08	9.473;
    130	8.48E-09	12.636;
    140	3.85E-09	16.149;
    150	2.07E-09	22.523;
    180	5.46E-10	29.74;
    200	2.79E-10	37.105;
    250	7.25E-01	45.546;
    300	2.42E-11	53.628;
    350	9.52E-12	53.298;
    400	3.73E-12	58.515;
    450	1.59E-12	60.828;
    500	6.97E-13	63.822;
    600	1.45E-13	71.835;
    700	3.61E-14	88.667;
    800	1.17E-14	124.64;
    900	5.25E-15	181.05;
    1000 3.02E-15	268];

% Compute Altitude
for i = 1:length(xB)
    rmag(i)     = norm(xB(i,:))-Re;
    param(i,:)  = data(max(find(data(:,1) < rmag(i))),:);
end

% Compute Drag
vB = vB*1e3;        % Convert to m/s
for i = 1:length(xB)
    V(i,1)   = norm(vB(i,:)); % m/s
    Vsq(i,1) = V(i)^2;        
end

% Atmospheric Density
rho  = param(:,2).*exp(-(rmag'-param(:,1))./param(:,3));    % kg/m^3

% Compute Drag
drag(:,1) = -F.*rho.*Vsq.*vB(:,1)./V;
drag(:,2) = -F.*rho.*Vsq.*vB(:,2)./V;
drag(:,3) = -F.*rho.*Vsq.*vB(:,3)./V;

return



