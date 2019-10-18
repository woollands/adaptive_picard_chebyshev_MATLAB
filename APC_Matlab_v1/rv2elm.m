% Robyn Woollands 2016
% Texas A&M University - Department of Aerospace Engineering
% File name     : rv2elm.m
% Description   : Convert Cartesian to Keplerian elements
% Date Written  : June 15, 2016
% Date Modified : June 15, 2016
% Reference     : Vallado (p. 120 , Algorithm 9)
%
% Inputs: r   -- Cartesian Position (km)
%         v   -- Cartesian Velocity (km/s)
%         mu  -- Gravitational Parameter (km^3 / s^2)
%
% Outputs: p  -- Semilatus Rectum (km)
%          a  -- Semimajor Axis (km)
%          e  -- Eccentricity
%          i  -- Inclination (rad)
%          Om -- Right Ascension of Ascending Node (rad)
%          w  -- Argument of Perigee (rad)
%          f  -- True Anomaly (rad)
%          E  -- Eccentric Anomaly (rad)
%          M  -- Mean Anomaly (rad)
%          s  -- Special case location of perigee (rad)
%                    -- Longitude of Perigee
%                    -- Argument of Latitude
%                    -- True Longitude 
% ================================================================

function elm = rv2elm(r,v,mu,tol)

% Default Tolerance
% tol     = 1e-10;

% Correct Vector Dimensions
if isrow(r) ~= 1
    r = r';
end
if isrow(v) ~= 1
    v = v';
end

% Position & Velocity Magnitudes
R       = norm(r);
V       = norm(v);

% Angular Momentum Vector
h       = cross(r,v);
H       = norm(h);

% Line of Nodes Vector
nvec    = cross([0 0 1]',h);
n       = norm(nvec);

% Eccentricity Vector
evec    = 1/mu*((V^2 - mu/R)*r - (r*v')*v);
e       = norm(evec);

% Energy
xi      = (V^2)/2 - mu/R;

% Semimajor Axis (a) & Semillatus Rectum (p)
if abs(1-e) < tol
    a   = Inf;
    p   = (H^2)/mu;
else
    a   = -mu/2/xi;
    p   = a*(1 - e^2);
end

% Inclination
i       = acos(h(3)/H);

% Right Ascension of Ascending Node
Om      = acos(nvec(1)/n);
if (nvec(2) < 0)
    Om = 2*pi - Om;
end

% Argument of Perigee
w = acos((nvec*evec')/n/e);
if (evec(3) < 0)
    w = 2*pi - w;
end

% True Anomaly
f = real(acos((evec*r')/R/e));

if (r*v' < 0)
    f = 2*pi - f;
end

% Mean Anomaly & Eccentric Anomaly
E = 2*atan2(sqrt(1-e)*tan(f/2),sqrt(1+e));
if E < 0
    E = 2*pi + E;
end
M = E - e*sin(E);
if M < 0
    M = 2*pi + M;
end

% Special Cases
% Initialize s
s = 0;

% Elliptical Equatorial (ascending node undefined)
if (i < tol) && (e >= tol)
    s = acos(evec(1)/e);
    if (evec(2) < 0)
        s = 2*pi - s;   % Longitude of Perigee
    end
% Circular Inclined (perigee undefined)
elseif (i >= tol) && (e < tol)
    s = acos((nvec*r')/R/n);    % Argument of Latitude
    if (r(3) < 0)
        s = 2*pi - s;
    end
% Circular Equatorial (perigee & ascending node undefined)
elseif (i < tol) && (e < tol)
    s = acos(r(1)/R);
    if (r(2) < 0)
        s = 2*pi - s;    % True Longitude
    end
end

% Output
elm = [p a e i Om w f E M s];

end

