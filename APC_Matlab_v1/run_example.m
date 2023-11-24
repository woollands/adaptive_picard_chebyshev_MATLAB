% VERSION 1
% Texas A&M University - Department of Aerospace Engineering
% Author        : Robyn Woollands
% File name     : run_example.m
% Description   : Sets up initial conditions for propagating a perturbed
% two-body orbit and integrates using adaptive_picard.m to find the
% solution.
% Date Written  : May, 2017
% Date Modified : Oct, 2019
%================================================================

clear
close all
clc

% Specifiy Constants
const;

%% BEGIN USER INPUT

%%% Specify Orbit Elements OR Cartesian State %%%
% Orbit Elements
a           = 6878+550;             % Semimajor Axis (km)
e           = 0.0;                 % Eccentricity
inc         = 7*pi/180;             % Inclination (radians)
w           = 0*pi/180;             % Argument of Perigee (radians)
Om          = 0*pi/180;             % Right Ascension of Ascending Node (radians)
p           = a*(1 - e^2);          % Semilatus Rectum (km)
M0          = 0*pi/180;             % Mean Anomaly (radians)
s           = 0*pi/180;             % "Special Parameter" (see Vallado) (radians)
Period      = 2*pi*sqrt(a^3/mu);    % Period (s)

% Cartesian State
r0          = [];                   % Initial Position Vector (km)
v0          = [];                   % Initial Velocity Vector (km/s)

% Specify Time Interval
t0          = 0;                    % Initial Time (s)
tf          = 10*Period;             % Final Time (s)
dt          = 60;                   % Output Time Interval (s)

% Specify Solution Accuracy
Deg         = 70;                   % Degree of Gravity Model (max = 100)
tol         = 1e-15;                % Tolerance (max = 1e-15)

% Drag Parameters
m           = 295;                 % Mass (kg)
Cd          = 2;                    % Drag Coefficient
A           = 2;                   % Area (m^2)

% Specify Algorithm Specifications
var_fid     = 1;                    % Variable Fidelity Gravity Model (YES = 1, NO = 0)
rad_grav    = 1;                    % Radially Adaptive Gravity (YES = 1, NO = 0)
quasi_lin   = 0;                    % Acceleration via Quasilinearization (YES = 1, NO = 0)
hot_start   = 0;                    % Hot Start for subsequent orbits (YES = 1, NO = 0)
atm_drag    = 0;                    % Exponential Drag (see Vallado) (YES = 1, NO = 0)

% Specify Output
inform      = 0;                    % Print the output information to screen (YES = 1, NO = 0)
show_figs   = 0;                    % Plot trajectory and Hamiltonian (YES = 1, NO = 0)

%%% END USER INPUT

%% Propagate using Adaptive-Picard-Chebyshev

% Compute r & v from elements
if isempty(r0) && isempty(v0) == 1
    [r0,v0] = elm2rv(a,e,inc,Om,w,M0,s,mu);
end

% Check Orbit Feasibility (Earth collisions and non-elliptic orbits)
res = orb_feasibility(r0,v0);
if res == 0
    return
end

% Build Input Structure
input = struct('r0',r0,'v0',v0,'t0',t0,'tf',tf,'dt',dt);
    
% Build Physical Parameters Structure
if atm_drag == 0
    params = struct('Deg',Deg,'tol',tol);
elseif atm_drag == 1
    params = struct('Deg',Deg,'tol',tol,...
        'm',m,'Cd',Cd,'A',A);
end

% Build Algorithm Specfics Structure
specs = struct('var_fid',var_fid,...
    'rad_grav',rad_grav,...
    'quasi_lin',quasi_lin,...
    'hot_start',hot_start,...
    'atm_drag',atm_drag,...
    'show_figs',show_figs);

% Call adaptive algorithm (Junkins & Woollands, JGCD, Jan 2019: https://arc.aiaa.org/doi/10.2514/1.G003318)
tic
[time,X,V,seg,N,Teval] = adaptive_picard_chebyshev(input,params,specs);
toc

%% Output

% Compute elements from r & v
if isempty(r0) && isempty(v0) == 0
    elm    = rv2elm(r,v,mu,tol);
    Period = 2*pi*sqrt(elm(2)^3/mu);    % Period (s)
end

% Compute Hamiltonian
H = jacobi_integral(time,X,V,Period,1);

% Inform User 
if inform == 1
    inform_user
end
