% Robyn Woollands 2017
% Texas A&M University - Department of Aerospace Engineering
% File name     : inform_user.m
% Description   : Prints information to screen
% Date Written  : March 27, 2017
% Date Modified : March 27, 2017

disp('=============================================')
disp('     ADAPTIVE PICARD-CHEBYSHEV ITERATION')
disp('        (perturbed two-body problem)')
disp('            TEXAS A&M UNIVERSITY')
disp('                   (2017)')
disp('=============================================')
fprintf('\n')
disp('                    INPUT')
disp('---------------------------------------------')
disp(['Orbit Elements: a=',num2str(a), 'km, ',...
    'e=',num2str(e),', ',...
    'i=',num2str(inc*180/pi),' (deg), '])
disp(['Omega=',num2str(Om*180/pi),' (deg), ',...
    'w=',num2str(w*180/pi),' (deg), ',...
    'M0=',num2str(M0*180/pi),' (deg).'])
fprintf('\n')
disp(['Time: t0 = ',num2str(t0),' (s), tf = ',num2str(tf),' (s).'])
disp(['Parameters: Gravity Degree = ',num2str(Deg),', tol = ',num2str(tol),','])
if atm_drag == 1
    disp('Atmospheric Drag = ON.')
else
    disp('Atmospheric Drag = OFF.')
end
fprintf('\n')
disp('                  ALGORITHM')
disp('---------------------------------------------')
if var_fid == 1
    disp('Variable Fidelity Gravity Model = YES')
else
    disp('Grav                            = NO')
end
if rad_grav == 1
    disp('Radially Adaptive Gravity       = YES')
else
    disp('Radially Adaptive Gravity       = NO')
end
if quasi_lin == 1
    disp('Quasilinearization              = YES')
else
    disp('Quasilinearization              = NO')
end
if hot_start == 1
    disp('Hot Start (1+ orbits)           = YES')
else
    disp('Hot Start (1+ orbits)           = NO')
end
fprintf('\n')
disp('                    OUTPUT')
disp('---------------------------------------------')
disp(['Segments per orbit    = ',num2str(seg)])
disp(['Nodes per segment     = ',num2str(N),])
disp(['Function Evaluations  = ',num2str(Teval)])
disp(['Max Hamiltonian Error = ',num2str(max(H))])
disp('=============================================')