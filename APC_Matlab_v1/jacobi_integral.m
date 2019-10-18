%================={ Copyright (c) Ahmad Bani-Younes 2013 }=================
%                Texas A&M University - AEROSPACE department
% 
% File name     : Jacobi_Integral.m
% Subject       : Computes the Hamiltonian.
% Description   : To validate orbit propagator accuracy.
%    
% Sub-files     : egm2008GPsphericalharmonic.m
%
% Compiler      : MATLAB 7.11.0 (R2010b)
% Date Modified : 03/27/2017 (Robyn Woollands)
%==========================================================================

function H = jacobi_integral(time,X,V,Period,show_figs)

global GM Re omega Deg C S

load('aeroegm2008.mat') % [GM, Re, degree, C, S]

% ECI to ECEF
[xB,vB] = eci2ecef(time,X,V,omega);

% Convert to m and m/s
xB     = xB.*1e3;
vB     = vB.*1e3;

term   = 0.5*omega*omega.*(xB(:,1).^2 + xB(:,2).^2);
KE     = 0.5*( vB(:,1).^2 + vB(:,2).^2 + vB(:,3).^2 );
V      = -egm2008GPsphericalharmonic(xB, Deg) ;
Enrgy  = V + KE - term;
Enrgyo = Enrgy(1); 
dEnrgy = abs((Enrgy - Enrgyo)/Enrgyo);
H      = dEnrgy;

if show_figs == 1
    figure;
    semilogy(time./Period,dEnrgy,'ro-','linewidth',2);
    grid on
    set(gca, 'FontName', 'Helvetica','FontSize',16)
    title('Hamiltonian')
    xlabel('Time (Periods)')
    ylabel('| E - E_o | / | E_o |')
end

return