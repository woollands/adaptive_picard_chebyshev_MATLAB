% Robyn Woollands 2017
% Texas A&M University - Department of Aerospace Engineering
% File name     : picard_error_feedback.m
% Description   : Generates linear correction term for acceleration
% Date Written  : Apr 29, 2017
% Date Modified : Apr 29, 2017
%
% Input:  X     -- State: position (km) and velocity (km/s)
%         mu    -- Gravitational parameter (km^3 / s^2)
%         del_X -- Position Error (km)
%         del_V -- Velocity Error (km/s)
%
% Output: del_A -- Acceleration correction (km/s^2)
%================================================================

function del_A = picard_error_feedback(X,mu,del_X,del_V)

R3   = sqrt(X(:,1).^2 + X(:,2).^2 + X(:,3).^2).^3;
R5   = sqrt(X(:,1).^2 + X(:,2).^2 + X(:,3).^2).^5;

for i=1:length(X(:,1))
    
    % Two-body Jacobian
    J(1:3,1:3) = zeros(3,3);
    J(1:3,4:6) = eye(3,3);
    J(4,1)     = (3*mu*X(i,1)^2)/R5(i) - mu/R3(i);
    J(4,2)     = (3*mu*X(i,1)*X(i,2))/R5(i);
    J(4,3)     = (3*mu*X(i,1)*X(i,3))/R5(i);
    J(5,1)     = (3*mu*X(i,2)*X(i,1))/R5(i);
    J(5,2)     = (3*mu*X(i,2)^2)/R5(i) - mu/R3(i);
    J(5,3)     = (3*mu*X(i,2)*X(i,3))/R5(i);
    J(6,1)     = (3*mu*X(i,3)*X(i,1))/R5(i);
    J(6,2)     = (3*mu*X(i,3)*X(i,2))/R5(i);
    J(6,3)     = (3*mu*X(i,3)^2)/R5(i) - mu/R3(i);
    J(4:6,4:6) = zeros(3,3);
    
    % Acceleration Error Feedback
    del_A(i,:) = J(4:6,1:3)*del_X(i,:)';
    
end

return