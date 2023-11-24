% VERSION 1
% Texas A&M University - Department of Aerospace Engineering
% Author        : Robyn Woollands
% File name     : adaptive_picard_chebyshev.m
% Description   : Solves perturbed two-body problem with automated adaptive
% tuning and error feedback
% Date Written  : Mar 8, 2017
% Date Modified : Oct 17, 2019
% References    : R. Woollands and J. Junkins, "Nonlinear Differential
% Equations Solvers via Adaptive Picard-Chebyshev Iteration: Applications
% in Astrodynamics", Journal of Guidance, Control and Navigation, Jan 2019.

% Inputs: input  -- Struct containing input parameters
%         params -- Struct containing phsyical parameters
%         specs  -- Struct containing algorithm specifics
%
% Outputs: Tout  -- Times (s)
%          Xout  -- Position Solution (km)
%          Vout  -- Velocity Solution (km/s)
%          seg   -- Number of segments per orbit
%          N     -- Number of nodes per segment
%          Teval -- Total number of equivalent function evaluations (see paper)

%==========================================================================

function [Tout,Xout,Vout,seg,N,Teval] = adaptive_picard_chebyshev(input,params,specs)

global Re C S mu GM omega Deg tol itr grav_vec Jeval Feval hot DU TU VU AU

%% PART 1: Assign Variables, Compute Tolerances & Set Counters

% Inputs
r0          = input.r0;
v0          = input.v0;
t0          = input.t0;
t_final     = input.tf;
dt          = input.dt;

% Algorithm Specifics
var_fid     = specs.var_fid;
rad_grav    = specs.rad_grav;
quasi_lin   = specs.quasi_lin;
hot_start   = specs.hot_start;
atm_drag    = specs.atm_drag;
show_figs   = specs.show_figs;

% Physical Parameters
Deg         = params.Deg;
tol         = params.tol;

if atm_drag == 1
    m       = params.m;
    Cd      = params.Cd;
    A       = params.A;
    Dconst  = 0.5*Cd*A/m;
end

% Tolerance Computation
% coeff_tol   = tol/100;
% fit_tol     = tol/10;
coeff_tol   = tol/5;
fit_tol     = tol;

% Function Evaluation Counters
Jeval       = 0;  % Approximate Gravity Model
Feval       = 0;  % Full Gravity Model

if show_figs == 1
    color = {'k-o','m-o','c-o','g-o','r-o','b-o','y-o',...
        'y-+','m-*','c-*','g-*','r-*','b-*',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+',...
        'y-+','m-+','c-+','g-+','r-+','b-+'};
end
%% PART 2: Determine Segmentation Scheme

% In PART 2 we determine how many segments per orbit and how many nodes per
% segment are required to produce a solution that satisfies the user
% desired tolerance giving the minimum overall number of function
% evaluations. See README.txt for more details.

% Compute Keplerian Orbit Period
elm     = rv2elm(r0,v0,mu,tol);
a       = elm(2);
e       = elm(3);
rp      = a*(1 - e);
Period  = 2*pi*sqrt(a^3 / mu);

% Compute Time Vector for 1 Orbit Period
w1      = Period/2;
w2      = Period/2;
tau     = -cos((0:100)*pi/100);
time    = w2*tau + w1;
% Compute F&G Analytical Solution
for i = 1:length(time)
    [r, v]   = FnG(0,time(i),r0,v0,mu);
    X(i,:)   = r';
    V(i,:)   = v';
    normX(i) = norm(X(i,:));
end

% Compute State and Time at Keplerian Perigee
ind     = min(find(abs(normX'-min(normX)) < 1e-10));
rp      = X(ind,:);
vp      = V(ind,:);
tp      = time(ind);

% If user specified ICs were not at perigee (i.e. ind~=1) then we need a
% finer grid to find perigee. I put 100 sample points in the vacinity of
% the perigee point and search again for a more exact value.
if ind ~= 1
    time_test = linspace(time(ind-1),time(ind+1),100);
    for i = 1:length(time_test)
        [r, v]     = FnG(time_test(1),time_test(i),X(ind-1,:),V(ind-1,:),mu);
        Xtest(i,:) = r';
        Vtest(i,:) = v';
        normX(i)   = norm(Xtest(i,:));
    end
    
    % Perigee
    ind = min(find(abs(normX'-min(normX)) < 1e-10));
    rp  = Xtest(ind,:);
    vp  = Vtest(ind,:);
    tp  = time_test(ind);
end

fit_check = 0;    % Loop check condition
seg       = 3;    % Minimum number of segments per orbit
coeff     = 3;    % Value of last 3 coefficients must be below tolerance

while fit_check <= coeff
    clear G GG Gprev Yp errA
    f       = linspace(0,360,seg+1).*pi/180; % Equal segments in true anomaly
    E       = 2*atan2(tan(0.5*f(2))*sqrt(1-e),sqrt(1+e)); % Eccentric Anomaly
    if E < 0
        E = 2*pi + E;
    end
    Mf      = E - e*sin(E);    % Mean Anomaly
    n       = 2*pi/Period;     % Mean Motion
    t0      = 0;               % Initial Time
    tf      = Mf/n;            % Final Time
    % Scaling parameters & time
    w1      = (t0 + tf)/2;
    w2      = (tf - t0)/2;
    
    Nint   = 10;
    jmax   = 3;
    Nvec   = [Nint 2*Nint 4*Nint 8*Nint];   % Node doubling (see paper)
    
    % Loop through different values for N
    for j  = 1:jmax
        % Initialize indices
        pcnt = 1;
        gcnt = 1;
        clear X V G tau time
        % Generate discrete cosine nodes
        if j == 1
            N    = Nvec(j);
            tau  = -cos((0:N)*pi/N);
            time = (w2*tau + w1);
        elseif j > 1
            N    = Nvec(j);
            tau  = -cos((1:2:N)*pi/N);
            time = (w2*tau + w1);
        end
        % Compute F&G Analytical Solution
        for i = 1:length(time)
            [r, v] = FnG(t0,time(i),rp,vp,mu);
            X(i,:) = r';
            V(i,:) = v';
        end
        % Compute Gravity Radial Adaptation
        grav_vec = radial_gravity(X./DU,tol,Deg,rad_grav);
        % Compute Acceleration
        for i = 1:length(time)
            G(i,:) = EGMGravMex([X(i,1), X(i,2), X(i,3)].*1e3,grav_vec(i))./1e3;
            Feval  = Feval + grav_vec(i)^2 / Deg^2;  % Function Evaluations
        end
        % Build acceleration vector (include new points as they are
        % computed)
        if j == 1
            Gprev = G;
        elseif j > 1
            for i=1:N+1
                if mod(i,2) == 1
                    GG(i,:) = Gprev(pcnt,:);
                    pcnt = pcnt + 1;
                elseif mod(i,2) == 0
                    GG(i,:) = G(gcnt,:);
                    gcnt = gcnt + 1;
                end
            end
            Gprev = GG;
        end
        
        % Least Sqaures Acceleration Approximation
        M = N;
        [A,T]    = lsq_chebyshev_fit(1,N-1,M);  % Least Squares Operator (A)
        gamma    = A*Gprev;                     % Least Squares Coefficients
        G_approx = T*gamma;                     % Approximated Acceleration
        
        % Normalized Residuals
        Err = max(max(abs(Gprev - G_approx)))./AU;
        
        % Check the magnitude of the last few coefficients
        clear max_gamma
        for i = 1:length(gamma)
            N = i;
            max_gamma = max(abs(gamma(i,:)));
            if abs(max_gamma) < coeff_tol
                fit_check = fit_check + 1;
                if fit_check == coeff
                    break
                end
            end
        end
%                 disp(['seg = ',num2str(seg),', N = ',num2str(N),', Coefficient Magnitude (4th last) = ',...
%                     num2str(max_gamma(i)),', Acceleration Fit = ',num2str(Err)])
        if fit_check == coeff
            break
        end
        fit_check = 0;
    end
    
    clear gamma
    seg = seg + 2;
    
    if fit_check == coeff
        break
    end
end

seg = seg - 2;

% Reduce N from above until least squares acceleration fit satisfies tolerance
while Err < coeff_tol
    N = N-1;
    % Least Sqaures Acceleration Approximation
    [A,T]    = lsq_chebyshev_fit(1,N,M);  % Least Squares Operator (A)
    zeta     = A*Gprev;                   % Least Squares Coefficients
    G_approx = T*zeta;                    % Approximated Acceleration
    
    % Normalized Residuals
    for i = 1:M+1
        fit_err(i) = max(abs(Gprev(i,:) - G_approx(i,:)))./AU;
    end
    Err = max(fit_err);
end

% gamma = zeta;
% 
% figure(223)
% subplot 311
% semilogy(abs(gamma(1:N,1)),'k-o','Linewidth',2)
% hold on
% grid on
% plot(1e-7.*ones(length(gamma(1:N,1)),1),'k--','Linewidth',2)
% plot(1e-15.*ones(length(gamma(1:N,1)),1),'k-.','Linewidth',2)
% ylabel('x')
% subplot 312
% semilogy(abs(gamma(1:N,2)),'k-o','Linewidth',2)
% hold on
% grid on
% plot(1e-7.*ones(length(gamma(1:N,1)),1),'k--','Linewidth',2)
% plot(1e-15.*ones(length(gamma(1:N,1)),1),'k-.','Linewidth',2)
% ylabel('y')
% legend('Magnitude','Tolerance 1x10^{-7}','Tolerance 1x10^{-15}')
% subplot 313
% semilogy(abs(gamma(1:N,3)),'k-o','Linewidth',2)
% hold on
% grid on
% plot(1e-7.*ones(length(gamma(1:N,1)),1),'k--','Linewidth',2)
% plot(1e-15.*ones(length(gamma(1:N,1)),1),'k-.','Linewidth',2)
% xlabel('Number of Nodes')
% ylabel('z')
   
% Plot Earth
if show_figs == 1
    % NOTE: The following code (17 lines) was obtained from MathWorks
    % online file exchange (Ryan Gray).
    load('topo.mat','topo','topomap1');
    colormap(topomap1);
    % Create the surface.
    [x,y,z] = sphere(50);
    props.AmbientStrength           = 0.1;
    props.DiffuseStrength           = 1;
    props.SpecularColorReflectance  = .5;
    props.SpecularExponent          = 20;
    props.SpecularStrength          = 1;
    props.FaceColor                 = 'texture';
    props.EdgeColor                 = 'none';
    props.FaceLighting              = 'phong';
    props.Cdata                     = topo;
    figure(1)
    surface(x,y,z,props);
%     set(gca,'color','black')
    axis equal
end

%% PART 3: Adaptive Picard-Chebyshev Iteration Initial Setup

% In PART 3 we compute and store the begin and end times for each segment
% (based on true anomaly segmentation). We also build the Clenshaw-Curtis
% quadrature constant matrices.

% Compute Time Vector (based on Keplerian true anomaly segments)
for i=1:seg+1
    E(i) = 2*atan2(tan(0.5*f(i))*sqrt(1-e),sqrt(1+e));  % Eccentric Anomaly
    if E(i) < 0
        E(i) = 2*pi + E(i);
    end
    Mf(i) = E(i) - e*sin(E(i));    % Mean Anomaly
end
% Time vector
t_orig = Mf./n;
tvec   = t_orig;

% Short first segment if user specified IC's do not coincide with segment
% break
ts     = Period - tp;

% Update tvec to include short first segment if required
if abs(Period - ts) > 1e-5 || abs(tp) > 1e-5
    tvec = [t0 t_orig(find(t_orig >= ts))-ts];
end

% Hot Start switch condition
prep_HS = length(tvec)-1;
if prep_HS == seg
    prep_HS = -1;
end

% User specified time vector for output
time_out = [t0:dt:t_final];

% Initialize Arrays and Variables
TimeVEC = 0;   % Initialize Time vector
ITR     = [];  % Store number of iterations for plotting
BETA    = [];  % Store Solution Velocity Coefficients
ALPHA   = [];  % Store Solution Position Coefficients
fvec    = [];  % Store True Anomaly Vector for Perigee Passage Calculation
W       = [];  % Store w1 and w2 values for each segment for interpolation
k       = 1;   % Segments per orbit counter
bloop   = 0;   % Break loop condition
cnt     = 1;   % Total Segment Counter
hot     = 0;   % Hot start condition

% Build Constant Matrices
% This can be done once, for all N and loaded to save CPU time. Not done
% here but it is done in the C code. Not concerned with "timming" for this
% MATLAB code.
[T2,P2,T1,P1,Ta,A] = clenshaw_curtis_ivpII(N,N);

%% PART 4: Adaptive Picard-Chebyshev Iteration (loop through all segments)

% In PART 4 we loop through and iterate each segment. There are two
% important points to note:
% PART 4A: This is the Picard Iteration for one segment. It includes all
% the algorithm enhancements if selected by the user in run_example.m.
% PART 4B: This is the computation and re-osculation of perigee after
% propagating each orbit.

XOUT = [];
VOUT = [];
TOUT = [];
while bloop == 0
    
    clear X V Xo Vo WS time G Gout
    
    % Compute cosine time vector for given segment
    t0   = tvec(k);
    tf   = tvec(k+1);
    if tf > t_final
        tf    = t_final;
    end
    w1   = (t0 + tf)/2;
    w2   = (tf - t0)/2;
    W    = [W; [w1 w2]];
    tau  = -cos((0:N)*pi/N);
    time = (w2*tau' + w1);
    
    % Compute F&G Analytical Solution (Warm Start)
    for i = 1:length(time)
        [r, v] = FnG(t0,time(i),r0,v0,mu);
        X(i,:) = r';
        V(i,:) = v';
    end
    WS = [X V];

    % Initial Conditions
    Xo = [X(1,:); zeros(length(X)-1,3)];
    Vo = [V(1,:); zeros(length(X)-2,3)];
    
    % Re-initialize variables
    itr    = 0;     % Picard-Chebyshev Iteration Counter
    MaxIt  = 30;    % Max Iterations
    errX   = 10;    % Initialize Error
    errVec = [];    % Store error (useful for plotting convergence)
    MODEL  = [];    % Store particular gravity model used on each iteration
    
    % Compute Hot Start (after 1+ orbits)
    if hot_start == 1
        if hot == 1
            X    = X + HotX((k-1)*(length(X))+1:k*(length(X)),:);
            V    = V + HotV((k-1)*(length(X))+1:k*(length(X)),:);
            errX = 1e-3;
        end
    end
    
    % Compute Gravity Radial Adaptation
    grav_vec = radial_gravity(X./DU,tol,Deg,rad_grav);
    
    %% PART 4A: Picard Iteration
    while(errX > tol)
        
%         R3 = sqrt(X(:,1).^2 + X(:,2).^2 + X(:,3).^2).^3;
        
        % ECI to ECEF
        [xB,vB] = eci2ecef(time,X,V,omega);
        xB_prev = xB;
        x_prev = X;
        
        % Compute Spherical Harmonic Gravity
        if var_fid == 0
            for i = 1:length(xB(:,1))
                Gout(i,:) = EGMGravMex([xB(i,1), xB(i,2), xB(i,3)].*1e3,Deg)./1e3;
                Feval     = Feval + 1;
                model     = 3;
            end
        elseif var_fid == 1
            [Gout,model] = grav_var_fid(xB.*1e3,errX);  % Variable Fidelity
            Gout         = Gout./1e3;
        end
        % Compute Atmospheric Drag
        if atm_drag == 1
            atm_drg = exp_drag(xB,vB,Dconst,Re);
            Gout    = Gout + atm_drg;
        end
        
        % ECEF to ECI
        G = ecef2eci(time,Gout,omega);
        
        % Velocity
        beta   = w2.*P1*A*G + Vo;
        Vorig  = T1*beta;
        
        % Position
        alpha  = w2.*P2*beta + Xo;
        Xorig  = T2*alpha;
        
        % ECI to ECEF
        [xB,~] = eci2ecef(time,Xorig,0,omega);
        
        % Accelerated Picard Linear Error Correction
        if quasi_lin == 1
            % Linear Error Correction Position
%             del_X = xB - xB_prev;
            del_X = Xorig - x_prev;
            
            % Linear Error Correction Acceleration
            del_a = picard_error_feedback(xB,mu,del_X,0);
%             del_a = picard_error_feedback(Xorig,mu,del_X,0);
            
            % Convert to Inertial
%             del_a = ecef2eci(time,del_a,omega);
            
            % Linear Error Correction Velocity Coefficients
            gamma = w2.*P1*A*del_a;
        end
        
        % Corrected Velocity
        if quasi_lin == 1
            Beta   = beta + gamma;
        else
            Beta = beta;
        end
        Vnew   = T1*Beta;
        
        % Corrected Position
        if quasi_lin == 1
            kappa  = w2.*P2*gamma;
            Alpha  = alpha + kappa;
        else
            Alpha = alpha;
        end
        Xnew   = T2*Alpha;
        
        % Nondimensional Error
        errX   = max(max(abs([X(:,1:3)./DU V(:,1:3)./VU] - [Xnew(:,1:3)./DU Vnew(:,1:3)./VU])));
%         if itr > 1
% %             err1 = max([max(max(abs(Alpha-Alpha_old)))./max(max(abs(Xnew))) max(max(abs(Beta-Beta_old)))./max(max(abs(Vnew)))]);
%             err1 = max(max(abs(A*G-A*G_old)))./max(max(abs(G)));
%             [errX err1]
% %             errX = err1;
%         end
        
        % Update
        X      = Xnew;
        V      = Vnew;
        Alpha_old = Alpha;
        Beta_old = Beta;
        G_old = G;
        
        % Iteration Counter
        if (itr < MaxIt)
            if var_fid == 0
                itr = itr+1;
            end
        else
            itr = itr-1;
            break;
        end
        
        % Store Convergence Error and Gravity Model
        errVec = [errVec; errX];
        MODEL  = [MODEL; model];
    end
    
    % Convergence Figure
%     figure(800)
%     semilogy([1:length(errVec)],errVec,'k-','Linewidth',2)
%     hold on; grid on;
%     for i = 1:length(MODEL)
%        if MODEL(i) == 1
%            semilogy(i,errVec(i),'bo','Markersize',10,'Linewidth',2)
%        end
%        if MODEL(i) == 2
%            semilogy(i,errVec(i),'go','Markersize',10,'Linewidth',2)
%        end
%        if MODEL(i) == 3
%            semilogy(i,errVec(i),'mo','Markersize',10,'Linewidth',2)
%        end
%     end
%     legend('Path','Zonals (J2:J6)',['Full SH',num2str(Deg)],'Approx')
%     title(['Elm: a=',num2str(a),', e=',num2str(e),', Deg=',num2str(Deg)])
%     xlabel('Iteration Number')
%     ylabel('Relative Convergence Error')
%     set(gca, 'FontName', 'Helvetica','FontSize',16)
    
    XOUT = [XOUT; X];
    VOUT = [VOUT; V];
    TOUT = [TOUT; time];
    
    if show_figs == 1
        figure(1)
        hold on
        plot3(X(:,1)./Re,X(:,2)./Re,X(:,3)./Re,color{cnt})
        plot3(X(end,1)./Re,X(end,2)./Re,X(end,3)./Re,'k.','MarkerSize',15)
        xlabel('Earth Radii')
        ylabel('Earth Radii')
        axis equal
        set(gca, 'FontName', 'Helvetica','FontSize',16)
    end
    
    % Loop exit condition
    if abs(tf - t_final)/tf < 1e-12
        bloop = 1;
    end
    
    % Prepare Hot Start
    if prep_HS == -1
        HotX((k-1)*(length(X))+1:k*(length(X)),:) = X - WS(:,1:3);
        HotV((k-1)*(length(X))+1:k*(length(X)),:) = V - WS(:,4:6);
        if k == seg
            hot = 1;
        end
    end
    
    % Assign new initial conditions (next segment)
    r0 = X(end,:);
    v0 = V(end,:);
    
    % Compute True Anomaly (along trajectory)
    for i = 1:length(X)
        elm     = rv2elm(X(i,:),V(i,:),mu,tol);
        fvec(i) = elm(7);
    end
    
    %% PART 4B: Perigee Passage
    % Re-osculate Keplerian perigee after each orbit propagation. If this
    % is not done then the precomputed segment times no longer align with
    % true segment breaks (due to perturbations) and the precomputed segment
    % scheme fails to a solution that satisfies the required tolerance. This
    % effect increases with increasing eccentricity.
    peri_check = 0;
    if k == seg || k == prep_HS% && abs(tf - t_final)/tf < tol %abs(tf) > tol
        if abs(e) > 1e-5    % Skip for zero eccentricity (no need to re-osculate as perigee is undefined)
            for i = 1:length(fvec)
                % If passing through perigee
                if i > 2 && fvec(i) < fvec(i-1)

                    % Prepare for secant method
                    t1        = time(i-1,:);
                    t2        = time(i,:);
                    w1        = (time(end) + time(1))/2;
                    w2        = (time(end) - time(1))/2;
                    TAU_old   = (t1 - w1)/w2;
                    TAU       = (t2 - w1)/w2;
                    f_old     = fvec(i-1);
                    f_new     = fvec(i);
                    
                    % Initialize
                    itrf = 0;
                    err  = 10;
                    
                    % Secant method to find tau at f0
                    while err > tol
                        
                        for kk = 0:length(X)-1
                            TA(:,kk+1) = real(cos(kk*acos(TAU)));
                        end
                        for kk = 0:length(X)-2
                            TB(:,kk+1) = real(cos(kk*acos(TAU)));
                        end
                        r0 = TA*Alpha;
                        v0 = TB*Beta;
                        tf = TAU*w2 + w1;
                        
                        elm      = rv2elm(r0,v0,mu,tol);
                        f_new    = elm(7);
                        
                        if f_new > pi
                            f_new = f_new - 2*pi;
                        end
                        
                        % Compute new TAU
                        df_dtau = (f_new - f_old)/(TAU - TAU_old);
                        if df_dtau == 0
                            break;
                        end
                        TAU_new = TAU + (0 - f_new)/df_dtau;
                        
                        % Check error
                        err     = abs(f_new);
                        
                        % Update
                        f_old   = f_new;
                        TAU_old = TAU;
                        TAU     = TAU_new;
                        
                        itrf = itrf + 1;
                        if itrf > 10
                            break
                        end
                        
                    end
                    
                    % Compute time vector for next orbit
                    tvec = t_orig + tf;
                    k    = 0;
                    
                    peri_check = 1;
                    break
                end
            end
        end
        
        % Compute time vector for next orbit
        if peri_check == 0
            tvec = t_orig + tf;
        end
        
        % Reset counters
        prep_HS = -1;
        k       = 0;
    end
    
    
    % Segments per orbit counter
    k = k+1;
    
    % Store Trajectory Coefficients, Time and Iterations
    BETA    = [BETA; Beta];
    ALPHA   = [ALPHA; Alpha];
    TimeVEC = [TimeVEC; tf];
    ITR     = [ITR; itr];
    
    % Total segments counter
    cnt = cnt + 1;
    
end

% Compute Total Equivalent Function Evaluations
Jeval = Jeval*(N+1)*(6*6)/(Deg^2);
Teval = ceil(Feval+Jeval);

%% PART 5: Interpolate Solution

% IN PART 5 the Chebyshev coefficients from each of the orbit segments are
% used to compute the solution (r and v) at the user desired times.

Xout = [];
Vout = [];
Tout = [];
% Loop through all segments
for i = 1:length(TimeVEC)-1
    
    clear Tv Tx tau tt
    
    tt = [];
    % User desired output times for a given segment
    for j = 1:length(time_out)
        if time_out(j) >= TimeVEC(i) && time_out(j) < TimeVEC(i+1)
            tt = [tt; time_out(j)];
        end
    end
    
    if isempty(tt) == 0
        % Taus corresponding to user desired output times
        w1    = W(i,1);
        w2    = W(i,2);
        tau   = (tt - w1)./w2;

        % Chebyshev velocity matrix
        for kk = 0:N-1
            Tv(:,kk+1) = real(cos(kk*acos(tau')));
        end
        v_interp = Tv*BETA((i-1)*N+1:i*N,:);
        % Chebyshev position matrix
        for kk = 0:N
            Tx(:,kk+1) = real(cos(kk*acos(tau')));
        end
        x_interp = Tx*ALPHA((i-1)*(N+1)+1:i*(N+1),:);
        
        % Store solution at user specified outputs
        Xout = [Xout; x_interp];
        Vout = [Vout; v_interp];
        Tout = [Tout; tt];
    end
    
end

return