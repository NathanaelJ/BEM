
% Remember to start a parallel pool
% README.md available on GitHub (github.com/NathanaelJ/BEM)

% Nathanael Jenkins, Usmaan Yaqoob
% Imperial College London, 2021

%% Initialise variables (user-defined)
clear all
tic

% Design Specifications
simName = 'Sample Test'; % Simulation name (appears in results & waitbar)
R = 0.2;       % Max turbine radius
hub = 0.025;    % Hub diameter
V1 = 10;        % On-design speed
V2 = 5;         % Off-design speed
lambdaMin = 2;  % Minimum TSR
lambdaMax = 8;  % Maximum TSR
dT = 0.25;      % Initial TSR change
Bmin = 2;       % Minimum blade number
Bmax = 3;       % Maximum blade number
fMin = 1.0;     % Minimum chord scale factor
fMax = 1.04;    % Maximum chord scale factor
dF = 0.02;      % Initial chord factor change

% Optimisation/ BEM control
elem = 30;      % Number of blade element sections (must be > 5)
delta = 0.25;   % Optimisation factor (effects convergence)
outerMax = 10;  % Maximum outer iterations
tol = 5*10^-4;  % Optimisation solution tolerance
dLambda = 0.2;  % Verificaion plot control
RCF = 4.2;      % Root chord limit factor (effects blade shape at root)
TCF = 0.9;      % Position along blade where tip reduction begins (0 < TCF < 1)

% Aerofoil geometry
foils.N = 2;    % Number of aerofoil sections
foils.plots = ["S833.dat", "S834.dat"]; % Aerofoil geometry file name(s)
foils.files = ["S833_5e4.dat", "S834_5e4.dat"]; % Polar file name(s)
foils.dist = [0, 0.15]; % Aerofoil section distribution 

% Design limitations
defMax = 0.02;  % Maximum blade tip deflection (m)
lRoot = 0.03;   % Maximum axial length at root (m)
wRoot = 0.05;   % Maximum tangential width at root (m)
tChord = 0.007; % Tip chord (m)
cMax = 0.05;    % Maximum chord (m)
omegaMax = 3000;% Maximum angular velocity (rpm), used for structural limits
Vmax = 1.8*10^-5;   % Maximum volume for a single blade (m^3)

% Material properties
sYield = 37*10^6; % Yield stress
E = 1.82*10^9;    % Youngs modulus
rho = 1225;       % Material density

%% Calculated variables (programmed)
A = pi*(R+hub)^2;   % Rotor area (including hub)
for i = 1:elem
    final.r(i) = ((R-hub)/(elem))*(i-0.5)+hub;
    % ^Radial position at the center of each blade element
end
lambdaMinInit = lambdaMin;
lambdaMaxInit = lambdaMax;
TSR_Struct = (R/V1)*omegaMax*2*pi/60;
% ^ Structural limit tip speed ratio

%% Import polars
% Error check
if length(foils.files) ~= foils.N
    error(['Incorrect number of foils specified. Expected ', num2str(foils.N)])
end

for i = 1:foils.N
    temp = string(foils.files(i));
    try
        foils.polars(:, :, i) = importdata(temp); % Import file
    catch
        disp(temp)
        error('Aerofoil polar file error! Check file (named above).')
    end
    
    % Error check
    temp = foils.polars(end, 1, i);
    temp2 = foils.polars(1, 1, i);
    if (temp2 > -45 || temp < 45)
        warning('POLARS NOT EXTRAPOLATED. BEM may diverge.')
    end
end

% Calculates ClDes and aoaDes based on first polar (Cl and aoa at L/D_max)
% Only considers L/D in the range -10<aoa<30, to avoid errors later on
i1 = find(foils.polars(:, 1, 1)>-10, 1, 'first');
i2 = find(foils.polars(:, 1, 1)<30, 1, 'last');
temp = foils.polars(:, 2, 1)./foils.polars(:, 3, 1);
[~, temp] = max(temp(i1:i2));
ClDes = foils.polars(temp+i1, 2, 1);
aoaDes = foils.polars(temp+i1, 1, 1);

%% OPTIMISATION
% Initialisation
outer = 0;
on = true;
w = waitbar(0, 'Starting optimisation...');
disp(' ')
disp(' ')
disp('OPTIMISING')
disp('------------------------------')
disp('Power    | CF      | TSR')

while on
    outer = outer+1;

    %% (For first outer iteration only) compares range of blade numbers
    for B = Bmin:Bmax
        i1 = 0;

        %% Compares range of chord distributions
        for factor = fMin:dF:fMax
            i1 = i1+1;
            i2 = 0;
            nc = (factor-fMin)/(fMax-fMin);

            %% Compares range of TSRs
            for TSR = lambdaMin:dT:lambdaMax
                i2 = i2+1;

                % Generate twist & chord distributions
                opt(B-Bmin+1).beta = deg2rad(twist('Betz', TSR, R, final.r, aoaDes));
                opt(B-Bmin+1).c = chord('Schmitz', TSR, R, final.r, B, ClDes, factor);

                % Limit root chord
                while ((opt(B-Bmin+1).c(1)*cos(opt(B-Bmin+1).beta(1)) > lRoot) || (opt(B-Bmin+1).c(1)*sin(opt(B-Bmin+1).beta(1)) > wRoot))
                    [~, index] = max(opt(B-Bmin+1).c);
                    if index == length(opt(B-Bmin+1).c)
                        index = index-1;
                    end

                    for j = 1:index
                        opt(B-Bmin+1).c(j) = opt(B-Bmin+1).c(j)*RCF*sqrt(final.r(j));
                    end
                end

                % Limit maximum chord
                while max(opt(B-Bmin+1).c) > cMax
                    [~, index] = max(opt(B-Bmin+1).c);
                    if index < 2
                        opt(B-Bmin+1).c(index) = 0.95*opt(B-Bmin+1).c(index);
                        index = 2;
                    end
                    opt(B-Bmin+1).c(index) = min(opt(B-Bmin+1).c(index-1), opt(B-Bmin+1).c(index+1));
                    f = opt(B-Bmin+1).c(index)/opt(B-Bmin+1).c(index+1);
                    opt(B-Bmin+1).c(index:end) = f*opt(B-Bmin+1).c(index:end);
                end

                % Limit tip chord
                while opt(B-Bmin+1).c(end) > tChord
                    for j = round(0.9*length(opt(B-Bmin+1).c)):length(opt(B-Bmin+1).c)
                        opt(B-Bmin+1).c(j) = opt(B-Bmin+1).c(j)*sqrt(1-(final.r(j)-final.r(round(0.9*length(opt(B-Bmin+1).c))))^2);
                    end
                end

                % Calculate BEM power output (operational)
                [opt(B-Bmin+1).Fa1, opt(B-Bmin+1).Ft1, opt(B-Bmin+1).P1(i2, i1)] = BEM(opt(B-Bmin+1).beta, foils.polars, B, opt(B-Bmin+1).c, R, TSR, length(opt(B-Bmin+1).beta), V1, foils.dist, final.r);

                % Calculate FEM performance (operational - deflection)
                [opt(B-Bmin+1).BS1, opt(B-Bmin+1).CS1, opt(B-Bmin+1).def1, opt(B-Bmin+1).V] = structR(foils.plots, foils.dist, opt(B-Bmin+1).Fa1, opt(B-Bmin+1).Ft1, final.r, V1, TSR, opt(B-Bmin+1).c, opt(B-Bmin+1).beta, sYield, E, rho, B);
                
                % Calculate BEM power output (maximum angular velocity)
                [FaStruct, FtStruct, ~] = BEM(opt(B-Bmin+1).beta, foils.polars, B, opt(B-Bmin+1).c, R, TSR_Struct, length(opt(B-Bmin+1).beta), V1, foils.dist, final.r);   % On-design conditions
                
                % Calculate FEM performance (maximum angular velocity)
                [~, CSMax_Max, ~, ~] = structR(foils.plots, foils.dist, FaStruct, FtStruct, final.r, V1, TSR_Struct, opt(B-Bmin+1).c, opt(B-Bmin+1).beta, sYield, E, rho, B);  % Structure (on-design)

                % Tests case for failure
                if opt(B-Bmin+1).P1(i2, i1) > (0.5*1.225*(pi*(R)^2)*V1^3) || opt(B-Bmin+1).def1 > defMax || CSMax_Max > sYield || max(opt(B-Bmin+1).c) > 0.053 || opt(B-Bmin+1).V > Vmax
                    opt(B-Bmin+1).fail1(i2, i1) = true;   
                end
                
                % Update waitbar
                num = nc+dF/(fMax-fMin)*(TSR-lambdaMin)/(lambdaMax-lambdaMin);
                waitbar(num, w, [simName, ' | ', 'Pass ', num2str(outer), ' | ', num2str(B), ' Blades'])
            end
        end
    end

    %% Plot iteration results
    figure(outer)
    hold on
    if outer == 1
        for i = 1:(Bmax-Bmin)+1
            indexes = opt(i).fail1';
            [t, f] = meshgrid(lambdaMin:dT:lambdaMax, fMin:dF:fMax);
            z = opt(i).P1';
            z(indexes) = NaN;
            try
                surf(t, f, z)
            catch
                if isempty(find(opt(i).fail1 == 0))
                    close(w)
                    error('Blade fails in all tests!')
                end
            end
        end
    else
        indexes = opt(1).fail1';
        [t, f] = meshgrid(lambdaMin:dT:lambdaMax, fMin:dF:fMax);
        z = opt(1).P1';
        z(indexes) = NaN;
        try
            surf(t, f, z)
        catch
            if isempty(find(opt(i).fail1 == 0))
                close(w)
                error('Blade fails in all tests!')
            end
        end
    end
    view(3)
    hold off
   
    %% Check convergence / change for next iteration
    % Convergence checks
    if outer >= outerMax
        on = false;
        conv = false;
    elseif outer > 1 && max(max(max((opt(B-Bmin+1).P1))))-maxP2 < tol
        on = false;
        conv = true;
    % Changing for next iteration
    else
        % Optimise blade number on first iteration only
        if outer == 1
            B = Bmin;
            z1 = opt(1).P1;
            z1(opt(1).fail1) = NaN;
            maxP = max(max(max(z1)));
            for i=1:(Bmax-Bmin+1)
                z2 = opt(i).P1;
                z2(opt(i).fail1) = NaN;
                if max(max(max(z2))) > maxP
                    B = i+Bmin-1;
                    maxP = max(max(max(z2)));
                end
            end
            Bmin = B;
            Bmax = B;
        end
        
        % Refine TSR and chord factor
        z = opt(B-Bmin+1).P1;
        z(opt(B-Bmin+1).fail1) = NaN;
        [maxP2, cfM] = max(max(z));
        cfMaxP = fMin + (cfM-1)*dF;
        [i2M, ~] = find(z==maxP2);
        TSRMaxP = lambdaMin+(i2M-1)*dT;
        
        % New TSR range
        L = (lambdaMax-lambdaMin)/dT;
        dT = delta*dT;
        lambdaMax = TSRMaxP+dT*L*0.55;
        lambdaMin = TSRMaxP-dT*L*0.55;
        
        % New chord factor range
        L = (fMax-fMin)/dF;
        dF = dF*delta;
        fMin = cfMaxP-dF*L*0.55;
        fMax = cfMaxP+dF*L*0.55;
    end

    %% Error check
    if isempty(lambdaMin) || sum(opt(B-Bmin+1).P1(opt(B-Bmin+1).fail1)) == 0
        close (w)
        error('Blade fails in all tests!')
    end
    disp([num2str(maxP2, '%6.4f'), '   ', num2str(cfMaxP, '%6.4f'), '     ', num2str(TSRMaxP, '%6.4f')])

end
close(w)
disp('OPTIMISATION COMPLETE | PROCESSING RESULTS')

%% Optimum design calculations
% Optimum TSR & chord factor
z = opt(B-Bmin+1).P1;
z(opt(B-Bmin+1).fail1) = NaN;
[maxP, cfM] = max(max(z));
factor = fMin + (cfM-1)*dF;
[i2M, ~] = find(z==maxP);
final.lambda = lambdaMin+(i2M-1)*dT;

%% Final BEM/ FEM script
 % Generate twist & chord distributions
final.beta = deg2rad(twist('Betz', final.lambda, R, final.r, aoaDes));
final.c = chord('Schmitz', final.lambda, R, final.r, B, ClDes, factor);

% Limit root chord
while ((opt(B-Bmin+1).c(1)*cos(opt(B-Bmin+1).beta(1)) > lRoot) || (opt(B-Bmin+1).c(1)*sin(opt(B-Bmin+1).beta(1)) > wRoot))
    [~, index] = max(opt(B-Bmin+1).c);
    if index == length(opt(B-Bmin+1).c)
        index = index-1;
    end
    for j = 1:index
        opt(B-Bmin+1).c(j) = opt(B-Bmin+1).c(j)*RCF*sqrt(final.r(j));
    end
end

% Limit tip chord
while opt(B-Bmin+1).c(end) > tChord
    for j = round(TCF*length(opt(B-Bmin+1).c)):length(opt(B-Bmin+1).c)
        opt(B-Bmin+1).c(j) = opt(B-Bmin+1).c(j)*sqrt(1-(final.r(j)-final.r(round(TCF*length(opt(B-Bmin+1).c))))^2);
    end
end

% % Limit maximum chord (disabled)
% while max(opt(B-Bmin+1).c) > 0.053
%     opt(B-Bmin+1).c = opt(B-Bmin+1).c*factor./(opt(B-Bmin+1).c(1)/wRoot);
% end

[final.Fa1, final.Ft1, final.P1] = BEM(final.beta, foils.polars, B, final.c, R, final.lambda, length(final.beta), V1, foils.dist, final.r);   % On-design conditions
[final.Fa2, final.Ft2, final.P2] = BEM(final.beta, foils.polars, B, final.c, R, final.lambda, length(final.beta), V2, foils.dist, final.r);   % Off-design conditions

[final.BS1, final.CS1, final.def1, ~] = structR(foils.plots, foils.dist, final.Fa1, final.Ft1, final.r, V1, final.lambda, final.c, final.beta, sYield, E, rho, B, true);  % Structure (on-design)
[final.BS2, final.CS2, final.def2, ~] = structR(foils.plots, foils.dist, final.Fa2, final.Ft2, final.r, V2, final.lambda, final.c, final.beta, sYield, E, rho, B);        % Structure (off-design)
% ^ Note true (line 271) generates structural plot for on-design conditions

[FaStruct, FtStruct, ~] = BEM(final.beta, foils.polars, B, final.c, R, TSR_Struct, length(final.beta), V1, foils.dist, final.r);   % On-design conditions (maximum structural load)
[FaStruct2, FtStruct2, ~] = BEM(final.beta, foils.polars, B, final.c, R, TSR_Struct, length(final.beta), V2, foils.dist, final.r); % Off-design conditions (maximum structural load)

[~, CSMax_Max, ~, ~] = structR(foils.plots, foils.dist, FaStruct, FtStruct, final.r, V1, TSR_Struct, final.c, final.beta, sYield, E, rho, B);       % Structure (on-design maximum)
[~, CSMax_Max2, ~, ~] = structR(foils.plots, foils.dist, FaStruct2, FtStruct2, final.r, V2, TSR_Struct, final.c, final.beta, sYield, E, rho, B);    % Structure (off-design maximum)

% Tests final cases for failure
if final.P1 > (0.5*1.225*(pi*(R)^2)*V1^3) || final.def1 > defMax || CSMax_Max > sYield || max(final.c) > 0.053
    final.fail1 = true;
    warning('Structural failure in on-design solution')
end
if final.P2 > (0.5*1.225*(pi*(R)^2)*V2^3) || final.def2 > defMax || CSMax_Max2 > sYield
    final.fail2 = true;
    warning('Structural failure in off-design solution')
end

%% Lambda verification plot
% Initialisation
disp(' ')
disp(' ')
disp('CALCULATING LAMBDA')
disp('---------------------------------------')
disp('ITERATION | POWER    | DIFF      | TSR')
lam.lambda(1) = lambdaMinInit;
lam.lambda(2) = lambdaMinInit + dLambda;
lam.beta = final.beta;
lam.c = final.c;

% Generate validation plot (twist & chord distribution fixed)
for iter = 1:(lambdaMaxInit-lambdaMinInit)/dLambda
    %% Calculate BEM power output
    [lam.Fa1, lam.Ft1, lam.P1(iter)] = BEM(lam.beta, foils.polars, B, lam.c, R, lam.lambda(iter), length(lam.beta), V1, foils.dist, final.r);   % On-design conditions
    [lam.Fa2, lam.Ft2, lam.P2(iter)] = BEM(lam.beta, foils.polars, B, lam.c, R, lam.lambda(iter), length(lam.beta), V2, foils.dist, final.r);   % Off-design conditions
    
    %% Prepare for another iteration
    lam.lambda(iter+1) = lam.lambda(iter)+dLambda;

    if iter>1
        disp([num2str(iter, '%6.2u'), '          ', num2str(lam.P1(iter), '%+6.4f'), '   ', num2str(lam.P1(iter) - lam.P1(iter-1), '%+6.6f'), '  ', num2str(lam.lambda(iter), '%6.6f')])
    end
end

% Plot final design
figure
hold on
title('Final blade design')
plot(final.lambda, final.P1, 'ro')
plot(final.lambda, final.P2, 'co')
plot(lam.lambda(1:(length(lam.P1))), lam.P1(1:(length(lam.P1))), 'b')
plot(lam.lambda(1:(length(lam.P2))), lam.P2(1:(length(lam.P2))), 'g')
legend('On-Design optimisation', 'Off-design Optimisation', 'On-Design Curve', 'Off-Design Curve', 'Pmax', 'location', 'best')
ylim([0 inf])
set(gca, 'XAxisLocation', 'origin')
xlabel('TSR')
ylabel('Power out (W)')
title(['Optimised power output = ', num2str(final.P1(end)), 'W at TSR = ', num2str(final.lambda(end))])
hold off

%% Results overview
disp(' ')
disp(' ')
disp('-------------------------------------')
disp('HAWT OPTIMISATION')
disp(simName)
disp(['TIME: ', num2str(toc)])
if conv
    disp('CONVERGED')
else
    disp('FAILED TO CONVERGE')
end
disp('-------------------------------------')
disp('AERODYNAMIC OPTIMISATION')
disp(['OPTIMISED POWER:       ', num2str(maxP)])
disp(['OPTIMISED TSR:         ', num2str(final.lambda)])
disp('-------------------------------------')
disp(['ROOT CHORD:            ', num2str(final.c(1))])
disp(['TIP CHORD:             ', num2str(final.c(end))])
disp(['TIP DEFLECTION:        ', num2str(final.def1)])
disp('-------------------------------------')

%% Ask user to continue to blade generation
cont = questdlg('Inspect plot. Continue or stop?', 'Computation paused', 'Generate blade', 'Stop', 'Generate blade');

%% Export CAD Parts
% ! NOTE - this method may fail to work properly when using multiple aerofoils to define the blade.
if cont == "Generate blade"
    disp(' ')
    disp(' ')
    disp('LOFTING')
    disp('---------------------------------------')
    loft2(final.c, -rad2deg(final.beta), foils.plots, foils.dist, R, hub)
end
