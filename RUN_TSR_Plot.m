
% Remember to start a parallel pool BEFORE running
% README.md available on GitHub (github.com/NathanaelJ/BEM)

% Nathanael Jenkins, Usmaan Yaqoob
% Imperial College London, 2021

%% Initialise variables (user-defined)
clear all
tic

% Design Specifications
simName = 'Sample Test'; % Simulation name (appears in results & waitbar)
R = 0.25;       % Max turbine radius
hub = 0.025;    % Hub diameter
V1 = 10;        % On-design speed
V2 = 5;         % Off-design speed
lambdaMin = 3;  % Min TSR
lambda = 6.2;   % Actual TSR
dLambda = 0.1;  % Change in TSR per step
lambdaMax = 10;  % Max TSR
lDmin = 4;      % Min design TSR
dT = 0.5;         % Change in design TSR per step
lDmax = 8;      % Max design TSR
B = 2;          % Number of blades
cFactor = 1.0;  % Chord scale factor

% BEM control
elem = 30;      % Number of blade element sections (must be > 5)
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

%% Calculated variables (programmed)
A = pi*(R+hub)^2;   % Rotor area (including hub)
for i = 1:elem
    result.r(i) = ((R-hub)/(elem))*(i-0.5)+hub;
    % ^Radial position at the center of each blade element
end
lamLim = (lambdaMax-lambdaMin)/dLambda;

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

%% Generate data

figure(1)
hold on
index = 0;
for lDes = lDmin:dT:lDmax
    index = index + 1;
    disp(' ')
    disp(' ')
    disp('CALCULATING LAMBDA')
    disp('---------------------------------------')
    disp('ITERATION | POWER    | DIFF     | TSR')
    opt(index).lambda(1) = lambdaMin;
    opt(index).lambda(2) = lambdaMin + dLambda;
    opt(index).beta = deg2rad(twist('Schmitz', lDes, R, result.r, aoaDes));
    opt(index).c = chord('Schmitz', lDes, R, result.r, B, ClDes);
    tic

    % Generate validation plot (twist & chord distribution fixed)
    for iter = 1:lamLim
        
        %% Calculate BEM power output (on-design only, by default)
        [opt(index).Fa1, opt(index).Ft1, opt(index).P1(iter)] = BEM(opt(index).beta, foils.polars, B, opt(index).c, R, opt(index).lambda(iter), length(opt(index).beta), V1, foils.dist, result.r);   % On-design conditions
        % [opt(index).Fa2, opt(index).Ft2, opt(index).P2(iter)] = BEM(opt(index).beta, foils.polars, B, opt(index).c, R, opt(index).lambda(iter), length(opt(index).beta), V2, foils.dist, result.r);   % Off-design conditions

        %% Prepare for another iteration
        opt(index).lambda(iter+1) = opt(index).lambda(iter)+dLambda;
        if iter>1
            disp([num2str(iter, '%6.2u'), '          ', num2str(opt(index).P1(iter), '%+6.4f'), '   ', num2str(opt(index).P1(iter) - opt(index).P1(iter-1), '%+6.6f'), '  ', num2str(opt(index).lambda(iter), '%6.6f')])
        end
        
        if opt(index).P1(iter) < -50
            disp('--Maximum power input. New pass started.--')
            break
        end
    end

    % Plot results (on-design only, by default)
    plot(opt(index).lambda(1:(length(opt(index).P1))), opt(index).P1(1:(length(opt(index).P1))))
    % plot(opt(index).lambda(1:(length(opt(index).P2))), opt(index).P2(1:(length(opt(index).P2))))
    drawnow
end

% Plot formatting
ylim([0 inf])
title('Effect of TSR and blade number on Power output')
hold off