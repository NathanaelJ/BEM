
% Remember to start a parallel pool BEFORE running
% README.md available on GitHub (LINK)

% Nathanael Jenkins, Usmaan Yaqoob
% Imperial College London, 2021

%% Initialise variables (user-defined)
clear all

% Fixed variables (USER DEFINED)
rho = 1.225;    % Air density
R = 0.25;       % Max radius
hub = 0.025;    % Hub diameter
V1 = 10;        % On-design speed
V2 = 5;         % Off-design speed
lambdaMin = 1.5;  % Min TSR                                                  ***
lambdaMax = 14; % Max TSR                                                  ***
B = 2;          % Number of blades                                         ***
elem = 25;      % Number of blade element sections                         ***
foils.N = 1;    % Number of aerofoil sections                              ***
foils.files = "QbladeAerofoil.dat"; % Polar file name(s)                   ***
foils.plots = "SG6043.txt"; % Aerofoil geometry file name(s)               ***
foils.dist = 0; % Aerofoil section distribution                            ***
delta = 0.1;    % Change in lambda per step                                ***

%% Calculated variables (programmed)
A = pi*(R+hub)^2;   % Rotor area (including hub)
for i = 1:elem
    result.r(i) = ((R-hub)/(elem))*(i-0.5)+hub;
    % ^Radial position at the center of each blade element
end
opt.lambda(1) = lambdaMin;
opt.lambda(2) = lambdaMin+delta;
lamLim = round((lambdaMax-lambdaMin)/delta);

%% Import polars
% Error check
if length(foils.files) ~= foils.N
    error(['Incorrect number of foils specified. Expected ', num2str(foils.N)])
end

for i = 1:foils.N
    temp = string(foils.files(i));
    foils.polars(:, :, i) = importdata(temp); % Import file
    
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

%% Generate varying TSR plot
% Initialisation

figure(1)
hold on
% lDes MULTIPLIED BY TWO FOR INDICES
for lDes = 16:32
    disp(' ')
    disp(' ')
    disp('CALCULATING LAMBDA')
    disp('---------------------------------------')
    disp('ITERATION | POWER    | DIFF     | LAMBDA')
    opt(lDes).lambda(1) = lambdaMin;
    opt(lDes).lambda(2) = lambdaMin + delta;
    opt(lDes).beta = deg2rad(twist('Schmitz', lDes/4, R, result.r, aoaDes));
    opt(lDes).c = chord('Schmitz', lDes/4, R, result.r, B, ClDes);
    tic

    % Generate validation plot (twist & chord distribution fixed)
    for iter = 1:lamLim
        
        %% Calculate BEM power output
        [opt(lDes).Fa1, opt(lDes).Ft1, opt(lDes).P1(iter)] = BEM(opt(lDes).beta, foils.polars, B, opt(lDes).c, R, opt(lDes).lambda(iter), length(opt(lDes).beta), V1, foils.dist, result.r);   % On-design conditions
        % [opt(lDes).Fa2, opt(lDes).Ft2, opt(lDes).P2(iter)] = BEM(opt(lDes).beta, foils.polars, B, opt(lDes).c, R, opt(lDes).lambda(iter), length(opt(lDes).beta), V2, foils.dist, result.r);   % Off-design conditions

        %% Calculate FEM performance
        % (WAITING)
        % [] = structR()    % On-design conditions
        % [] = structR()    % Off-design conditions

        %% Prepare for another iteration
        opt(lDes).lambda(iter+1) = opt(lDes).lambda(iter)+delta;
        if iter>1
            disp([num2str(iter, '%6.2u'), '          ', num2str(opt(lDes).P1(iter), '%+6.4f'), '   ', num2str(opt(lDes).P1(iter) - opt(lDes).P1(iter-1), '%+6.6f'), '  ', num2str(opt(lDes).lambda(iter), '%6.6f')])
        end
        
        if opt(lDes).P1(iter) < -50
            disp('--Power input max--')
            break
        end
    end

    % Plot results
    plot(opt(lDes).lambda(1:(length(opt(lDes).P1))), opt(lDes).P1(1:(length(opt(lDes).P1))))
    % plot(opt(lDes).lambda(1:(length(opt(lDes).P2))), opt(lDes).P2(1:(length(opt(lDes).P2))))
    drawnow
end

% Plot formatting
% str = {};
% for i=8:0.5:8
%     str(end+1) = {sprintf(num2str(i), 'On-design')};
%     % str(end+1) = {sprintf(num2str(i), 'Off-Design')};
% end
ylim([0 inf])
% legend(str, 'location', 'best')
legend('4', '4.25', '4.5', '4.75', '5', '5.25', '5.5', '5.75', '6', '6.25', '6.5', '6.75', '7', '7.25', '7.5', '7.75', '8')
title('Effect of TSR and blade number on Power output')
hold off