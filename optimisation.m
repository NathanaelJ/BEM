
% Remember to start a parallel pool BEFORE running
% README.md available on GitHub (github.com/NathanaelJ/BEM)

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
lambdaMin = 2;  % Minimum TSR                                              ***
lambdaMax = 13; % Maximum TSR                                              ***
Bmin = 2;       % Min number of blades                                     ***
Bmax = 6;       % Max number of blades                                     ***
elem = 25;      % Number of blade element sections                         ***
tMax = 240;     % Max time for module 2 (optimisation) in seconds          ***
foils.N = 1;    % Number of aerofoil sections                              ***
foils.files = "QbladeAerofoil.dat"; % Polar file name(s)                   ***
foils.plots = "SG6043.txt"; % Aerofoil geometry file name(s)               ***

% Optimised variables (USER DEFINED)
lambda = 6.3;     % Initial TSR                                            ***
foils.dist = 0;   % Aerofoil section distribution                          ***
delta = 0.1;      % Change in lambda per step                              ***

%% Calculated variables (programmed)
A = pi*(R+hub)^2;   % Rotor area (including hub)
for i = 1:elem
    result.r(i) = ((R-hub)/(elem))*(i-0.5)+hub;
    % ^Radial position at the center of each blade element
end
result.lambda(1) = lambda;
result.lambda(2) = lambda+delta;
opt.lambda(1) = lambdaMin;
opt.lambda(2) = lambdaMin+delta;
lamLim = round((lambdaMax-lambdaMin)/delta);
deltaInit = delta;

%% Import polars
% Error check
if length(foils.files) ~= foils.N
    error(['Incorrect number of foils specified. Expected ', num2str(foils.N)])
end

w = waitbar(0, 'Importing polars');
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

%% Optimise for blade number
% Initialisation
waitbar(0.01, w, 'Optimising | Module 1 | Iter_I_n_n_e_r = 1 | B_m_a_x_P = 0')
disp(' ')
disp(' ')
disp('OPTIMISATION MODULE 1')
disp('---------------------------------------')
disp('ITERATION | POWER    | DIFF     | LAMBDA  | BLADES ')
Boptim = [0, 0, lambda];

% Figure formatting
figure(1)
hold on
xlabel('TSR')
ylabel('Power out (W)')
ylim([0 inf])
title('Optimising blade number')

% Optimisation
for B = Bmin:Bmax
    for iter = 1:lamLim
        %% Generate optimal twist & chord distributions
        opt.beta = deg2rad(twist('Schmitz', opt.lambda(iter), R, result.r, aoaDes));
        opt.c = chord('Schmitz', opt.lambda(iter), R, result.r, B, ClDes);
            
        % Error check
        if opt.beta(opt.beta > (pi/2)) ~= 0
            warning('Twist angles > 90º have been generated. Review polars')
        end
        if length((opt.beta)) ~= elem
            error(['Twist values differ from elements. Expected length ', num2str(elem)])
        elseif length((opt.c)) ~= elem
            error(['Chord values differ from elements. Expected length ', num2str(elem)])
        end

        %% Calculate BEM power output
        [opt.Fa1, opt.Ft1, opt.P1(iter)] = BEM(opt.beta, foils.polars, B, opt.c, R, opt.lambda(iter), length(opt.beta), V1, foils.dist, result.r);   % On-design conditions
        [opt.Fa2, opt.Ft2, opt.P2(iter)] = BEM(opt.beta, foils.polars, B, opt.c, R, opt.lambda(iter), length(opt.beta), V2, foils.dist, result.r);   % Off-design conditions

        %% Calculate FEM performance
        % (WAITING)
        % [] = structR()    % On-design conditions
        % [] = structR()    % Off-design conditions

        %% Prepare for another iteration
        opt.lambda(iter+1) = opt.lambda(iter)+delta;
        progress = 0.01+(0.69/(Bmax-Bmin))*(B-Bmin)/(Bmax-Bmin)+(0.69/(Bmax-Bmin))*iter/lamLim;
        waitbar(progress, w, ['Optimising | Module 1 | Iter_I_n_n_e_r = ', num2str(iter), ' | B_m_a_x_P = ', num2str(Boptim(1))])
        if iter > 1
            disp([num2str(iter, '%6.2u'), '          ', num2str(opt.P1(iter), '%+6.4f'), '   ', num2str(opt.P1(iter) - opt.P1(iter-1), '%+6.6f'), '  ', num2str(opt.lambda(iter), '%6.6f'), '  ', num2str(B, '%6.2u')])
        end
    end
    
    % Compare performance for each blade number
    [temp1, temp2] = max(opt.P1);
    if temp1 > Boptim(2)
        Boptim = [B, temp1, opt.lambda(temp2)];
    end
    
    % Update plot
    plot(opt.lambda(1:(length(opt.P1))), opt.P1(1:(length(opt.P1))))
    plot(opt.lambda(1:(length(opt.P2))), opt.P2(1:(length(opt.P2))))
end

% Plot formatting
str = {};
for i=Bmin:Bmax
    str(end+1) = {sprintf(num2str(i), 'On-design')};
    str(end+1) = {sprintf(num2str(i), 'Off-Design')};
end
legend(str, 'location', 'best')
title(['Optimum blade number: ', num2str(Boptim(1)), ' | P_m_a_x = ', num2str(Boptim(2))])
hold off

%% Optimise at selected blade number
% Initialisation
waitbar(0.01, w, 'Optimising | Module 2 | Iter = 1')
disp(' ')
disp(' ')
disp('OPTIMISTING')
disp('---------------------------------------')
disp('ITERATION | POWER    | DIFF     | LAMBDA')
on = true;
iter = 0;
result.lambda(1) = Boptim(3);   % Uses previously calculated value of optimum TSR
result.lambda(2) = result.lambda(1)+delta;
B = Boptim(1);
tic

% Optimisation
while on == true
    iter = iter+1;
    
    %% Generate optimal twist & chord distributions
    result.beta = deg2rad(twist('Schmitz', result.lambda(iter), R, result.r, aoaDes));
    result.c = chord('Schmitz', result.lambda(iter), R, result.r, B, ClDes);
    
    % Error check
    if result.beta(result.beta > (pi/2)) ~= 0
        warning('Twist angles > 90º have been generated. Review polars')
    end
    if length((result.beta)) ~= elem
        error(['Twist values differ from elements. Expected length ', num2str(elem)])
    elseif length((result.c)) ~= elem
        error(['Chord values differ from elements. Expected length ', num2str(elem)])
    end
    
    %% Calculate BEM power output
    [result.Fa1, result.Ft1, result.P1(iter)] = BEM(result.beta, foils.polars, B, result.c, R, result.lambda(iter), length(result.beta), V1, foils.dist, result.r);   % On-design conditions
    [result.Fa2, result.Ft2, result.P2(iter)] = BEM(result.beta, foils.polars, B, result.c, R, result.lambda(iter), length(result.beta), V2, foils.dist, result.r);   % Off-design conditions
    
    %% Calculate FEM performance
    % (WAITING)
    % [] = structR()    % On-design conditions
    % [] = structR()    % Off-design conditions
    
    %% Prepare for another iteration
    % Adjusts lambda depending on power change (optimisation procedure)
    if iter > 1
        disp([num2str(iter, '%6.2u'), '          ', num2str(opt.P1(iter), '%+6.4f'), '   ', num2str(opt.P1(iter) - opt.P1(iter-1), '%+6.6f'), '  ', num2str(opt.lambda(iter), '%6.6f')])
        if result.P1(iter) < result.P1(iter-1)
            if result.lambda(iter) < result.lambda(iter-1)
                result.lambda(iter+1) = result.lambda(iter)+delta;
            else
                result.lambda(iter+1) = result.lambda(iter)-delta;
            end
        else
            if result.lambda(iter) < result.lambda(iter-1)
                result.lambda(iter+1) = result.lambda(iter)-delta;
            else
                result.lambda(iter+1) = result.lambda(iter)+delta;
            end
        end
    
        if abs(result.P1(iter) - result.P1(iter-1)) < 10^-6
            on = false;
            disp('----CONVERGED----')
        elseif toc > tMax
            on = false;
            disp('Maximum time exceeded. Optimisation may not have converged')
        elseif abs(result.P1(iter)-result.P1(iter-1)) < delta*100
            delta = 0.7*delta*(0.99+0.02*rand);
        end
        % Randomisation to prevent oscillation
        delta = delta*(0.99+0.02*rand);
        
        % Update waitbar (assumes ~ 20 iterations for convergence)
        if iter <= 20
            num = 0.7+0.1*(iter/10)^0.4;
            waitbar(num, w, ['Optimising | Module 2 | Iter = ', num2str(iter)])
        else
            waitbar(0.8, w, ['Optimising | Module 2 | Iter = ', num2str(iter)])
        end
    end
end
% Plot result
figure(2)
hold on
title('Optimising lambda')
plot(result.lambda(1:length(result.P1)), result.P1(1:length(result.P1)), 'ro')
plot(result.lambda(1:length(result.P2)), result.P2(1:length(result.P2)), 'co')
xlabel('TSR')
ylabel('Power out (W)')
title(['Optimised power output = ', num2str(result.P1(end)), 'W at TSR = ', num2str(result.lambda(end))])
hold off

%% Generate verification plot
% Initialisation
disp(' ')
disp(' ')
waitbar(0.41, w, 'Optimising | Module 2 | Calculating Lambda')
disp('CALCULATING LAMBDA')
disp('---------------------------------------')
disp('ITERATION | POWER    | DIFF     | LAMBDA')
delta = 0.5*deltaInit;          % Reset delta
opt.lambda(1) = lambdaMin;
opt.lambda(2) = lambdaMin + delta;
opt.beta = result.beta;
opt.c = result.c;
tic

% Generate validation plot (twist & chord distribution fixed)
for iter = 1:lamLim
    %% Calculate BEM power output
    [opt.Fa1, opt.Ft1, opt.P1(iter)] = BEM(opt.beta, foils.polars, B, opt.c, R, opt.lambda(iter), length(opt.beta), V1, foils.dist, result.r);   % On-design conditions
    [opt.Fa2, opt.Ft2, opt.P2(iter)] = BEM(opt.beta, foils.polars, B, opt.c, R, opt.lambda(iter), length(opt.beta), V2, foils.dist, result.r);   % Off-design conditions

    %% Calculate FEM performance
    % (WAITING)
    % [] = structR()    % On-design conditions
    % [] = structR()    % Off-design conditions
    
    %% Prepare for another iteration
    opt.lambda(iter+1) = opt.lambda(iter)+delta;
    num = 0.8+0.2*iter/lamLim;
    waitbar(num, w, 'Optimising | Module 2 | Calculating Lambda')
    if iter>1
        disp([num2str(iter, '%6.2u'), '          ', num2str(opt.P1(iter), '%+6.4f'), '   ', num2str(opt.P1(iter) - opt.P1(iter-1), '%+6.6f'), '  ', num2str(opt.lambda(iter), '%6.6f')])
    end
end
close(w)

% Plot results
figure(2)
hold on
plot(opt.lambda(1:(length(opt.P1))), opt.P1(1:(length(opt.P1))), 'b')
plot(opt.lambda(1:(length(opt.P2))), opt.P2(1:(length(opt.P2))), 'g')
plot([result.lambda(end) result.lambda(end)], [0 result.P1(end)], 'k-')
ylim([0 inf])
set(gca, 'XAxisLocation', 'origin')
legend('On-Design optimisation', 'Off-design Optimisation', 'On-Design Curve', 'Off-Design Curve', 'Pmax', 'location', 'best')
hold off

%% OPTIONAL User-defined optimisation
disp(' ')
disp(' ')
cont = questdlg('Run user-defined optimisation?', 'Optimisation', 'Yes', 'No', 'No');

if cont == "Yes"
    disp(' ')
    disp(' ')
    disp('OPTIMISTING | MODULE 3 (USER-LED)')
    disp('---------------------------------------')
    disp('ITERATION | POWER    | DIFF     | LAMBDA')
    figure(3)
    xlabel('TSR')
    ylabel('Power out (W)')
    title(['Optimised power output = ', num2str(result.P1(end)), ' at TSR = ', num2str(result.lambda(end))])
    hold on
    l1 = animatedline;
    l2 = animatedline;
    on = true;
    letrun = false;
    iter = 0;
    tic

    % Optimisation
    while on == true
        iter = iter+1;

        %% Generate optimal twist & chord distributions
        result.beta = deg2rad(twist('Schmitz', result.lambda(iter), R, result.r, aoaDes));
        result.c = chord('Schmitz', result.lambda(iter), R, result.r, B, ClDes);

        % Error check
        if result.beta(result.beta > (pi/2)) ~= 0
            warning('Twist angles > 90º have been generated. Review polars')
        end
        if length((result.beta)) ~= elem
            error(['Twist values differ from elements. Expected length ', num2str(elem)])
        elseif length((result.c)) ~= elem
            error(['Chord values differ from elements. Expected length ', num2str(elem)])
        end

        %% Calculate BEM power output
        [result.Fa1, result.Ft1, result.P1(iter)] = BEM(result.beta, foils.polars, B, result.c, R, result.lambda(iter), length(result.beta), V1, foils.dist, result.r);   % On-design conditions
        [result.Fa2, result.Ft2, result.P2(iter)] = BEM(result.beta, foils.polars, B, result.c, R, result.lambda(iter), length(result.beta), V2, foils.dist, result.r);   % Off-design conditions

        %% Calculate FEM performance
        % (WAITING)
        % [] = structR()    % On-design conditions
        % [] = structR()    % Off-design conditions

        %% Prepare for another iteration
        % Plot
        addpoints(l1, result.lambda(end), result.P1(end));
        addpoints(l2, result.lambda(end), result.P2(end));
        drawnow update
        
        % Adjust optimisation variables
        if iter > 1
            disp([num2str(iter, '%6.2u'), '          ', num2str(opt.P1(iter), '%+6.4f'), '   ', num2str(opt.P1(iter) - opt.P1(iter-1), '%+6.6f'), '  ', num2str(opt.lambda(iter), '%6.6f')])
        end
        if letrun == true
            if result.P1(iter) < result.P1(iter-1)
                if result.lambda(iter) < result.lambda(iter-1)
                    result.lambda(iter+1) = result.lambda(iter)+delta;
                else
                    result.lambda(iter+1) = result.lambda(iter)-delta;
                end
            else
                if result.lambda(iter) < result.lambda(iter-1)
                    result.lambda(iter+1) = result.lambda(iter)-delta;
                else
                    result.lambda(iter+1) = result.lambda(iter)+delta;
                end
            end

            if abs(result.P1(iter) - result.P1(iter-1)) < 10^-6
                on = false;
                disp('----CONVERGED----')
            elseif toc > tMax
                on = false;
                disp('Maximum time exceeded. Optimisation may not have converged')
            elseif abs(result.P1(iter)-result.P1(iter-1)) < delta*50
                delta = 0.7*delta*(0.99+0.02*rand);
            end
            % Randomisation to prevent oscillation
            delta = delta*(0.99+0.02*rand);

        else
            change = listdlg('ListString', {'B', 'lambda', 'aoaDes', 'none'}, 'PromptString', ['Change variables? | B = ', num2str(B), ' | lambda = ', num2str(result.lambda(end)), ' | aoaDes = ', num2str(aoaDes)], 'SelectionMode', 'multiple', 'Name', 'Optimisation', 'CancelString', 'End user-defined sequence', 'ListSize', [300 70]);
            if sum(ismember(change, 1)) == 1
                new = inputdlg('Enter a new number of blades (INT)', 'Optimisation module 3', 1, {num2str(B)});
                B = cell2mat(new);
            end
            if sum(ismember(change, 2)) == 1
                new = inputdlg('Enter a new TSR', 'Optimisation module 3', 1, {num2str(result.lambda(iter))});
                result.lambda(iter+1) = str2double(new);
            end
            if sum(ismember(change, 3)) == 1
                new = inputdlg('Enter a new aoaDes', 'Optimisation module 3', 1, {num2str(aoaDes)});
                aoaDes = str2double(new);
                index = find(foils.polars(:, 1, 1)>aoaDes, 1, 'first');
                ClDes = foils.polars(index, 2, 1);
            end
            if isempty(change)
                letrun = true;
                tic
            end
        end
    end
    hold off
end

disp(' ')
disp(' ')
disp('---------------------------')
disp('OPTIMISATION COMPLETE')

%% Ask user to continue to blade generation
disp(' ')
disp(' ')
cont = questdlg('Inspect plot. Continue or stop?', 'Computation paused', 'Generate blade', 'Stop', 'Generate blade');

%% Export CAD Parts
if cont == "Generate blade"
    disp(' ')
    disp(' ')
    disp('LOFTING')
    disp('---------------------------------------')
    loft2(result.c, rad2deg(result.beta), foils.plots, foils.dist, R, hub)
end