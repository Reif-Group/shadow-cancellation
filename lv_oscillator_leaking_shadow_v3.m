%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Rajiv Nagipogu
% Affiliation: Duke University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;
concX = 4; concY = 2; stopTime = 3600; nCatalytic = 3;

% compare ideal lotka-volterra simulation with approximation
% ideal_lv_oscillator(concX, concY, stopTime, 1);
% ideal_cat_lv_oscillator(concX, concY, nCatalytic, stopTime, 1);
% approx_cat_lv_oscillator(concX, concY, nCatalytic, stopTime, 1);
% approx_link_cat_lv_oscillator(concX, concY, nCatalytic, stopTime, 1);
approx_link_cat_lv_oscillator_leaky_shadow(concX, concY, nCatalytic, stopTime, 1);

function ideal_lv_oscillator(concA, concB, stopTime, figNo)
    % enter A and B concentration assuming all gates are excess
    %
    % stopTime is the total simulation time
    % figNo is the display figure number
    
    % rate of reactions 
    rate = 1/60;
    % allowed error rate tolerance value in pico molars
    absTol = 1e-12;
    relTol = 1e-12;
    
    % Create the amplification model
    model = sbiomodel('Ideal CRNs for LV oscillator');

    r = cell(1, 3);
    k = cell(1, length(r));
    p = cell(1, length(r));
         
    % Add reaction set X -> 2X to the solver
    r{1} = addreaction(model,'X -> X + X');
    r{2} = addreaction(model,'Y -> W');
    r{3} = addreaction(model,'X + Y -> Y + Y');
        
    % Set the Kinetic Law for Reactions.
    for i = 1:length(r)
        k{i} = addkineticlaw(r{i}, 'MassAction'); 
    end
    
    p{1} = addparameter(k{1}, 'c1', 'Value', rate);
    p{2} = addparameter(k{2}, 'c2', 'Value', rate);
    p{3} = addparameter(k{3}, 'c3', 'Value', rate);
    
    % Set the Kinetic Law for Reactions.
    for i = 1:length(r)
        k{i}.ParameterVariableNames = {strcat('c',num2str(i))};
    end
    
    % Set initial amounts for species 
    r{1}.Reactants(1).InitialAmount = concA;        % A
    r{2}.Reactants(1).InitialAmount = concB;        % B

    % Display the Completed Model Objects
    model

    % Display the Reaction Objects
    model.Reactions

    % Display the Species Objects
    model.Species

    % Simulate ODE and plot for totalTime
    cs = getconfigset(model,'active');
    set(cs.SolverOptions, 'AbsoluteTolerance', absTol);
    set(cs.SolverOptions, 'RelativeTolerance', relTol);
    cs.SolverType = 'ode15s';
    cs.StopTime = stopTime;
    [t,X] = sbiosimulate(model);

    % plot X over the time
    figure(figNo); subplot(2, 1, 1);
    hold on; box on;
    plot(X(:,1), X(:,2), 'LineWidth', 2.0);
    ylabel('Y'); xlabel('X');
    set(gca, 'LineWidth', 2.0); 
    
    subplot(2, 1, 2);
    box on; hold on;
    plot(t/60.0, X(:,1), 'LineWidth', 2.0); % CHECK: changed the time scale
    plot(t/60.0, X(:,2), 'LineWidth', 2.0); % CHECK: changed the time scale
    ylabel('Concentration'); xlabel('Time');
    set(gca, 'LineWidth', 2.0); 
    legend('X', 'Y')
end

function ideal_cat_lv_oscillator(initX, initY, N, stopTime, figNo)
    % implements catalytic version of an ideal lotka volterra oscillator
    % reaction
    %
    % initX is the initial concentration of X
    % initY is the initial concentration of Y
    % N is the number of catalytic reactions
    % stopTime is simulation time (eg. 100)
    % figNo is the figure # for plot

    % set all the default params and initialize model
    % default rate and concentration
    rate = 1/60;
    base = 1;
    
    % allowed error rate tolerance value in pico molars
    absTol = 1e-12;
    relTol = 1e-12;
    
    % create a sym bio model for simulation
    model = sbiomodel('Approximate CRNs for LV oscillator');
    r = cell(1, 3*N);
    k = cell(1, length(r));
    p = cell(1, length(r));
    
    % add catalytic approximations for X -> X + X step
    for i = 1 : N - 1
        r{i} = addreaction(model, strcat('x', int2str(i), ' -> x', int2str(i),' + x', int2str(i+1))); 
    end
    r{N} = addreaction(model, strcat('x', int2str(N), ' -> x', int2str(N), + ' + x1')); 
    
    % add catalytic approximations for X + Y -> Y + Y step
    for i = N + 1 : 2 * N - 1
        r{i} = addreaction(model, strcat('x', int2str(i-N), ' + y', int2str(i-N),...
            ' -> y', int2str(i-N),' + y', int2str(i-N+1))); 
    end
    r{2*N} = addreaction(model, strcat('x', int2str(N), ' + y', int2str(N),...
            ' -> y', int2str(N),' + y1')); 
        
    % add one ideal waste reaction for LV
    for i = 2 * N + 1 : 3 * N
        r{i} = addreaction(model,strcat('y',  int2str(i-2*N), ' -> waste')); 
    end
   
    % add all the rate laws as mass action kinetics
    for i = 1 : N
        k{i} = addkineticlaw(r{i}, 'MassAction'); 
        k{i}.ParameterVariableNames = {strcat('c', int2str(i))};
        p{i} = addparameter(k{i}, strcat('c', int2str(i)), 'Value', rate);
    end
    for i = 1+N : 2*N
        k{i} = addkineticlaw(r{i}, 'MassAction'); 
        k{i}.ParameterVariableNames = {strcat('c', int2str(i))};
        p{i} = addparameter(k{i}, strcat('c', int2str(i)), 'Value', N*rate);
    end
    for i = 2*N+1 : 3*N
        k{i} = addkineticlaw(r{i}, 'MassAction'); 
        k{i}.ParameterVariableNames = {strcat('c', int2str(i))};
        p{i} = addparameter(k{i}, strcat('c', int2str(i)), 'Value', rate);
    end

    % set initial concentrations
    for i = 1:N
        r{N+i}.Reactants(1).InitialAmount = initX * base * (1/N);
    end
    for i = 1:N
        r{N+i}.Reactants(2).InitialAmount = initY * base * (1/N);
    end
    
    
    model
    model.Reactions
    model.Species

    % get solver config and set ode type, tolerance, stop time
    cs = getconfigset(model,'active');
    set(cs.SolverOptions, 'AbsoluteTolerance', absTol);
    set(cs.SolverOptions, 'RelativeTolerance', relTol);
    cs.SolverType = 'ode15s';
    cs.StopTime = stopTime;
    [t,X] = sbiosimulate(model);
    
    % plot
    figure(figNo); subplot(2, 1, 1);
    hold on; box on;
    plot(sum(X(:, 1:N), 2), sum(X(:, N+1:2*N), 2), ':', 'LineWidth', 2.0);
    ylabel('Y'); xlabel('X');
    set(gca, 'LineWidth', 2.0); 
    
    subplot(2, 1, 2);
    box on; hold on;
    plot(t./60, sum(X(:, 1:N), 2), 'LineWidth', 2.0);
    plot(t./60, sum(X(:, N+1:2*N), 2), 'LineWidth', 2.0);
    ylabel('Concentration'); xlabel('Time');
    set(gca, 'LineWidth', 2.0); legend('X', 'Y'); %ylim([0 5])
end

function approx_cat_lv_oscillator(initX, initY, N, stopTime, figNo)
 % implements a DNA version of catalytic lotka volterra oscillator 
    %
    % initX is the initial concentration of X
    % initY is the initial concentration of Y
    % N is the number of catalytic reactions
    % stopTime is simulation time (eg. 100)
    % figNo is the figure # for plot

    % set all the default params and initialize model
    % default rate and concentration
    base = 1;
    % slow down all reactions by multiplying with a scaling factor 1e-3
    rate = (1/60);
    % fastest reaction are assumed to occur 
    infRate = 1*(rate);
    % inf concentration of gates will be 1000 uM
    infConc = 1e3*(base);
    
    % allowed error rate tolerance value in pico molars
    absTol = 1e-12;
    relTol = 1e-12;
    
    % create a sym bio model for simulation
    model = sbiomodel('Approximate CRNs for LV oscillator');
    r = cell(1, 5*N);
    k = cell(1, length(r));
    p = cell(1, length(r));
    
    % add two-step catalytic approximations for X -> X + X step
    % For N catalytic reaction, we have 2N reactions
    for i = 1 : N - 1
        r{2*i-1} = addreaction(model, strcat('x', int2str(i), ' + Gx', int2str(i), ...
            '1 -> Ix', int2str(i))); 
        r{2*i} = addreaction(model, strcat('Ix', int2str(i), ...
            ' + Gx', int2str(i),'2 -> x', int2str(i),' + x', int2str(i+1))); 
    end
    r{2*N-1} = addreaction(model, strcat('x', int2str(N), ' + Gx', int2str(N), ...
            '1 -> Ix', int2str(N)));
    r{2*N} = addreaction(model, strcat('Ix', int2str(N), ...
            ' + Gx', int2str(N),'2 -> x', int2str(N),' + x1')); 
    
    % add catalytic approximations for X + Y -> Y + Y step
    for i = 1 : N - 1
        r{2*N+2*i-1} = addreaction(model, strcat('x', int2str(i), ' + y', int2str(i),...
            ' -> Iy', int2str(i))); 
        r{2*N+2*i} = addreaction(model, strcat('Iy', int2str(i), ' + Gy', int2str(i),...
            '2 -> y', int2str(i),' + y', int2str(i+1))); 
    end
    r{4*N-1} = addreaction(model, strcat('x', int2str(N), ' + y', int2str(N),...
            ' -> Iy', int2str(N))); 
    r{4*N} = addreaction(model, strcat('Iy', int2str(N), ' + Gy', int2str(N),...
            '2 -> y', int2str(N),' + y1')); 
        
    % add one ideal waste reaction for LV, and linker reactions
    for i = 4 * N + 1 : 5 * N
        r{i} = addreaction(model,strcat('y',  int2str(i-4*N), ' + Gw', int2str(i-4*N),' -> waste')); 
    end
   
    for i = 1 : 5*N
        k{i} = addkineticlaw(r{i}, 'MassAction'); 
        k{i}.ParameterVariableNames = {strcat('c', int2str(i))};
    end
    
    
    q1 = rate/(infConc); q2 = infRate;
    q3 = rate*N; q4 = infRate;
    for i = 1 : N
        % add rates for X -> 2X
        p{2*i-1} = addparameter(k{2*i-1}, strcat('c', int2str(2*i-1)), 'Value', q1);
        p{2*i} = addparameter(k{2*i}, strcat('c', int2str(2*i)), 'Value', q2);
        % add rates for X + Y -> 2Y
        p{2*N+2*i-1} = addparameter(k{2*N+2*i-1}, strcat('c', int2str(2*N+2*i-1)), 'Value', q3);
        p{2*N+2*i} = addparameter(k{2*N+2*i}, strcat('c', int2str(2*N+2*i)), 'Value', q4);
        % add rates for Y -> phi
        p{4*N+i} = addparameter(k{4*N+i}, strcat('c', int2str(4*N+i)), 'Value', q1);
    end
    

    % set initial concentrations for x and y
    for i = 1 : N
        r{2*N+2*i-1}.Reactants(1).InitialAmount = initX * base * (1/N); % X
        r{2*N+2*i-1}.Reactants(2).InitialAmount = initY * base * (1/N); % Gy
    end
    % set initial concentrations for excess gates
    for i = 1 : N
        r{2*i-1}.Reactants(2).InitialAmount = infConc;
        r{2*i}.Reactants(2).InitialAmount = infConc;
        r{2*N+2*i}.Reactants(2).InitialAmount = infConc;
        r{4*N+i}.Reactants(2).InitialAmount = infConc;
    end
    
    model
    model.Reactions
    model.Species

    % get solver config and set ode type, tolerance, stop time
    cs = getconfigset(model,'active');
    set(cs.SolverOptions, 'AbsoluteTolerance', absTol);
    set(cs.SolverOptions, 'RelativeTolerance', relTol);
    cs.SolverType = 'ode15s';
    cs.StopTime = stopTime;
    [t,X] = sbiosimulate(model);
    
    % plot
    x =  sum(X(:, 1:4:N*4), 2); y = sum(X(:, 4*N+1:3:4*N+3*N), 2);
    figure(figNo); subplot(2, 1, 1);
    hold on; box on;
    plot(x, y, 'LineWidth', 2.0);
    ylabel('Y'); xlabel('X');
    set(gca, 'LineWidth', 2.0); 
    
    subplot(2, 1, 2);
    box on; hold on;
    plot(t./60, x, ':','LineWidth', 2.0);
    plot(t./60, y, ':','LineWidth', 2.0);
    ylabel('Concentration (nM)'); xlabel('Time (mins)');
    set(gca, 'LineWidth', 2.0);
    legend('cat. X', 'cat. Y', 'approx. X', 'approx. Y'); 
end

function approx_link_cat_lv_oscillator(initX, initY, N, stopTime, figNo)
    % implements a DNA version of catalytic lotka volterra oscillator 
    %
    % initX is the initial concentration of X
    % initY is the initial concentration of Y
    % N is the number of catalytic reactions
    % stopTime is simulation time (eg. 100)
    % figNo is the figure # for plot

    % set all the default params and initialize model
    % default rate and concentration
    base = 1;
    % slow down all reactions by multiplying with a scaling factor 1e-3
    rate = (1/60);
    % fastest reaction are assumed to occur 
    infRate = 1e5*(rate); % different from approx cat
    % inf concentration of gates will be 1000 uM
    infConc = 1e3*(base);
    
    % allowed error rate tolerance value in pico molars
    absTol = 1e-12;
    relTol = 1e-12;
    
    % create a sym bio model for simulation
    model = sbiomodel('Approximate CRNs for LV oscillator');
    r = cell(1, 7*N);
    k = cell(1, length(r));
    p = cell(1, length(r));
    
    % add two-step catalytic approximations for X -> X + X step
    % For N catalytic reaction, we have 2N reactions
    for i = 1 : N - 1
        r{2*i-1} = addreaction(model, strcat('x', int2str(i), ' + Gx', int2str(i), ...
            '1 -> Ix', int2str(i))); 
        r{2*i} = addreaction(model, strcat('Ix', int2str(i), ...
            ' + Gx', int2str(i),'2 -> x', int2str(i),' + x', int2str(i+1))); 
    end
    r{2*N-1} = addreaction(model, strcat('x', int2str(N), ' + Gx', int2str(N), ...
            '1 -> Ix', int2str(N)));
    r{2*N} = addreaction(model, strcat('Ix', int2str(N), ...
            ' + Gx', int2str(N),'2 -> x', int2str(N),' + x1')); 
    
    % add catalytic approximations for X + Y -> Y + Y step
    for i = 1 : N - 1
        r{2*N+2*i-1} = addreaction(model, strcat('x', int2str(i), ' + Gy', int2str(i),...
            '1 -> Iy', int2str(i))); 
        r{2*N+2*i} = addreaction(model, strcat('Iy', int2str(i), ' + Gy', int2str(i),...
            '2 -> y', int2str(i),' + y', int2str(i+1))); 
    end
    r{4*N-1} = addreaction(model, strcat('x', int2str(N), ' + Gy', int2str(N),...
            '1 -> Iy', int2str(N))); 
    r{4*N} = addreaction(model, strcat('Iy', int2str(N), ' + Gy', int2str(N),...
            '2 -> y', int2str(N),' + y1')); 
        
    % add one ideal waste reaction for LV, and linker reactions
    for i = 4 * N + 1 : 5 * N
        r{i} = addreaction(model,strcat('y',  int2str(i-4*N), ' + Gw', int2str(i-4*N),' -> waste')); 
        r{N+i} = addreaction(model,strcat('Gy',  int2str(i-4*N), '1 + Gw', int2str(i-4*N),' -> waste')); 
    end
   
    % add linker reactions
    for i = 6*N+1:7*N
        r{i} = addreaction(model,strcat('y',  int2str(i-6*N),' -> Gy', int2str(i-6*N), '1')); 
        r{i+N} = addreaction(model,strcat('Gy',  int2str(i-6*N),'1 -> y', int2str(i-6*N)));
    end
   
    for i = 1 : 8*N
        k{i} = addkineticlaw(r{i}, 'MassAction'); 
        k{i}.ParameterVariableNames = {strcat('c', int2str(i))};
    end
    
    q1 = rate/(infConc); q2 = infRate;
    q3 = rate*N; q4 = infRate;
    for i = 1 : N
        % add rates for X -> 2X
        p{2*i-1} = addparameter(k{2*i-1}, strcat('c', int2str(2*i-1)), 'Value', q1);
        p{2*i} = addparameter(k{2*i}, strcat('c', int2str(2*i)), 'Value', q2);
        % add rates for X + Y -> 2Y
        p{2*N+2*i-1} = addparameter(k{2*N+2*i-1}, strcat('c', int2str(2*N+2*i-1)), 'Value', q3);
        p{2*N+2*i} = addparameter(k{2*N+2*i}, strcat('c', int2str(2*N+2*i)), 'Value', q4);
        % add rates for Y -> phi
        p{4*N+i} = addparameter(k{4*N+i}, strcat('c', int2str(4*N+i)), 'Value', q1);
        p{5*N+i} = addparameter(k{5*N+i}, strcat('c', int2str(5*N+i)), 'Value', q1);
        % add rates for Y <-> Gy
        p{6*N+i} = addparameter(k{6*N+i}, strcat('c', int2str(6*N+i)), 'Value', q4);
        p{7*N+i} = addparameter(k{7*N+i}, strcat('c', int2str(7*N+i)), 'Value', q4);
    end
    

    % set initial concentrations for x and y
    for i = 1 : N
        r{2*N+2*i-1}.Reactants(1).InitialAmount = initX * base * (1/N); % X
        r{2*N+2*i-1}.Reactants(2).InitialAmount = initY * base * (2/N); % Gy
        r{4*N+i}.Reactants(1).InitialAmount = initY * base * (2/N); % y
    end
    % set initial concentrations for excess gates
    for i = 1 : N
        r{2*i-1}.Reactants(2).InitialAmount = infConc;
        r{2*i}.Reactants(2).InitialAmount = infConc;
        r{2*N+2*i}.Reactants(2).InitialAmount = infConc;
        r{4*N+i}.Reactants(2).InitialAmount = infConc;
    end
    
    model
    model.Reactions
    model.Species

    % get solver config and set ode type, tolerance, stop time
    cs = getconfigset(model,'active');
    set(cs.SolverOptions, 'AbsoluteTolerance', absTol);
    set(cs.SolverOptions, 'RelativeTolerance', relTol);
    cs.SolverType = 'ode15s';
    cs.StopTime = stopTime;
    [t,X] = sbiosimulate(model);
    
    % plot
    x =  sum(X(:, 1:4:N*4), 2); y = X(:, 4*N+1) + sum(X(:, 4*N+6:4:4*N+4*N), 2);
    gy = X(:, 4*N+4)+ sum(X(:, 4*N+5:4:4*N+4*N), 2);
    figure(figNo); subplot(2, 1, 1);
    hold on; box on;
    plot(x./2, (y+gy)./2, 'LineWidth', 2.0);
    ylabel('Y'); xlabel('X');
    set(gca, 'LineWidth', 2.0); 
    
    subplot(2, 1, 2);
    box on; hold on;
    plot(t./60, x./2, 'LineWidth', 2.0);
    plot(t./60, (y+gy)./2, 'LineWidth', 2.0);
    ylabel('Concentration (nM)'); xlabel('Time (mins)');
    set(gca, 'LineWidth', 2.0); 
    legend('cat. X', 'cat. Y', 'approx. X', 'approx. Y'); 
end


function approx_link_cat_lv_oscillator_leaky_shadow(initX, initY, N, stopTime, figNo)
    % implements a DNA version of catalytic lotka volterra oscillator 
    %
    % initX is the initial concentration of X
    % initY is the initial concentration of Y
    % N is the number of catalytic reactions
    % stopTime is simulation time (eg. 100)
    % figNo is the figure # for plot

    % set all the default params and initialize model
    % default rate and concentration
    base = 1;

    % alpha and beta for scaling
    alpha = 1e-3;
    beta = 1e-8;

    % slow down all reactions by multiplying with a scaling factor 1e-3
    rate = (1/60);
    % fastest reaction are assumed to occur 
    infRate = 1e5*(rate); % different from approx cat
    % inf concentration of gates will be 1000 uM
    infConc = 1e5*(base);
    % Annihilation rate
    annRate = 1e6*(base); 
    
    % leak rate
    leak_rate = 1e-6;

    % leak and shadow flags
    leak = 1;
    shadow = 1;
    regular = 1;
    deactivate = 0; % Keep this 0 always!

    % allowed error rate tolerance value in pico molars
    absTol = 1e-12;
    relTol= 1e-12;
    
    % create a sym bio model for simulation
    model = sbiomodel('Approximate CRNs for LV oscillator with Leaks and Shadowing');
    r = cell(1, 28*N);
    k = cell(1, length(r));
    p = cell(1, length(r));
    
    %% Reactions
    % X + Y -> 2 Y
    % X -> 2 X
    % Y -> W
    for i = 1:N
        istr = int2str(i);
        inxt = int2str(max(1, rem(i+1, N+1)));
        % X + Y -> 2 Y
        r{2*i-1} = addreaction(model, strcat('x', istr, ' + Yg', istr, ' -> Iy', istr));
        r{2*i} = addreaction(model, strcat('Iy', istr, ' + Gy', istr, ' -> y', istr, ' + y', inxt));
        % X -> 2 X
        r{2*N + 2*i-1} = addreaction(model, strcat('x', istr, ' + Gxi', istr, ' -> Ix', istr));
        r{2*N + 2*i} = addreaction(model, strcat('Ix', istr, ' + Gix', istr, ' -> x', istr, ' + x', inxt));
        % Y -> W
        r{4*N + i} = addreaction(model, strcat('y', istr, ' + Gw1', istr, ' -> w'));
        r{5*N + i} = addreaction(model, strcat('Yg', istr, ' + Gw2', istr, ' -> w'));
        % Linker reactions
        r{6*N + i} = addreaction(model, strcat('y', istr, ' -> Yg', istr));
        r{7*N + i} = addreaction(model, strcat('Yg', istr, ' -> y', istr));
    end

    %% Shadow Circuit reactions
    for i = 1:N
        istr = int2str(i);
        inxt = int2str(max(1, rem(i+1, N+1)));
        % X + Y -> 2 Y
        r{8*N + 2*i-1} = addreaction(model, strcat('sx', istr, ' + sYg', istr, ' -> sIy', istr));
        r{8*N + 2*i} = addreaction(model, strcat('sIy', istr, ' + sGy', istr, ' -> sy', istr, ' + sy', inxt));
        % X -> 2 X
        r{10*N + 2*i-1} = addreaction(model, strcat('sx', istr, ' + sGxi', istr, ' -> sIx', istr));
        r{10*N + 2*i} = addreaction(model, strcat('sIx', istr, ' + sGix', istr, ' -> sx', istr, ' + sx', inxt));
        % Y -> W
        r{12*N + i} = addreaction(model, strcat('sy', istr, ' + sGw1', istr, ' -> w'));
        r{13*N + i} = addreaction(model, strcat('sYg', istr, ' + sGw2', istr, ' -> w'));
        % Linker reactions
        r{14*N + i} = addreaction(model, strcat('sy', istr, ' -> sYg', istr));
        r{15*N + i} = addreaction(model, strcat('sYg', istr, ' -> sy', istr));
    end

    %% Leak reactions
    for i = 1:N
        istr = int2str(i);
        inxt = int2str(max(1, rem(i+1, N+1)));
        % Yg1 -> Iy1
        r{16*N + i} = addreaction(model, strcat('Yg', istr, ' -> Iy', istr));
        r{17*N + i} = addreaction(model, strcat('sYg', istr, ' -> sIy', istr));
        r{18*N + i} = addreaction(model, strcat('Iy', istr, ' + sIy', istr, ' -> w'));
        
        % Gy1 -> y1 + y2
        r{19*N + i} = addreaction(model, strcat('Gy', istr, ' -> y', istr, ' + y', inxt));
        r{20*N + i} = addreaction(model, strcat('sGy', istr, ' -> sy', istr, ' + sy', inxt));
        r{21*N + i} = addreaction(model, strcat('y', istr, ' + sy', istr, ' -> w'));

        % Gix1 -> x1 + x2
        r{22*N + i} = addreaction(model, strcat('Gix', istr, ' -> x', istr, ' + x', inxt));
        r{23*N + i} = addreaction(model, strcat('sGix', istr, ' -> sx', istr, ' + sx', inxt));
        r{24*N + i} = addreaction(model, strcat('x', istr, ' + sx', istr, ' -> w'));

        % Gxi1 -> Ix1
        r{25*N + i} = addreaction(model, strcat('Gxi', istr, ' -> Ix', istr));
        r{26*N + i} = addreaction(model, strcat('sGxi', istr, ' -> sIx', istr));
        r{27*N + i} = addreaction(model, strcat('Ix', istr, ' + sIx', istr, ' -> w'));
 
    end
    

    %% Kinetic laws
    for i = 1:length(r)
        k{i} = addkineticlaw(r{i}, 'MassAction');
        k{i}.ParameterVariableNames = {strcat('c', int2str(i))};
    end

    %% Rate constants
    gamma = 2;
    bimol_scaled_rate = N*rate; % When xi = x/N, yi = y/N
    bimol = gamma*bimol_scaled_rate; % Accounting for the presence of y strand in both gate form and strand form
    unimol_scaled_rate = rate; % when xi = x/N, yi = y/N
    unimol = unimol_scaled_rate/infConc;

    for i = 1:N
        % X + Y -> 2 Y
        p{2*i-1} = addparameter(k{2*i-1}, strcat('c', int2str(2*i-1)), 'Value', regular*bimol);
        p{2*i} = addparameter(k{2*i}, strcat('c', int2str(2*i)), 'Value', regular*infRate);
        % X -> 2 X
        p{2*N + 2*i-1} = addparameter(k{2*N + 2*i-1}, strcat('c', int2str(2*N + 2*i-1)), 'Value', regular*unimol);
        p{2*N + 2*i} = addparameter(k{2*N + 2*i}, strcat('c', int2str(2*N + 2*i)), 'Value', regular*infRate);
        % Y -> w
        p{4*N + i} = addparameter(k{4*N + i}, strcat('c', int2str(4*N + i)), 'Value', regular*unimol);
        p{5*N + i} = addparameter(k{5*N + i}, strcat('c', int2str(5*N + i)), 'Value', regular*unimol);
        % Linker reactions
        p{6*N + i} = addparameter(k{6*N + i}, strcat('c', int2str(6*N + i)), 'Value', regular*infRate);
        p{7*N + i} = addparameter(k{7*N + i}, strcat('c', int2str(7*N + i)), 'Value', regular*infRate);

        % SX + SY -> 2 SY
        p{8*N + 2*i-1} = addparameter(k{8*N + 2*i-1}, strcat('c', int2str(8*N + 2*i-1)), 'Value', shadow*bimol);
        p{8*N + 2*i} = addparameter(k{8*N + 2*i}, strcat('c', int2str(8*N + 2*i)), 'Value', shadow*infRate);
        % SX -> 2 SX
        p{10*N + 2*i-1} = addparameter(k{10*N + 2*i-1}, strcat('c', int2str(10*N + 2*i-1)), 'Value', shadow*unimol);
        p{10*N + 2*i} = addparameter(k{10*N + 2*i}, strcat('c', int2str(10*N + 2*i)), 'Value', shadow*infRate);
        % SY -> w
        p{12*N + i} = addparameter(k{12*N + i}, strcat('c', int2str(12*N + i)), 'Value', shadow*unimol);
        p{13*N + i} = addparameter(k{13*N + i}, strcat('c', int2str(13*N + i)), 'Value', shadow*unimol);
        % Linker reactions
        p{14*N + i} = addparameter(k{14*N + i}, strcat('c', int2str(14*N + i)), 'Value', shadow*infRate);
        p{15*N + i} = addparameter(k{15*N + i}, strcat('c', int2str(15*N + i)), 'Value', shadow*infRate);

        % Leaks
        % Yg -> Iy
        p{16*N + i} = addparameter(k{16*N + i}, strcat('c', int2str(16*N + i)), 'Value', leak*regular*leak_rate);
        p{17*N + i} = addparameter(k{17*N + i}, strcat('c', int2str(17*N + i)), 'Value', leak*shadow*leak_rate);
        p{18*N + i} = addparameter(k{18*N + i}, strcat('c', int2str(18*N + i)), 'Value', leak*regular*shadow*infRate*1e1);
        % Gy1 -> y1 + y2
        p{19*N + i} = addparameter(k{19*N + i}, strcat('c', int2str(19*N + i)), 'Value', leak*regular*leak_rate);
        p{20*N + i} = addparameter(k{20*N + i}, strcat('c', int2str(20*N + i)), 'Value', leak*shadow*leak_rate);
        p{21*N + i} = addparameter(k{21*N + i}, strcat('c', int2str(21*N + i)), 'Value', leak*regular*shadow*infRate*1e1); % CHECK!
        % Gix1 -> x1 + x2
        p{22*N + i} = addparameter(k{22*N + i}, strcat('c', int2str(22*N + i)), 'Value', leak*regular*leak_rate);
        p{23*N + i} = addparameter(k{23*N + i}, strcat('c', int2str(23*N + i)), 'Value', leak*shadow*leak_rate);
        p{24*N + i} = addparameter(k{24*N + i}, strcat('c', int2str(24*N + i)), 'Value', leak*regular*shadow*infRate*1e1); % CHECK!
        % Gxi1 -> Ix1
        p{25*N + i} = addparameter(k{25*N + i}, strcat('c', int2str(25*N + i)), 'Value', leak*regular*leak_rate);
        p{26*N + i} = addparameter(k{26*N + i}, strcat('c', int2str(26*N + i)), 'Value', leak*regular*leak_rate);
        p{27*N + i} = addparameter(k{27*N + i}, strcat('c', int2str(27*N + i)), 'Value', leak*regular*shadow*infRate*1e1); % CHECK!

    end

    %% Initial concentrations
    for i = 1:N
        r{2*i-1}.Reactants(1).InitialAmount = initX/N; % x
        r{2*i}.Reactants(2).InitialAmount = infConc; % Gy

        r{2*N + 2*i-1}.Reactants(2).InitialAmount = infConc; % Gxi
        r{2*N + 2*i}.Reactants(2).InitialAmount = infConc; % Gix

        r{4*N + i}.Reactants(1).InitialAmount = initY/N; % y_i
        r{4*N + i}.Reactants(2).InitialAmount = infConc; % Gw1
        r{5*N + i}.Reactants(2).InitialAmount = infConc; % Gw2

        % Shadow circuit
        r{8*N + 2*i}.Reactants(2).InitialAmount = infConc; % sGy
        r{10*N + 2*i-1}.Reactants(2).InitialAmount = infConc; % sGxi
        r{10*N + 2*i}.Reactants(2).InitialAmount = infConc; % sGix
        r{12*N + i}.Reactants(2).InitialAmount = infConc; % sGw1
        r{13*N + i}.Reactants(2).InitialAmount = infConc; % sGw2

    end

    
    model
    model.Reactions
    model.Species

    % get solver config and set ode type, tolerance, stop time
    cs = getconfigset(model,'active');
    set(cs.SolverOptions, 'AbsoluteTolerance', absTol);
    set(cs.SolverOptions, 'RelativeTolerance', relTol);
    cs.SolverType = 'ode15s';
    cs.StopTime = stopTime;
    [t,X, names] = sbiosimulate(model);
    
    function [val] = get_plot_values(X, names, s, n)
        val = zeros(length(t), 1);
        
        for iter = 1:length(names)
            if strncmp(names{iter}, s, n)
                disp(strcat(names{iter}, ' should be ', s));
                val = val + X(:, iter);
            end
        end
    end
    
    x = get_plot_values(X, names, 'x', 1);
    y = get_plot_values(X, names, 'y', 1);
    yg = get_plot_values(X, names, 'Yg', 2);
    sx = get_plot_values(X, names, 'sx', 2);
    sy = get_plot_values(X, names, 'sy', 2);
    syg = get_plot_values(X, names, 'syg', 3);


    figure(figNo); subplot(2, 1, 1);
    hold on; box on;
    plot(x, (y+yg), 'LineWidth', 2.0);
    ylabel('Y'); xlabel('X');
    set(gca, 'LineWidth', 2.0); 
    
    subplot(2, 1, 2);
    box on; hold on;
    plot(t./60, x, '-', 'LineWidth', 2.0);
    plot(t./60, (y+yg), ':', 'LineWidth', 2.0);
%     plot(t./60, y, '-', 'LineWidth', 2.0);
%     plot(t./60, yg, ':', 'LineWidth', 2.0);
%     plot(t./60, sx, '-','LineWidth', 2.0);
%     plot(t./60, (sy+syg), ':', 'LineWidth', 2.0);
    ylabel('Concentration (nM)'); xlabel('Time (mins)');
    set(gca, 'LineWidth', 2.0); 
    legend('approx. X', 'approx. Y'); 
end