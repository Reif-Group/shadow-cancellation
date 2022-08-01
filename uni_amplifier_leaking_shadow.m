%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Shalin Shah
% Date created: 04/01/2020
% Affiliation: Duke University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
concX = 2; stopTime = 6000; nCatalytic = 5;

% compare ideal lotka-volterra simulation with approximation
% ideal_amp(concX, stopTime, 1)
% cat_amp(concX, nCatalytic, stopTime, 1);
approx_cat_amp(concX, nCatalytic, stopTime, 1);

function ideal_amp(concA, stopTime, figNo)
    % enter A concentration assuming all gates are excess
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

    % Add reaction set X -> 2X to the solver
    r1 = addreaction(model,'X -> X + X');
        
    % Set the Kinetic Law for Reactions.
    k1 = addkineticlaw(r1, 'MassAction'); 
    
    p1 = addparameter(k1, 'c1', 'Value', rate);
    
    % Set the Kinetic Law for Reactions.
    k1.ParameterVariableNames = {'c1'};
    
    % Set initial amounts for species 
    r1.Reactants(1).InitialAmount = concA;        % A

    % Display the Completed Model Objects
    model;

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

%     % plot X over the time
    figure(figNo);
    box on; hold on;
    plot(t, X, 'LineWidth', 2.0);
    ylabel('Concentration'); xlabel('Time');
    set(gca, 'LineWidth', 2.0); 
end

function cat_amp(initX, N, stopTime, figNo)
    % implements an ideal lotka volterra oscillator with one modified
    % reaction
    %
    % initX is the initial concentration of X
    % N is the number of catalytic reactions
    % stopTime is simulation time (eg. 100)
    % figNo is the figure # for plot

    % default rate and concentration
    rate = 1/60; % per second (or 1 per minute)
    base = 1; % nano moles
    
    % allowed error rate tolerance value in pico molars
    absTol = 1e-12;
    relTol = 1e-12;
    
    % create a sym bio model for simulation
    model = sbiomodel('Approximate CRNs for LV oscillator');
    r = cell(1, N);
    k = cell(1, length(r));
    p = cell(1, length(r));
    
    % add catalytic approximations for last step
    for i = 1 : N-1
        r{i} = addreaction(model, strcat('x', int2str(i), ' -> x', int2str(i),' + x', int2str(i+1))); 
    end
    r{N} = addreaction(model, strcat('x', int2str(N), ' -> ', 'x', int2str(N), ' + x1')); 
    
    % add all the rate laws as mass action kinetics
    for i = 1:N
        k{i} = addkineticlaw(r{i}, 'MassAction'); 
        k{i}.ParameterVariableNames = {strcat('c', int2str(i))};
        p{i} = addparameter(k{i}, strcat('c', int2str(i)), 'Value', rate);
    end
    
    % set initial concentrations
    for i = 1:N
        r{i}.Reactants(1).InitialAmount = initX * base * (1/N);
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
%     
%     % plot
    figure(figNo);
    box on; hold on;
    plot(t, sum(X, 2), ':', 'LineWidth', 2);
    ylabel('Concentration (nM)'); xlabel('Time (mins)'); 
    set(gca, 'LineWidth', 2.0); 
end

function approx_cat_amp(initX, N, stopTime, figNo)
    
    % default rate and concentration
    base = 1;
    % slow down all reactions
    rate = 1/60;
    % fastest rate reactions
    infRate = 1e2*(rate);
    % inf concentration of the gates will be 1000uM
    infConc = 1e3*(base);

    % allowed error rate tolerance value in pico molars
    absTol = 1e-12
    relTol = 1e-12;

    % leak rate
    leak_rate = 1e-4;

    % switch shadow circuit flag. 1 is ON, 0 is OFF
    shadow = 1; 

    % switch leak flag. 1 is ON, 0 is OFF
    leak = 1;
    
    % create a symbio model for simulation
    model = sbiomodel('Approximate CRNs for a Autocatalytic Unimolecular amplifier');
    
    % 2N reactions for the regular circuit
    % N reactions for the leaks
    % 2N reactions for the shadow circuit
    % N reactions for shadow leaks
    % N reactions for annihilation
    r = cell(1, 9*N)
    k = cell(1, length(r));
    p = cell(1, length(r));

    % DNA implementations for each catalytic approximation reaction
    % Xi + Gxi1 -> Ixi
    % Ixi + Gxi2 -> Xi + X_(i+1)
    
    % Reactions in the regular circuit
    for i = 1: N-1
        r{2*i - 1} = addreaction(model, strcat('x', int2str(i), ' + Gx' ...
            , int2str(i), '1 -> Ix', int2str(i)));
        r{2*i} = addreaction(model, strcat('Ix', int2str(i), ' + Gx', ...
            int2str(i), '2 -> x', int2str(i), ' + x', int2str(i+1)));
    end
    r{2*N - 1} = addreaction(model, strcat('x', int2str(N), ' + Gx', ...
        int2str(N), '1 -> Ix', int2str(N)));
    r{2*N} = addreaction(model, strcat('Ix', int2str(N), ' + Gx', ...
        int2str(N), '2 -> x', int2str(N), ' + x1'));

    % Leak reactions in the regular circuit
    for i = 1:N-1
        r{2*N + 2*i - 1} = addreaction(model, strcat('Gx', int2str(i), ...
            '1 -> Ix', int2str(i)));
        r{2*N + 2*i} = addreaction(model, strcat('Gx', int2str(i), ...
            '2 -> x', int2str(i),' + x', int2str(i+1)));
    end
    r{2*N + 2*N - 1} = addreaction(model, strcat('Gx', int2str(N), ...
        '1 -> Ix', int2str(N)));
    r{2*N + 2*N} = addreaction(model, strcat('Gx', int2str(N), ...
        '2 -> x', int2str(N), ' + x1'))'
    
    % Reactions in the shadow circuit
    for i = 1: N-1
        r{4*N + 2*i - 1} = addreaction(model, strcat('s', int2str(i), ...
            ' + Gs', int2str(i), '1 -> Is', int2str(i)));
        r{4*N + 2*i} = addreaction(model, strcat('Is', int2str(i), ...
            ' + Gs', int2str(i), '2 -> s', int2str(i), ' + s', int2str(i+1)));
    end
    r{4*N + 2*N - 1} = addreaction(model, strcat('s', int2str(N), ' + Gs', ...
        int2str(N), '1 -> Is', int2str(N)));
    r{4*N + 2*N} = addreaction(model, strcat('Is', int2str(N), ' + Gs', ...
        int2str(N), '2 -> s', int2str(N), ' + s1'));
    
    % Leak reactions in the shadow circuit
    for i = 1:N-1
        r{6*N + 2*i - 1} = addreaction(model, strcat('Gs', int2str(i), ...
            '1 -> Is', int2str(i)));
        r{6*N + 2*i} = addreaction(model, strcat('Gs', int2str(i), ...
            '2 -> s', int2str(i),' + s', int2str(i+1)));
    end
    r{6*N + 2*N - 1} = addreaction(model, strcat('Gs', int2str(N), ...
        '1 -> Is', int2str(N)));
    r{6*N + 2*N} = addreaction(model, strcat('Gs', int2str(N), ...
        '2 -> s', int2str(N), ' + s1'))'

    % Annihilation reactions
    for i = 1: N
        r{8*N + i} = addreaction(model, strcat('x', int2str(i), ...
            ' + s', int2str(i), ' -> W'));
    end

    % Kinetics
    for i = 1: 9*N
        k{i} = addkineticlaw(r{i}, 'MassAction');
        k{i}.ParameterVariableNames = {strcat('c', int2str(i))};
    end

    % Assigning Rate parameters to all the reactions
    q1 = rate/infConc; q2 = infRate;
    
    % Rate params for the regular circuit
    for i = 1:2:2*N-1
        p{i} = addparameter(k{i}, strcat('c', int2str(i)), 'Value', q1);
        p{i+1} = addparameter(k{i+1}, strcat('c', int2str(i+1)), 'Value', q2); % check
    end

    % Rate params for the leaks in the regular circuit
    for i = 2*N + 1:2:4*N - 1
        p{i} = addparameter(k{i}, strcat('c', int2str(i)), 'Value', leak*leak_rate*q1);
        p{i+1} = addparameter(k{i+1}, strcat('c', int2str(i + 1)), 'Value', leak*leak_rate*q2);
    end


    % Rate for the shadow circuit
    for i = 4*N + 1: 2: 6*N - 1
        istr = int2str(i);
        p{i} = addparameter(k{i}, strcat('c', istr), 'Value', shadow*q1);
        p{i+1} = addparameter(k{i+1}, strcat('c', int2str(i+1)), 'Value', shadow*q2); % check
    end

    % Rate for the shadow leaks
      for i = 6*N + 1: 2: 8*N - 1
        p{i} = addparameter(k{i}, strcat('c', int2str(i)), 'Value', leak*shadow*leak_rate*q1);
        p{i+1} = addparameter(k{i+1}, strcat('c', int2str(i + 1)), 'Value', leak*shadow*leak_rate*q2);
      end


    % Rate for annihilation reactions
    for i = 8*N + 1: 9*N
        istr = int2str(i);
        p{i} = addparameter(k{i}, strcat('c', istr), 'Value', shadow*infRate); % check
    end

    % Set initial concentration of x
    for i = 1: 2: 2*N
        r{i}.Reactants(1).InitialAmount = initX*(1/N);
%         r{i}.Reactants(1).InitialAmount = 0.0; % check.
    end
    
    % Set initial concentration for regular gates 
    for i = 1: 2*N
        r{i}.Reactants(2).InitialAmount = infConc;
    end

    % Set initial concentrations for s
    for i = 4*N + 1: 2 : 6*N
        r{i}.Reactants(1).InitialAmount = 0.0;
    end
    % Set initial concentration for shadow gates
    for i = 4*N + 1: 6*N
        r{i}.Reactants(2).InitialAmount = infConc*1e-2 % check. This is good!
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
    disp(names);
    % plot
    figure(figNo);
    box on; hold on;
    s = sum(X(:, 21:4:20+N*4), 2); % changed, recheck.
    x = sum(X(:, 1:4:N*4), 2);
    plot(t, x, '-','LineWidth', 3.0);
    plot(t, s, ':', LineWidth=3.0)
    ylabel('Concentration (nM)'); xlabel('Time (mins)');
    set(gca, 'LineWidth', 2.0); %legend('auto. X', 'cat. X', 'approx. X', 'leak. S');    
%     set(gca, 'LineWidth', 2.0); legend('cat.X', 'leak.S');    
    
end
