%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Shalin Shah
% Date created: 04/01/2020
% Affiliation: Duke University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
concX = 4; concY = 3; stopTime = 1000; nCatalytic = 3;

% compare ideal lotka-volterra simulation with approximation
ideal_amp(concX, concY, stopTime, 1);
cat_amp(concX, concY, nCatalytic, stopTime, 1);
%approx_cat_amp(concX, concY, nCatalytic, stopTime, 1);

function ideal_amp(concA, concB, stopTime, figNo)
    % enter A concentration assuming all gates are excess
    %
    % stopTime is the total simulation time
    % figNo is the display figure number
    
    % rate of reactions 
    rate = (1/60);
    % allowed error rate tolerance value in pico molars
    absTol = 1e-12;
    relTol = 1e-12;
    
    % Create the amplification model
    model = sbiomodel('Ideal CRNs for LV oscillator');

    % Add reaction set X + Y -> 2Y to the solver
    r1 = addreaction(model,'X + Y -> Y + Y');
        
    % Set the Kinetic Law for Reactions.
    k1 = addkineticlaw(r1, 'MassAction'); 
    
    p1 = addparameter(k1, 'c1', 'Value', rate);
    
    % Set the Kinetic Law for Reactions.
    k1.ParameterVariableNames = {'c1'};
    
    % Set initial amounts for species 
    r1.Reactants(1).InitialAmount = concA;        % A
    r1.Reactants(2).InitialAmount = concB;        % B

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
    figure(figNo);
    box on; hold on;
    plot(t, X(:, 1), 'LineWidth', 2.0);
    plot(t, X(:, 2), 'LineWidth', 2.0);
%     plot(t, X(:, 3), 'LineWidth', 2.0);
%     plot(t, sum(X(:, 1:2), 2), 'LineWidth', 2.0);
    ylabel('Concentration'); xlabel('Time');
    set(gca, 'LineWidth', 2.0); legend('X', 'Y');
end

function cat_amp(initX, initY, N, stopTime, figNo)
    % implements an ideal lotka volterra oscillator with one modified
    % reaction
    %
    % initX is the initial concentration of X
    % initY is the initial concentration of Y
    % N is the number of catalytic reactions
    % stopTime is simulation time (eg. 100)
    % figNo is the figure # for plot

    % HOW TO SCALE FOR THE CATALYTIC REACTIONS
    % 1- set all X init concentration to initX/N
    % 2- set end time to to t * N
    % 3- divide all time units to t/N
    
    % default rate and concentration
    rate = (1/60); % one per minute
    base = 1;
    
    % allowed error rate tolerance value in pico molars
    absTol = 1e-12;
    relTol = 1e-12;
    
    % create a sym bio model for simulation
    model = sbiomodel('Approximate CRNs for LV oscillator');
    r = cell(1, 3*N);
    k = cell(1, length(r));
    p = cell(1, length(r));
    
    % add catalytic approximations for X + Y -> 2Y
    for i = 1:N
        r{i} = addreaction(model, strcat('x', int2str(i), ' -> x', int2str(i)));
    end
    for i = N+1:2*N
        r{i} = addreaction(model, strcat('y', int2str(i-N), ' -> y', int2str(i-N)));
    end
    
    for i = 2*N+1:3*N - 1
        r{i} = addreaction(model, strcat('x', int2str(i-2*N), ' + y', int2str(i-2*N), ...
            ' -> y', int2str(i-2*N),' + y', int2str(i+1-2*N))); 
    end
    r{3*N} = addreaction(model, strcat('x', int2str(N), ' + y', int2str(N), ...
            ' -> y', int2str(N),' + y', int2str(1))); 
    
    % add all the rate laws as mass action kinetics
    for i = 1:N
        k{i} = addkineticlaw(r{i}, 'MassAction'); 
        k{i}.ParameterVariableNames = {strcat('c', int2str(i))};
        p{i} = addparameter(k{i}, strcat('c', int2str(i)), 'Value', rate * 1);
    end
    
    for i = N+1:3*N
        k{i} = addkineticlaw(r{i}, 'MassAction'); 
        k{i}.ParameterVariableNames = {strcat('c', int2str(i))};
        p{i} = addparameter(k{i}, strcat('c', int2str(i)), 'Value', rate);
    end
    
    % set SCALED initial concentrations of X, Y
    for i = 2*N+1:3*N
        r{i}.Reactants(1).InitialAmount = initX * base * (1/N);
        r{i}.Reactants(2).InitialAmount = initY * base * (1/N);
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
    figure(figNo);
    box on; hold on;
    plot(t, sum(X(:, 1:N), 2), ':', 'LineWidth', 2.0);
    plot(t, sum(X(:, N+1:2*N), 2), ':', 'LineWidth', 2.0);
    ylabel('Concentration'); xlabel('Time');
    set(gca, 'LineWidth', 2.0); legend('X', 'Y', 'cat_X', 'cat_Y');
end

function approx_cat_amp(initX, initY, N, stopTime, figNo)
    
    % default rate and concentration
    base = 1;
    % slow down all reactions
    rate = 1/60;
    % fastest reaction are assumed to occur 
    infRate = 0.5*(rate);
    % inf concentration of gates will be 1000 uM
    infConc = 10*(base);

    % allowed error rate tolerance value in pico molars
    absTol = 1e-12;
    relTol = 1e-12;

    % leak rate
    leak_rate = 5e-5;

    % switch shadow circuit flag. 1 is ON, 0 is OFF
    shadow = 1;

    % switch leak flag. 1 is ON, 0 is OFF
    leak = 1;
    
    % create a symbio model for simulation
    model = sbiomodel('Approximate CRNs for a Autocatalytic Bimolecular amplifier');
    
    % 2N reactions for the regular circuit
    % N reactions for the leaks
    % 2N reactions for the shadow circuit
    % N reactions for shadow leaks
    % N reactions for annihilation
    r = cell(1, 7*N);
    k = cell(1, length(r));
    p = cell(1, length(r));

    % DNA implementations for each catalytic approximation reaction
    % Xi + Yi -> Ixyi
    % Ixyi + Gxyi -> Yi + Y_(i+1)
    
    % Reactions in the regular circuit
    for i = 1: N-1
        r{2*i - 1} = addreaction(model, strcat('y', int2str(i), ' + x' ...
            , int2str(i), ' -> Ixy', int2str(i)));
        r{2*i} = addreaction(model, strcat('Ixy', int2str(i), ' + Gxy', ...
            int2str(i), ' -> y', int2str(i), ' + y', int2str(i+1)));
    end
    r{2*N - 1} = addreaction(model, strcat('y', int2str(N), ' + x', ...
        int2str(N), ' -> Ixy', int2str(N)));
    r{2*N} = addreaction(model, strcat('Ixy', int2str(N), ' + Gxy', ...
        int2str(N), ' -> y', int2str(N), ' + y1'));

    % Leak reactions in the regular circuit
    for i = 1:N-1 %2N+1:3N
        r{2*N + i} = addreaction(model, strcat('Gxy', int2str(i), ...
            ' -> y', int2str(i),' + y', int2str(i+1)));
    end
    r{3*N} = addreaction(model, strcat('Gxy', int2str(N), ...
        ' -> y', int2str(N), ' + y1'));
    
    % Reactions in the shadow circuit
    for i = 1: N-1 %3N+1:5N
        r{3*N + 2*i -1} = addreaction(model, strcat('sy', int2str(i), ' + sx' ...
            , int2str(i), ' -> Isxy', int2str(i)));
        r{3*N + 2*i} = addreaction(model, strcat('Isxy', int2str(i), ' + Gsxy', ...
            int2str(i), ' -> sy', int2str(i), ' + sy', int2str(i+1)));
    end
    r{5*N - 1} = addreaction(model, strcat('sy', int2str(N), ' + sx', ...
        int2str(N), ' -> Isxy', int2str(N)));
    r{5*N} = addreaction(model, strcat('Isxy', int2str(N), ' + Gsxy', ...
        int2str(N), ' -> sy', int2str(N), ' + sy1'));
    
    % Leak reactions in the shadow circuit
    for i = 1:N-1 %5N+1:6N
        r{5*N + i} = addreaction(model, strcat('Gsxy', int2str(i), ...
            ' -> sy', int2str(i),' + sy', int2str(i+1)));
    end
    r{6*N} = addreaction(model, strcat('Gsxy', int2str(N), ...
        ' -> sy', int2str(N), ' + sy1'));

    % Annihilation reactions
    for i = 1: N %6N+1:7N
        r{6*N + i} = addreaction(model, strcat('y', int2str(i), ...
            ' + sy', int2str(i), ' -> W',int2str(i)));
    end

    % Kinetics
    for i = 1: 7*N
        k{i} = addkineticlaw(r{i}, 'MassAction');
        k{i}.ParameterVariableNames = {strcat('c', int2str(i))};
    end

    % Assigning Rate parameters to all the reactions
    q1 = N*rate; q2 = infRate;
    
    % Rate params for the regular circuit
    for i = 1:2:2*N-1 %1, 3, 5, 
        p{i} = addparameter(k{i}, strcat('c', int2str(i)), 'Value', q1);
        p{i+1} = addparameter(k{i+1}, strcat('c', int2str(i+1)), 'Value', q2); % check
    end

    % Rate params for the leaks in the regular circuit
    for i = 2*N+1:1:3*N
        p{i} = addparameter(k{i}, strcat('c', int2str(i)), 'Value', leak*leak_rate);
    end

    % Rate for the shadow circuit
    for i = 3*N + 1: 2: 5*N - 1
        istr = int2str(i);
        p{i} = addparameter(k{i}, strcat('c', istr), 'Value', shadow*q1);
        p{i+1} = addparameter(k{i+1}, strcat('c', int2str(i+1)), 'Value', shadow*q2); % check
    end

    % Rate for the shadow leaks
      for i = 5*N + 1: 1: 6*N
        p{i} = addparameter(k{i}, strcat('c', int2str(i)), 'Value', leak*shadow*leak_rate);
      end


    % Rate for annihilation reactions
    for i = 6*N + 1: 7*N
        istr = int2str(i);
        p{i} = addparameter(k{i}, strcat('c', istr), 'Value', shadow*infRate); % check
    end

    % Set initial concentration of x and y
    for i = 1: 2: 2*N-1
        % y + x -> Ixy = 1, 3, 5, ...
        r{i}.Reactants(1).InitialAmount = initY*(1/N);
        r{i}.Reactants(2).InitialAmount = initX*(1/N);
%       r{i}.Reactants(1).InitialAmount = 0.0; % check.
    end
    
    % Set initial concentration for regular gates 
    for i = 2: 2: 2*N
        % Ixy + Gxy -> y + y = 2, 4, 6, ... 
        r{i}.Reactants(2).InitialAmount = infConc;
    end

    % Set initial concentrations for s
    for i = 3*N+1: 2 : 5*N-1
        % sy + sx -> Isxy = 3N+1, 3N+3, ... 5N-1
        r{i}.Reactants(1).InitialAmount = 0.0;
        r{i}.Reactants(2).InitialAmount = 0.0;
    end
    % Set initial concentration for shadow gates
    for i = 3*N+2: 2: 5*N
        % Isxy + Gsxy -> sy + sy = 3N+2, 3N+4, ... 5N
        r{i}.Reactants(2).InitialAmount = infConc; % check. This is good!
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
    x = sum(X(:, 2:4:11), 2);
    y = sum(X(:, 1:4:10), 2);
    %sy = sum(X(:, 13:4:22), 2);
    W = sum(X(:, 25:27), 2);
    plot(t, x, ':', LineWidth=3.0);
    plot(t, y, '-','LineWidth', 3.0);
    %plot(t, sy, '-','LineWidth', 3.0);
    plot(t, W, '-','LineWidth', 3.0);
    ylabel('Concentration (nM)'); xlabel('Time (mins)');
    set(gca, 'LineWidth', 2.0); legend('ideal.X', 'ideal.Y', 'cat.X', 'cat.Y','dna.X','dna.Y','dna.W');    
    
end