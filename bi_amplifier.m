%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Shalin Shah
% Date created: 04/01/2020
% Affiliation: Duke University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
concX = 4; concY = 3; stopTime = 300; nCatalytic = 3;

% compare ideal lotka-volterra simulation with approximation
ideal_amp(concX, concY, stopTime, 1)
cat_amp(concX, concY, nCatalytic, stopTime, 1);
% approx_cat_amp(concX, concY, nCatalytic, stopTime, 1);

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
    set(gca, 'LineWidth', 2.0); legend('X', 'Y', 'X+Y')
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
    set(gca, 'LineWidth', 2.0); legend('X', 'Y', 'X+Y')
end

function approx_cat_amp(initX, initY, N, stopTime, figNo)
    % implements a DNA version of catalytic lotka volterra oscillator 
    %
    % initX is the initial concentration of X
    % N is the number of catalytic reactions
    % stopTime is simulation time (eg. 100)
    % figNo is the figure # for plot

    % set all the default params and initialize model
    % default rate and concentration
    base = 1; % (or 1 nano molar)
    crazy = 1e-5;
    % slow down all reactions by multiplying with a scaling factor 1e-3
    rate = (1/60); % (or 1 per minute)
    % fastest reaction are assumed to occur 
    infRate = 1e1*(rate);
    % inf concentration of gates will be 1000 uM
    infConc = 1e1*base;
    
    % allowed error rate tolerance value in pico molars
    absTol = 1e-12;
    relTol = 1e-12;
    
    % create a sym bio model for simulation
    model = sbiomodel('Approximate CRNs for autocatalytic CRNs');
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
    
    % add two-step catalytic approximations for X -> X + X step
    % For N catalytic reaction, we have 2N reactions
    for i = 1 : N-1
        r{2*N+2*i-1} = addreaction(model, strcat('x', int2str(i), ' + y', int2str(i), ...
            ' -> Ixy', int2str(i))); 
        r{2*N+2*i} = addreaction(model, strcat('Ixy', int2str(i), ...
            ' + Gxy', int2str(i),' -> y', int2str(i),' + y', int2str(i+1))); 
    end
    r{4*N-1} = addreaction(model, strcat('x', int2str(N), ' + y', int2str(N), ...
            ' -> Ixy', int2str(N)));
    r{4*N} = addreaction(model, strcat('Ixy', int2str(N), ...
            ' + Gxy', int2str(N),' -> y', int2str(N),' + y1')); 
    
    % add kinetic rates, first 2N are useless 
    for i = 1 : 4*N
        k{i} = addkineticlaw(r{i}, 'MassAction'); 
        k{i}.ParameterVariableNames = {strcat('c', int2str(i))};
    end
    
    % add useless params
    for i = 1 : 2*N
        p{i} = addparameter(k{i}, strcat('c', int2str(i)), 'Value', crazy);
    end
    % add the real params
    q1 = N*rate; q2 = infRate;
    for i = 1 : N
        p{2*N+2*i-1} = addparameter(k{2*N+2*i-1}, strcat('c', int2str(2*N+2*i-1)), 'Value', q1);
        p{2*N+2*i} = addparameter(k{2*N+2*i}, strcat('c', int2str(2*N+2*i)), 'Value', q2);
    end

    % set initial concentrations for x
    for i = 2*N + 1 : 2 : 4*N
        r{i}.Reactants(1).InitialAmount = initX * (1/N);
        r{i}.Reactants(2).InitialAmount = initY * (1/N);
    end
    % set initial concentrations for excess gates
    for i = 1 : 2 : 2*N
        r{2*N+i+1}.Reactants(2).InitialAmount = infConc;
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
    plot(t, sum(X(:, 1:N), 2), ':','LineWidth', 2.0);   
    plot(t, sum(X(:, 1+N:2*N), 2), ':','LineWidth', 2.0);   
    ylabel('Concentration (nM)'); xlabel('Time (mins)');
    set(gca, 'LineWidth', 2.0);  legend('X', 'Y', 'cat. X', 'cat. Y', 'approx. X', 'approx. Y')
end