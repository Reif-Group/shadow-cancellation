clear; clc;
concX = 2; stopTime = 300; nCatalytic = 5;

% compare ideal lotka-volterra simulation with approximation
ideal_amp(concX, stopTime, 1)
cat_amp(concX, nCatalytic, stopTime, 1);
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
    
    % plot
    figure(figNo);
    box on; hold on;
    plot(t, sum(X, 2), ':', 'LineWidth', 2);
    ylabel('Concentration (nM)'); xlabel('Time (mins)'); 
    set(gca, 'LineWidth', 2.0); 
end

function approx_cat_amp(initX, N, stopTime, figNo)
    % implements a DNA version of catalytic lotka volterra oscillator 
    %
    % initX is the initial concentration of X
    % N is the number of catalytic reactions
    % stopTime is simulation time (eg. 100)
    % figNo is the figure # for plot

    % set all the default params and initialize model
    % default rate and concentration
    base = 1;
    % slow down all reactions by multiplying with a scaling factor 1e-3
    rate = (1/60); % (or 1 per minute)
    % fastest reaction are assumed to occur 
    infRate = 1e1*(rate);
    % inf concentration of gates will be 1000 uM
    infConc = 1e2*base;
    
    % allowed error rate tolerance value in pico molars
    absTol = 1e-12;
    relTol = 1e-12;
    
    % create a sym bio model for simulation
    model = sbiomodel('Approximate CRNs for autocatalytic CRNs');
    r = cell(1, 2*N);
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
    
    for i = 1 : 2*N
        k{i} = addkineticlaw(r{i}, 'MassAction'); 
        k{i}.ParameterVariableNames = {strcat('c', int2str(i))};
    end
    
    q1 = rate/infConc; q2 = infRate;
    for i = 1:N
        p{2*i-1} = addparameter(k{2*i-1}, strcat('c', int2str(2*i-1)), 'Value', q1);
        p{2*i} = addparameter(k{2*i}, strcat('c', int2str(2*i)), 'Value', q2);
    end

    % set initial concentrations for x
    for i = 1 : 2 : 2*N
        r{i}.Reactants(1).InitialAmount = initX * (1/N);
    end
    % set initial concentrations for excess gates
    for i = 1 : 2*N
        r{i}.Reactants(2).InitialAmount = infConc;
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
    x = sum(X(:, 1:4:N*4), 2);
    plot(t, x, ':','LineWidth', 2.0);   
    ylabel('Concentration (nM)'); xlabel('Time (mins)');
    set(gca, 'LineWidth', 2.0); legend('auto. X', 'cat. X', 'approx. X');
end