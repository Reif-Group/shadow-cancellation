clear; clc;
concX = 4; concY = 2; stopTime = 15; nCatalytic = 50;

% compare ideal lotka-volterra simulation with approximation
ideallvoscillator(concX, concY, stopTime, 1)
approxlvoscillator(concX, concY, nCatalytic, stopTime, 2);

function ideallvoscillator(concA, concB, stopTime, figNo)
    % enter A and B concentration assuming all gates are excess
    %
    % stopTime is the total simulation time
    % figNo is the display figure number
    
    % rate of reactions 
    rate = 1;
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
    plot(t, X(:,1), 'LineWidth', 2.0);
    plot(t, X(:,2), 'LineWidth', 2.0);
    ylabel('Concentration'); xlabel('Time');
    set(gca, 'LineWidth', 2.0); 
    legend('X', 'Y')
end

function approxlvoscillator(initX, initY, N, stopTime, figNo)
    % implements an ideal lotka volterra oscillator with one modified
    % reaction
    %
    % initX is the initial concentration of X
    % initY is the initial concentration of Y
    % N is the number of catalytic reactions
    % stopTime is simulation time (eg. 100)
    % figNo is the figure # for plot

    % default rate and concentration
    rate = 1;
    base = 1;
    
    % allowed error rate tolerance value in pico molars
    absTol = 1e-12;
    relTol = 1e-12;
    
    % create a sym bio model for simulation
    model = sbiomodel('Approximate CRNs for LV oscillator');
    r = cell(1, N+3);
    k = cell(1, length(r));
    p = cell(1, length(r));
    
    % add two ideal reactions for LV
    r{1} = addreaction(model,'X + Y -> Y + Y');
    r{2} = addreaction(model,'Y -> waste'); 
    
    % add catalytic approximations for last step
    for i = 1 : N-1
        r{i+2} = addreaction(model, strcat('x', int2str(i), ' -> x', int2str(i),' + x', int2str(i+1))); 
    end
    % add the final reaction to release X
    r{N+2} = addreaction(model, strcat('x', int2str(N), ' -> x', int2str(N),' + X'));
    r{N+3} = addreaction(model, strcat('x1 -> phi'));
   
    % add all the rate laws as mass action kinetics
    for i = 1:N+3
        k{i} = addkineticlaw(r{i}, 'MassAction'); 
        k{i}.ParameterVariableNames = {strcat('c', int2str(i))};
    end
    
    for i = 1:N+3
        p{i} = addparameter(k{i}, strcat('c', int2str(i)), 'Value', rate);
    end
%     p{N+3} = addparameter(k{N+3}, strcat('c', int2str(N+3)), 'Value', 1/10000*rate);

    % set initial concentrations
    r{1}.Reactants(1).InitialAmount = initX * base;
    r{1}.Reactants(2).InitialAmount = initY * base;
    r{3}.Reactants(1).InitialAmount = 10 * base;
    
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
    plot(X(:,1), X(:,2), ':', 'LineWidth', 2.0);
    ylabel('Y'); xlabel('X');
    set(gca, 'LineWidth', 2.0); 
    
    subplot(2, 1, 2);
    box on; hold on;
    plot(t, X(:,1), 'LineWidth', 2.0);
    plot(t, X(:,2), 'LineWidth', 2.0);
    ylabel('Concentration'); xlabel('Time');
    set(gca, 'LineWidth', 2.0); ylim([0 5]);
    legend('X', 'Y')
end