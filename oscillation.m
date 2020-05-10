clear; clc;
concX = 4; concY = 2; stopTime = 10; nCatalytic = 20;

% compare ideal lotka-volterra simulation with approximation
ideallvoscillator(concX, concY, stopTime, 1)
approxlvoscillator(concX, concY, nCatalytic, stopTime, 2);

function ideallvoscillator(concX, concY, stopTime, figNo)
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
    r{2} = addreaction(model,'Y -> phi');
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
    r{1}.Reactants(1).InitialAmount = concX;        % A
    r{2}.Reactants(1).InitialAmount = concY;        % B

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
    r = cell(1, N+4);
    k = cell(1, length(r));
    p = cell(1, length(r));
    
    % add two ideal reactions for LV
    r{1} = addreaction(model,'X + Y -> Y + Y');
    r{2} = addreaction(model,'Y -> waste'); 
    
    % add catalytic approximations for last step
    r{3} = addreaction(model, 'X -> X + x1');
    for i = 4 : N+3
        r{i} = addreaction(model, strcat('x', int2str(i-3), ' -> x', int2str(i-3),' + x', int2str(i-2))); 
    end
    % add the final reaction to release X
    r{N+4} = addreaction(model, strcat('x', int2str(N+1), ' -> x', int2str(N+1),' + X'));
   
    % add all the rate laws as mass action kinetics
    for i = 1:N+4
        k{i} = addkineticlaw(r{i}, 'MassAction'); 
        k{i}.ParameterVariableNames = {strcat('c', int2str(i))};
        p{i} = addparameter(k{i}, strcat('c', int2str(i)), 'Value', rate);
    end
    
    % set initial concentrations
    r{1}.Reactants(1).InitialAmount = initX * base;
    r{1}.Reactants(2).InitialAmount = initY * base;
    
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
    set(gca, 'LineWidth', 2.0); %ylim([0 20]);
    legend('X', 'Y')
end