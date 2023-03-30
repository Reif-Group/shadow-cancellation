%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Shalin Shah
% Date created: 04/01/2020
% Affiliation: Duke University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
concX = 1; concY = 0.5; stopTime = 2;

idealampnoncat(concX, concY, stopTime, 1)
idealampcat(concX, concY, stopTime, 1)
idealampautocat(concX, concY, stopTime, 1)

function idealampnoncat(concA, concB, stopTime, figNo)
    % enter A concentration assuming all gates are excess
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

    % Add reaction set X -> 2X to the solver
    r1 = addreaction(model,'Y -> Z');
        
    % Set the Kinetic Law for Reactions.
    k1 = addkineticlaw(r1, 'MassAction'); 
    
    p1 = addparameter(k1, 'c1', 'Value', rate);
    
    % Set the Kinetic Law for Reactions.
    k1.ParameterVariableNames = {'c1'};
    
    % Set initial amounts for species 
    r1.Reactants(1).InitialAmount = concA;        % A
%     r1.Reactants(2).InitialAmount = concB;        % B

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
    plot(t, X(:, 2), 'LineWidth', 2.0);
    ylabel('Concentration'); xlabel('Time');
    set(gca, 'LineWidth', 2.0); 
end

function idealampcat(concA, concB, stopTime, figNo)
    % enter A concentration assuming all gates are excess
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

    % Add reaction set X -> 2X to the solver
    r1 = addreaction(model,'Y -> Z + Y');
        
    % Set the Kinetic Law for Reactions.
    k1 = addkineticlaw(r1, 'MassAction'); 
    
    p1 = addparameter(k1, 'c1', 'Value', rate);
    
    % Set the Kinetic Law for Reactions.
    k1.ParameterVariableNames = {'c1'};
    
    % Set initial amounts for species 
    r1.Reactants(1).InitialAmount = concA;        % A
%     r1.Reactants(2).InitialAmount = concB;        % A

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
    plot(t, X(:, 2), 'LineWidth', 2.0);
    ylabel('Concentration'); xlabel('Time');
    set(gca, 'LineWidth', 2.0); 
end

function idealampautocat(concA, concB, stopTime, figNo)
    % enter A concentration assuming all gates are excess
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

    % Add reaction set X -> 2X to the solver
    r1 = addreaction(model,'Y -> Z + Y + Y');
        
    % Set the Kinetic Law for Reactions.
    k1 = addkineticlaw(r1, 'MassAction'); 
    
    p1 = addparameter(k1, 'c1', 'Value', rate);
    
    % Set the Kinetic Law for Reactions.
    k1.ParameterVariableNames = {'c1'};
    
    % Set initial amounts for species 
    r1.Reactants(1).InitialAmount = concA;        % A
%     r1.Reactants(2).InitialAmount = concB;        % A

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
    plot(t, X(:, 2), 'LineWidth', 2.0);
    ylabel('Concentration'); xlabel('Time');
    set(gca, 'LineWidth', 2.0); 
    
    plot([0:0.01:stopTime], [0:0.01:stopTime], ':', 'LineWidth', 2)
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel', a, 'fontsize', 16)

end
