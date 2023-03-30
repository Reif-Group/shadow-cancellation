%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Shalin Shah
% Date created: 04/01/2020
% Affiliation: Duke University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cat_consensus(2.7, 0.3, 1, 5, 5, 1);
% cat_consensus(0.6, 2.4, 1, 5, 5, 2);
% cat_consensus(1.3, 1.7, 1, 5, 5, 3);
% cat_consensus(1.6, 1.4, 1, 5, 5, 4);
ideal_consensus(27, 3, 30, 1, 1);

function ideal_consensus(initA, initB, initY, stopTime, figNo)
    % enter A and B concentration (in nM) assuming all gates are several
    % micro molars
    %
    % stopTime is the total simulation time
    %
    % figNo is the display figure number

    rate = 1;
    % Create the amplification model
    model = sbiomodel('Catalytic CRNs for consensus protocol');
    r = cell(1, 3);
    k = cell(1, length(r));
    p = cell(1, length(r));
    
    % Add reaction set to the solver
    r{1} = addreaction(model, 'A + B -> Y + Y'); 
    r{2} = addreaction(model, 'A + Y -> A + A'); 
    r{3} = addreaction(model, 'B + Y -> B + B'); 
    
    % Set the Kinetic Law for Reactions
    for i = 1 : 3
        k{i} = addkineticlaw(r{i}, 'MassAction'); 
        k{i}.ParameterVariableNames = {strcat('c', int2str(i))};
        p{i} = addparameter(k{i}, strcat('c', int2str(i)), 'Value', rate);
    end

    % Set initial concentrations
    r{1}.Reactants(1).InitialAmount = initA;
    r{1}.Reactants(2).InitialAmount = initB;
    r{2}.Reactants(2).InitialAmount = initY;
    
    % Display the Completed Model Objects
    model

    % Display the Reaction Objects
    model.Reactions

    % Display the Species Objects
    model.Species

    % Simulate ODE and plot for 20 minutes
    cs = getconfigset(model,'active');
    cs.SolverType = 'ode15s';
    cs.StopTime = stopTime;
    [t,X] = sbiosimulate(model);

    % plot X over the time
    figure(1); subplot(2, 2, figNo); hold on; box on;
    plot(t, X(:, 1), 'LineWidth', 2.0);
    plot(t, X(:, 2), 'LineWidth', 2.0);
    plot(t, X(:, 3), ':', 'LineWidth', 2.0);
    ylabel('concentration'); xlabel('time');
    set(gca, 'LineWidth', 2.0); legend('A', 'B', 'Y')
end

function cat_consensus(initA, initB, initY, N, stopTime, figNo)
    % enter A and B concentration (in nM) assuming all gates are several
    % micro molars
    %
    % stopTime is the total simulation time
    %
    % figNo is the display figure number

    % set all the default params and initialize model
    % default rate and concentration
    rate = 1;
    base = 1;
    
    % allowed error rate tolerance value in pico molars
    absTol = 1e-12;
    relTol = 1e-12;
    
    % Create the amplification model
    model = sbiomodel('Catalytic CRNs for consensus protocol');
    r = cell(1, 3*N);
    k = cell(1, length(r));
    p = cell(1, length(r));
    
    % Add reaction set A+B -> 2Y to the solver
    for i = 1 : N - 1
        r{i} = addreaction(model, strcat('a', int2str(i),' + b', int2str(i),...
            ' -> y', int2str(i),' + y', int2str(i+1))); 
    end
    r{N} = addreaction(model, strcat('a', int2str(N),' + b', int2str(N),...
            ' -> y', int2str(N),' + y1')); 
    
    % Add reaction set Y+A -> 2A to the solver
    for i = 1 : N - 1
        r{N+i} = addreaction(model, strcat('y', int2str(i),' + a', int2str(i),...
            ' -> a', int2str(i),' + a', int2str(i+1))); 
    end
    r{2*N} = addreaction(model, strcat('y', int2str(N),' + a', int2str(N),...
            ' -> a', int2str(N),' + a1')); 

    % Add reaction set Y+B -> 2B to the solver
    for i = 1 : N - 1
        r{2*N+i} = addreaction(model, strcat('y', int2str(i),' + b', int2str(i),...
            ' -> b', int2str(i),' + b', int2str(i+1))); 
    end
    r{3*N} = addreaction(model, strcat('y', int2str(N),' + b', int2str(N),...
            ' -> b', int2str(N),' + b1')); 
    
    % Set the Kinetic Law for Reactions
    for i = 1 : 3*N
        k{i} = addkineticlaw(r{i}, 'MassAction'); 
        k{i}.ParameterVariableNames = {strcat('c', int2str(i))};
        p{i} = addparameter(k{i}, strcat('c', int2str(i)), 'Value', rate);
    end

    % Set initial concentrations
    for i = 1:N
        r{i}.Reactants(1).InitialAmount = initA * base * 0;
    end
    for i = 1:N
        r{i}.Reactants(2).InitialAmount = initB * base;
        r{N+i}.Reactants(1).InitialAmount = initY * base;
    end

    % Display the Completed Model Objects
    model

    % Display the Reaction Objects
    model.Reactions

    % Display the Species Objects
    model.Species

    % Simulate ODE and plot for 20 minutes
    cs = getconfigset(model,'active');
    set(cs.SolverOptions, 'AbsoluteTolerance', absTol);
    set(cs.SolverOptions, 'RelativeTolerance', relTol);
    cs.SolverType = 'ode15s';
    cs.StopTime = stopTime;
    [t,X] = sbiosimulate(model);

    % plot A, B and Y over the time
    figure(2); subplot(2, 2, figNo); hold on; box on;
    a = X(:, 1) + sum(X(:, 5:3:3*N), 2);
    plot(t, a./N, 'LineWidth', 2.0);
    b = X(:, 2) + sum(X(:, 6:3:3*N), 2);
    plot(t, b./N, 'LineWidth', 2.0);
    y = sum(X(:, 3:4), 2) + sum(X(:, 7:3:3*N), 2);
    plot(t, y./N, ':', 'LineWidth', 2.0);
    ylabel('concentration'); xlabel('time');
    set(gca, 'LineWidth', 2.0); legend('A', 'B', 'Y')
end