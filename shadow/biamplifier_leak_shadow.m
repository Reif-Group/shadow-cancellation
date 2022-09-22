%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Rajiv Nagipogu Date created: 09/19/2022 Affiliation: Duke
% University %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation of a bi amplifier circuit X + Y --> 2Y


concX = 4; concY = 3; 
N = 3;
alpha = 1e-3;
beta = 1e-8;
% 
ideal_biamp(concX, concY, 100, 1);
ideal_cat_biamp(concX, concY, 100, N, 2);
dna_cat_biamp(concX, concY, 1e7, N, 3); 

function ideal_biamp(concX, concY, stopTime, figNo)
   
    
    % regular rate of each reaction
    rate = 1.0/60;
     
    % tolerance allowed
    absTol = 1e-12;
    relTol = 1e-12;

    % Create the amplification model
    model = sbiomodel('Ideal Biamplification CRN');

    % Define the variables such as reactions, kinetics, and rate constants
    r = cell(1, 1)
    k = cell(1, length(r))
    p = cell(1, length(r))

    % Set the reactions
    r{1} = addreaction(model, 'X + Y -> Y + Y');

    % Set the kinetics
    k{1} = addkineticlaw(r{1}, 'MassAction');

    % Set the rate parameter
    p{1} = addparameter(k{1}, 'c1', 'Value', rate);

    % Add parameter variable name
    k{1}.ParameterVariableNames = {'c1'};

    % Set initial concentrations
    r{1}.Reactants(1).InitialAmount = concX;
    r{1}.Reactants(2).InitialAmount = concY;
    
    % Display the model
    model

    % Display the Reactions
    model.Reactions

    % Display the Species involved
    model.Species

    % Simulate the plot for the ODE
    cs = getconfigset(model, 'active');
    set(cs.SolverOptions, 'AbsoluteTolerance', absTol);
    set(cs.SolverOptions, 'RelativeTolerance', relTol);
    cs.SolverType = 'ode15s';
    cs.StopTime = stopTime;
    [t, X, names] = sbiosimulate(model)

    % Plot the species concentrations w.r.t time
    figure(figNo);
    box on; hold on;
    plot(t, X(:, 1), '-', 'LineWidth', 2.0);
    plot(t, X(:, 2), '-', 'LineWidth', 2.0);
    xlabel('time');
    ylabel('Concentration');
end

function ideal_cat_biamp(concX, concY, stopTime, N, figNo)
    
    % Rate of the reactions
    rate = 1.0/60;

    % Base concentrations
    base = 1.0;

    % tolerance allowed
    absTol = 1e-12;
    relTol = 1e-12;


    % Reactions, Kinetics, and Rate parameters
    r = cell(1, N);
    k = cell(1, length(r));
    p = cell(1, length(r));

    % Define the model
    model = sbiomodel("Ideal Catalytic Biamplifier");

    % Add the reactions
    for i = 1:N
        r{i} = addreaction(model, strcat('x', int2str(i), ' + y', int2str(i), ' -> y', int2str(i), ' + y', int2str(max(1, rem(i+1, N+1)))));
    end

    % Add kinetics
    for i = 1:N
        k{i} = addkineticlaw(r{i}, 'MassAction');
        k{i}.ParameterVariableNames = {strcat('c', int2str(i))};
        p{i} = addparameter(k{i}, strcat('c', int2str(i)), 'Value', rate);
    end
   

    % Set the initial concentrations
    for i = 1:N
        r{i}.Reactants(1).InitialAmount = concX*base;
        r{i}.Reactants(2).InitialAmount = concY*base;
    end

    % Display the model
    model
    model.Reactions
    model.Species

    % Simulate the plot for the ODE
    cs = getconfigset(model, 'active');
    set(cs.SolverOptions, 'AbsoluteTolerance', absTol);
    set(cs.SolverOptions, 'RelativeTolerance', relTol);
    cs.SolverType = 'ode15s';
    cs.StopTime = stopTime;
    [t, X, names] = sbiosimulate(model)
    x = zeros(length(X));
    y = zeros(length(X));
    for i  = 1:length(names)
        if strncmpi(names{i}, 'x', 1)
            x = x + X(:, i);
        end 
        if strncmpi(names{i}, 'y', 1)
            y = y + X(:, i);
        end
        
    end
    x = x/N; 
    y = y/N;

    % Plot the species concentrations w.r.t time
    figure(figNo);
    box on; hold on;
    
    plot(t, x, ':', 'LineWidth', 3.0); % Plots X
    plot(t, y, ':', 'LineWidth', 3.0); % Plots Y
    xlabel('time');
    ylabel('Concentration');
    
end


function dna_cat_biamp(initX, initY, stopTime, N, figNo)
    
    % Base concentration of the reaction
    base = 1;

    % slow down all reactions by multiplying w/ a scaling factor of 1e-3
    alpha = 1e-3;
    % Scale to the nM.
    beta = 1e-8; 

    rate = (1.0/60)*alpha;

    % Fastest reactions
    infRate = 1e6*alpha;

    % Linker reactions
    gamma = 2;
    
    % 
    infConc = 6e1*base*beta;

    % allowed error rate tolerance value in pico molars
    absTol = 1e-12;
    relTol = 1e-12;

    % create a sym bio model for simulation
    model = sbiomodel('Approximate DNA cat biamp');

    r = cell(1, 4*N);
    k = cell(1, length(r));
    p = cell(1, length(r));

    % Add two-step catalytic approximations for X + Y -> 2Y
    for i = 1:N
        r{2*i-1} = addreaction(model, strcat('x', int2str(i), ' + Yg', ...
            int2str(i), ' -> Ixy',int2str(i)));
        r{2*i} = addreaction(model, strcat('Ixy', int2str(i), ...
            ' + Gy', int2str(i), ' -> y', int2str(i), ...
            ' + y', int2str(max(1, rem(i+1, N+1)))));
    end

    for i = 1:N
        r{2*N + 2*i-1} = addreaction(model, strcat('y', int2str(i), ' -> Yg', int2str(i)));
        r{2*N + 2*i} = addreaction(model, strcat('Yg', int2str(i), ' -> y', int2str(i)));
    end

    for i = 1:4*N
        k{i} = addkineticlaw(r{i}, 'MassAction');
        k{i}.ParameterVariableNames = {strcat('c', int2str(i))};
    end
    
    % Assign Rates
    q1 = gamma*N*rate/infConc;
    
    for i = 1:N
        p{2*i-1} = addparameter(k{2*i-1}, strcat('c', int2str(2*i-1)), 'Value', q1);
        p{2*i} = addparameter(k{2*i}, strcat('c', int2str(2*i)), 'Value', infRate);
        p{2*N + 2*i-1} = addparameter(k{2*N + 2*i-1}, strcat('c', int2str(2*N + 2*i-1)), 'Value', infRate);
        p{2*N + 2*i} = addparameter(k{2*N + 2*i}, strcat('c', int2str(2*N + 2*i)), 'Value', infRate);
    end
    
    % Set initial concentrations
    for i = 1:N
        r{2*i-1}.Reactants(1).InitialAmount = initX*(1/N)*beta;
        r{2*i}.Reactants(2).InitialAmount = infConc;
        r{2*N + 2*i-1}.Reactants(1).InitialAmount = initY*(1/N)*beta;
    end

    % Models
    model
    model.Reactions
    model.Species
    
    % get solver config and set ode type, tolerance, stop time
    cs = getconfigset(model,'active');
    set(cs.SolverOptions, 'AbsoluteTolerance', absTol);
    set(cs.SolverOptions, 'RelativeTolerance', relTol);
    cs.SolverType = 'ode15s';
    cs.StopTime = stopTime;
    [t, X, names] = sbiosimulate(model);
    
    x = zeros(length(X));
    y = zeros(length(X));
    yg = zeros(length(X));
    for i  = 1:length(names)
        if strncmpi(names{i}, 'x', 1)
            x = x + X(:, i);
        end
        if strncmpi(names{i}, 'Yg', 2)
            yg = yg + X(:, i);
            continue;
        end
        if strncmpi(names{i}, 'y', 1)
            y = y + X(:, i);
        end
    end
    x = x/beta;
    y = (y + yg)/beta;

    % Set figure
    figure(figNo);
    box on; hold on;
    plot(round(t, 3), x, ':', 'LineWidth', 2.0); % Plots X
    plot(round(t, 3), y, ':', 'LineWidth', 2.0); % Plots Y
    xlabel('time');
    ylabel('Concentration');

end






















































