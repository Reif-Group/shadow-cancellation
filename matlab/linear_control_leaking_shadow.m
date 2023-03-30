stopTime = 200;

% demonstrates consensus reactions reach majority
% ideal_rps(2, 10, 13, stopTime, 1)
dna_cat_rps(2, 10, 13, 3, stopTime, 1);

function ideal_linear_control(init, init, initC, stopTime, figNo)
    
    % Rate
    rate = 1/60;

    % Base concentration
    base = 1.0;
    
    % Set tolerance levels
    absTol = 1e-12;
    relTol = 1e-12;
    % Define the model
    model = sbiomodel('Ideal RPS');

    % Cells
    r = cell(1, 5);
    k = cell(1, length(r));
    p = cell(1, length(r));

    % Add reactions
    r{1} = addreaction(model, 'A -> B'); % k = 1 
    r{2} = addreaction(model, 'B -> I'); % k = 3
    r{3} = addreaction(model, 'I -> B'); % infRate
    r{4} = addreaction(model, 'B -> A'); % k = -3
    r{5} = addreaction(model, 'Buff -> A'); % k = sin(t)/Conc(Buff)

    % Add kinetic laws
    for i = 1:5
        k{i} = addkineticlaw(r{i}, 'MassAction');
        k{i}.ParameterVariableNames = {strcat('c', int2str(i))};
    end

    % Add rate parameters
    p{1} = addparameter(k{1}, 'c1', 'Value', 1);
    p{2} = addparameter(k{1}, 'c2', 'Value', 3);
    p{3} = addparameter(k{1}, 'c3', 'Value', infRate);
    p{4} = addparameter(k{1}, 'c4', 'Value', 1);
    p{5} = addparameter(k{1}, 'c5', 'Value', 1);
    

    % Initial concentrations
    r{1}.Reactants(1).InitialAmount = initA*base;
    r{1}.Reactants(2).InitialAmount = initB*base;
    r{2}.Reactants(1).InitialAmount = initC*base;

    % Display model
    model
    model.Reactions
    model.Species

    % ODE Solve
    cs = getconfigset(model,'active');
    set(cs.SolverOptions, 'AbsoluteTolerance', absTol);
    set(cs.SolverOptions, 'RelativeTolerance', relTol);
    cs.SolverType = 'ode15s';
    cs.StopTime = stopTime;
    [t,X, names] = sbiosimulate(model);
    disp(names);

    % Plots
    figure(figNo);
    hold on; box on;
    a = X(:, 1);
    b = X(:, 2);
    c = X(:, 3);


    plot(t, a, 'LineWidth', 2.0);
    plot(t, b, ':', 'LineWidth', 2.0);
    plot(t, c, ':', 'LineWidth', 2.0);
    ylabel('concentration'); 
    xlabel('time');
    set(gca, 'LineWidth', 2.0); legend('A', 'B', 'C')

end
