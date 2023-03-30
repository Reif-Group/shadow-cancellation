%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Rajiv Nagipogu
% Affiliation: Duke University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ideal_cat_consensus(2.7, 0.3, 1, 3, 200, 1);
% ideal_cat_consensus(0.6, 2.4, 1, 3, 200, 2);
% ideal_cat_consensus(1.3, 1.7, 1, 3, 200, 3);
% ideal_cat_consensus(1.6, 1.4, 1, 3, 200, 4);
% ideal_consensus(27, 3, 30, 30, 1);
% ideal_consensus(24, 6, 30., 30, 2);
% ideal_consensus(9, 21, 30., 30, 3);
% ideal_consensus(6, 24, 30., 30, 4);

% demonstrates consensus reactions reach majority
dna_cat_consensus(2.7, 0.3, 1, 3, 20000, 1);

function ideal_consensus(initA, initB, initY, stopTime, figNo)
    
    % Rate
    rate = 1/60;

    % Base concentration
    base = 1.0;
    
    % Set tolerance levels
    absTol = 1e-12;
    relTol = 1e-12;
    % Define the model
    model = sbiomodel('Ideal Consensus');

    % Cells
    r = cell(1, 3);
    k = cell(1, length(r));
    p = cell(1, length(r));

    % Add reactions
    
    % Add kinetic laws
    for i = 1:3
        k{i} = addkineticlaw(r{i}, 'MassAction');
        k{i}.ParameterVariableNames = {strcat('c', int2str(i))};
        p{i} = addparameter(k{i}, strcat('c', int2str(i)), 'Value', rate);
    end

    % Initial concentrations
    r{1}.Reactants(1).InitialAmount = initA*base;
    r{1}.Reactants(2).InitialAmount = initB*base;
    r{2}.Reactants(1).InitialAmount = initY*base;

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
    y = X(:, 3);


    plot(t, a, 'LineWidth', 2.0);
    plot(t, b, 'LineWidth', 2.0);
    plot(t, y, ':', 'LineWidth', 2.0);
    ylabel('concentration'); 
    xlabel('time');
    set(gca, 'LineWidth', 2.0); legend('A', 'B', 'Y')

end
