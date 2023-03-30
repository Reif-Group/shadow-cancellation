% demonstrates consensus reactions reach majority
crnconsensus(2.7, 0.3, 3, 3800, 1);
crnconsensus(2.4, 0.6, 3, 3800, 2);
crnconsensus(2.1, 0.9, 3, 3800, 3);
crnconsensus(0.9, 2.1, 3, 3800, 4);
crnconsensus(0.6, 2.4, 3, 3800, 5);
crnconsensus(0.3, 2.7, 3, 3800, 6);

function crnconsensus(concA, concB, concY, stopTime, figNo)
    % enter A and B concentration (in nM) assuming all gates are several
    % micro molars
    %
    % stopTime is the total simulation time
    %
    % figNo is the display figure number

    k = 1.0;
    alpha = 1e-3; beta = 1e-8; gamma = 2;
    % slow down all reactions by multiplying with a scaling factor 1e-3
    rateConst = k * (alpha);
    % fastest reaction are assumed to occur 
    infRate = 1e3 * (k);
    % inf concentration of gates will be 100 uM
    infConc = 1e3 * (beta);
    % allowed error rate tolerance value in pico molars
    absTol = 1e-8;
    relTol = 1e-8;
    
    % Create the amplification model
    model = sbiomodel('Polymerase-based strand displacement amplification CRN');

    r = cell(1, 12);
    k = cell(1, length(r));
    p = cell(1, length(r));
    
    % Add reaction set A+B -> 2Y to the solver
    r{1} = addreaction(model,'A + Gb -> I1');
    r{2} = addreaction(model,'I1 + Gyo -> Y + Y');
    r{3} = addreaction(model,'B + BfB1 -> Gb');
    r{4} = addreaction(model,'Gb + BfB2 -> B');

    % Add reaction set Y+A -> 2A to the solver
    r{5} = addreaction(model,'Y + Ga -> I2');
    r{6} = addreaction(model,'I2 + Gao -> A + A');
    r{7} = addreaction(model,'A + BfA1 -> Ga');
    r{8} = addreaction(model,'Ga + BfA2 -> A');


    % Add reaction set B+Y -> 2B to the solver
    r{9} = addreaction(model,'B + Gy -> I3');
    r{10} = addreaction(model,'I3 + Gbo -> B + B');
    r{11} = addreaction(model,'Y + BfY1 -> Gy');
    r{12} = addreaction(model,'Gy + BfY2 -> Y');
    
    % Set the Kinetic Law for Reactions
    for i = 1:length(r)
        k{i} = addkineticlaw(r{i}, 'MassAction');
    end
    
    % set the rate constants for each reaction
    p{1} = addparameter(k{1}, 'c1', 'Value', gamma*rateConst*(1/beta));
    p{2} = addparameter(k{2}, 'c2', 'Value', infRate);
    p{3} = addparameter(k{3}, 'c3', 'Value', infRate);
    p{4} = addparameter(k{4}, 'c4', 'Value', infRate);
   
    p{5} = addparameter(k{5}, 'c5', 'Value', gamma*rateConst*(1/beta));
    p{6} = addparameter(k{6}, 'c6', 'Value', infRate);
    p{7} = addparameter(k{7}, 'c7', 'Value', infRate);
    p{8} = addparameter(k{8}, 'c8', 'Value', infRate);
    
    p{9} = addparameter(k{9}, 'c9', 'Value', gamma*rateConst*(1/beta));
    p{10} = addparameter(k{10}, 'c10', 'Value', infRate);
    p{11} = addparameter(k{11}, 'c11', 'Value', infRate);
    p{12} = addparameter(k{12}, 'c12', 'Value', infRate);

    
    % Set the Kinetic Law for Reactions.
    for i = 1:length(r)
        k{i}.ParameterVariableNames = {strcat('c',num2str(i))};
    end

    % Set initial amounts for species 
    r{1}.Reactants(1).InitialAmount = concA*beta;      % A
    r{5}.Reactants(1).InitialAmount = concY*beta;      % B
    r{9}.Reactants(1).InitialAmount = concB*beta;      % Y

    % set gate concentrations in excess
    r{2}.Reactants(2).InitialAmount = infConc;              % Gba
    r{3}.Reactants(2).InitialAmount = infConc;              % Gab
    r{4}.Reactants(2).InitialAmount = infConc;              % Giy
    
    r{6}.Reactants(2).InitialAmount = infConc;              % Gba
    r{7}.Reactants(2).InitialAmount = infConc;              % Gab
    r{8}.Reactants(2).InitialAmount = infConc;              % Giy
    
    r{10}.Reactants(2).InitialAmount = infConc;              % Gba
    r{11}.Reactants(2).InitialAmount = infConc;              % Gab
    r{12}.Reactants(2).InitialAmount = infConc;              % Giy

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

    % plot X over the time
    subplot(2, 3, figNo); hold on; box on;
    plot(t./3600, X(:,1)./1e-9, 'LineWidth', 2.0);
    plot(t./3600, X(:,5)./1e-9, ':', 'LineWidth', 2.0);
    plot(t./3600, X(:,6)./1e-9, 'LineWidth', 2.0);
    ylim([0 35]); xlim([0 stopTime/3600])
    %title('Polymerase-based concensus network');
    ylabel('Concentration (nM)'); xlabel('Time (hours)');
    set(gca, 'LineWidth', 2.0); 
end