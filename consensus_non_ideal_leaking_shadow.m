%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Shalin Shah
% Date created: 04/01/2020
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
dna_cat_consensus(2.7, 0.3, 3, 3, 50000, 1);

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
    r{1} = addreaction(model, 'A + B -> 2 Y');
    r{2} = addreaction(model, 'Y + A -> 2 A');
    r{3} = addreaction(model, 'Y + B -> 2 B');

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

function ideal_cat_consensus(initA, initB, initY, N, stopTime, figNo)
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
    r = cell(1, 4*N);
    k = cell(1, length(r));
    p = cell(1, length(r));

    % Add reactions
    for i = 1:N
        i_nxt = max(1, rem(i+1, N+1));
        r{i} = addreaction(model, strcat('a', int2str(i), ' + b', ...
            int2str(i), ' -> y', int2str(i), ' + y', int2str(i_nxt)));
        r{N + i} = addreaction(model, strcat('y', int2str(i), ' + a', ...
            int2str(i), ' -> a', int2str(i), ' + a', int2str(i_nxt)));
        r{2*N + i} = addreaction(model, strcat('y', int2str(i), ' + b', ...
            int2str(i), ' -> b', int2str(i), ' + b', int2str(i_nxt)));
    end

    % Add kinetic laws and parameters
    for i = 1:3*N
        k{i} = addkineticlaw(r{i}, 'MassAction');
        k{i}.ParameterVariableNames = {strcat('c', int2str(i))};
        p{i} = addparameter(k{i}, strcat('c', int2str(i)), 'Value', rate);
    end
    
    % Add initial concentrations
    for i = 1:N
        i_nxt = max(1, rem(i+1, N+1));
        r{i}.Reactants(1).InitialAmount = initA*base;
        r{i}.Reactants(2).InitialAmount = initB*base;
        r{N+i}.Reactants(1).InitialAmount = initY*base; 
    end

    % Model display
    model
    model.Species

    % get solver config and set ode type, tolerance, stop time
    cs = getconfigset(model,'active');
    set(cs.SolverOptions, 'AbsoluteTolerance', absTol);
    set(cs.SolverOptions, 'RelativeTolerance', relTol);
    cs.SolverType = 'ode15s';
    cs.StopTime = stopTime;
    [t, X, names] = sbiosimulate(model);
    
    a = zeros(length(X));
    b = zeros(length(X));
    y = zeros(length(X));

    for i  = 1:length(names)
        if strncmpi(names{i}, 'a', 1)
            disp(strcat(names{i}, 'should be a'));
            a = a + X(:, i);
        end
        if strncmpi(names{i}, 'b', 1)
            disp(strcat(names{i}, 'should be b'));
            b = b + X(:, i);
        end
        if strncmpi(names{i}, 'y', 1)
            disp(strcat(names{i}, 'should be y'));
            y = y + X(:, i);
        end
    end
    disp(names);
    a = a/N;
    b = b/N;
    y = y/N;

    % Set figure
    figure(figNo);
    box on; hold on;
    plot(t, a, 'LineWidth', 2.0, 'Color', 'blue');
    plot(t, b, 'LineWidth', 2.0, 'Color', 'red');
    plot(t, y, ':', 'LineWidth', 2.0, 'Color', 'green');
     
    ylabel('concentration'); 
    xlabel('time');
    set(gca, 'LineWidth', 2.0); 

end

function dna_cat_consensus(initA, initB, initY, N, stopTime, figNo)
 
    % Base concentration
    base = 1.0;
    
    % Rate scaling
    alpha = 1e-3;

    % Concentration scaling 
    beta = 1e-8;

    % Set tolerance levels
    absTol = 1e-12;
    relTol = 1e-12;

    % Flags for leak and shadow
    leak = 1;
    shadow = 0;
    
    % Rate
    rateConstant = 1.0;
    rate = rateConstant*alpha;

    % Set infConc and infRate
    infConc = 1e4*beta;
    infRate = 1e4*alpha;

    % Set leak rate
    leak_rate = 1e-6*alpha;

    % Define the model
    model = sbiomodel('Ideal Consensus');
    
    % Cells
    r = cell(1, 12*N);
    k = cell(1, length(r));
    p = cell(1, length(r));
    
    % Regular circuit
    for i = 1:N
        istr = int2str(i);
        inxt = int2str(max(1, rem(i+1, N+1)));
        r{2*i-1} = addreaction(model, strcat('a', istr, ' + Bg', istr, ' -> Iy', istr));
        r{2*i} = addreaction(model, strcat('Iy', istr, ' + Gy', istr, ' -> y', istr, ' + y', inxt));
        
        r{2*N + 2*i - 1} = addreaction(model, strcat('b', istr, ' + Yg', istr, ' -> Ib', istr));
        r{2*N + 2*i} = addreaction(model, strcat('Ib', istr, ' + Gb', istr, ' -> b', istr, ' + b', inxt));
        
        r{4*N + 2*i - 1} = addreaction(model, strcat('y', istr, ' + Ag', istr, ' -> Ia', istr));
        r{4*N + 2*i} = addreaction(model, strcat('Ia', istr, ' + Ga', istr, ' -> a', istr, ' + a', inxt));
    end
    
    % Linker reactions
    for i = 1:N
        istr = int2str(i);
        inxt = int2str(max(1, rem(i+1, N+1)));
        r{6*N + 2*i - 1} = addreaction(model, strcat('a', istr, ' + La', istr, ' -> Ag', istr));
        r{6*N + 2*i} = addreaction(model, strcat('Ag', istr, ' + RLa', istr, ' -> a', istr));
        
        r{8*N + 2*i - 1} = addreaction(model, strcat('b', istr, ' + Lb', istr, ' -> Bg', istr));
        r{8*N + 2*i} = addreaction(model, strcat('Bg', istr, ' + RLb', istr, ' -> b', istr));
        
        r{10*N + 2*i - 1} = addreaction(model, strcat('y', istr, ' + Ly', istr, ' -> Yg', istr));
        r{10*N + 2*i} = addreaction(model, strcat('Yg', istr, ' + RLy', istr, ' -> y', istr));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% Leak Profile %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Leaks in regular reactions
    for i = 1:N
        istr = int2str(i);
        inxt = int2str(max(1, rem(i+1, N+1)));
        r{12*N + 2*i-1} = addreaction(model, strcat('Bg', istr, ' -> Iy', istr));
        r{12*N + 2*i} = addreaction(model, strcat('Gy', istr, ' -> y', istr, ' + y', inxt));
        
        r{14*N + 2*i - 1} = addreaction(model, strcat('Yg', istr, ' -> Ib', istr));
        r{14*N + 2*i} = addreaction(model, strcat('Gb', istr, ' -> b', istr, ' + b', inxt));
        
        r{16*N + 2*i - 1} = addreaction(model, strcat('Ag', istr, ' -> Ia', istr));
        r{16*N + 2*i} = addreaction(model, strcat('Ga', istr, ' -> a', istr, ' + a', inxt));
    end

    % Leaks in linker reactions
    for i = 1:N
        istr = int2str(i);
        inxt = int2str(max(1, rem(i+1, N+1)));
        r{18*N + 2*i - 1} = addreaction(model, strcat('La', istr, ' -> Ag', istr));
        r{18*N + 2*i} = addreaction(model, strcat('RLa', istr, ' -> a', istr));
        
        r{20*N + 2*i - 1} = addreaction(model, strcat('Lb', istr, ' -> Bg', istr));
        r{20*N + 2*i} = addreaction(model, strcat('RLb', istr, ' -> b', istr));
        
        r{22*N + 2*i - 1} = addreaction(model, strcat('Ly', istr, ' -> Yg', istr));
        r{22*N + 2*i} = addreaction(model, strcat('RLy', istr, ' -> y', istr));
    end


    % Add rate kinetics
    for i = 1:24*N
        k{i} = addkineticlaw(r{i}, 'MassAction');
        k{i}.ParameterVariableNames = {strcat('c', int2str(i))};
    end
    
    % Since the 
    q1 = 2*N*rate*(1/beta);
    for i = 1:N
        p{2*i-1} = addparameter(k{2*i-1}, strcat('c', int2str(2*i-1)), 'Value', q1);
        p{2*i} = addparameter(k{2*i}, strcat('c', int2str(2*i)), 'Value', infRate);

        p{2*N + 2*i-1} = addparameter(k{2*N + 2*i-1}, strcat('c', int2str(2*N + 2*i-1)), 'Value', q1);
        p{2*N + 2*i} = addparameter(k{2*N + 2*i}, strcat('c', int2str(2*N + 2*i)), 'Value', infRate);

        p{4*N + 2*i-1} = addparameter(k{4*N + 2*i-1}, strcat('c', int2str(4*N + 2*i-1)), 'Value', q1);
        p{4*N + 2*i} = addparameter(k{4*N + 2*i}, strcat('c', int2str(4*N + 2*i)), 'Value', infRate);
        
        % Rate constants of the linker reactions. Since the gates have
        % infinite concentration, I think setting their rate constant to
        % regular rates is fine.  Note: These rate parameters shouldn't
        % matter.
        p{6*N + 2*i-1} = addparameter(k{6*N + 2*i-1}, strcat('c', int2str(6*N + 2*i-1)), 'Value', infRate);
        p{6*N + 2*i} = addparameter(k{6*N + 2*i}, strcat('c', int2str(6*N + 2*i)), 'Value', infRate);

        p{8*N + 2*i-1} = addparameter(k{8*N + 2*i-1}, strcat('c', int2str(8*N + 2*i-1)), 'Value', infRate);
        p{8*N + 2*i} = addparameter(k{8*N + 2*i}, strcat('c', int2str(8*N + 2*i)), 'Value', infRate);

        p{10*N + 2*i-1} = addparameter(k{10*N + 2*i-1}, strcat('c', int2str(10*N + 2*i-1)), 'Value', infRate);
        p{10*N + 2*i} = addparameter(k{10*N + 2*i}, strcat('c', int2str(10*N + 2*i)), 'Value', infRate);
        
        p{12*N + 2*i-1} = addparameter(k{12*N + 2*i-1}, strcat('c', int2str(12*N + 2*i-1)), 'Value', leak*leak_rate);
        p{12*N + 2*i} = addparameter(k{12*N + 2*i}, strcat('c', int2str(12*N + 2*i)), 'Value', leak*leak_rate);
        
        p{14*N + 2*i-1} = addparameter(k{14*N + 2*i-1}, strcat('c', int2str(14*N + 2*i-1)), 'Value', leak*leak_rate);
        p{14*N + 2*i} = addparameter(k{14*N + 2*i}, strcat('c', int2str(14*N + 2*i)), 'Value', leak*leak_rate);

        p{16*N + 2*i-1} = addparameter(k{16*N + 2*i-1}, strcat('c', int2str(16*N + 2*i-1)), 'Value', leak*leak_rate);
        p{16*N + 2*i} = addparameter(k{16*N + 2*i}, strcat('c', int2str(16*N + 2*i)), 'Value', leak*leak_rate);

        p{18*N + 2*i-1} = addparameter(k{18*N + 2*i-1}, strcat('c', int2str(18*N + 2*i-1)), 'Value', leak*leak_rate);
        p{18*N + 2*i} = addparameter(k{18*N + 2*i}, strcat('c', int2str(18*N + 2*i)), 'Value', leak*leak_rate);

        p{20*N + 2*i-1} = addparameter(k{20*N + 2*i-1}, strcat('c', int2str(20*N + 2*i-1)), 'Value', leak*leak_rate);
        p{20*N + 2*i} = addparameter(k{20*N + 2*i}, strcat('c', int2str(20*N + 2*i)), 'Value', leak*leak_rate);

        p{22*N + 2*i-1} = addparameter(k{22*N + 2*i-1}, strcat('c', int2str(22*N + 2*i-1)), 'Value', leak*leak_rate);
        p{22*N + 2*i} = addparameter(k{22*N + 2*i}, strcat('c', int2str(22*N + 2*i)), 'Value', leak*leak_rate);

    end

    % Add initial amounts
    for i = 1:N
        % Substrates
        r{2*i-1}.Reactants(1).InitialAmount = beta*initA/N;
        r{2*N + 2*i-1}.Reactants(1).InitialAmount = beta*initB/N;
        r{4*N + 2*i-1}.Reactants(1).InitialAmount = beta*initY/N;

        % Regular gates
        r{2*i}.Reactants(2).InitialAmount = infConc;
        r{2*N + 2*i}.Reactants(2).InitialAmount = infConc;
        r{4*N + 2*i}.Reactants(2).InitialAmount = infConc;

        % Linker gates
        r{6*N + 2*i-1}.Reactants(2).InitialAmount = infConc;
        r{6*N + 2*i}.Reactants(2).InitialAmount = infConc;
        r{8*N + 2*i-1}.Reactants(2).InitialAmount = infConc;
        r{8*N + 2*i}.Reactants(2).InitialAmount = infConc;
        r{10*N + 2*i-1}.Reactants(2).InitialAmount = infConc;
        r{10*N + 2*i}.Reactants(2).InitialAmount = infConc;
    end

    % Model display
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
    
    a = zeros(length(t));
    b = zeros(length(t));
    y = zeros(length(t));
    ag = zeros(length(t));
    bg = zeros(length(t));
    yg = zeros(length(t));

    for i  = 1:length(names)
        if strncmpi(names{i}, 'Ag', 2)
            disp(strcat(names{i}, 'should be Ag'));
            ag = ag + X(:, i);
            continue;
        end
        if strncmpi(names{i}, 'Bg', 2)
            disp(strcat(names{i}, 'should be Bg'));
            bg = bg + X(:, i);
            continue;
        end
        if strncmpi(names{i}, 'Yg', 2)
            disp(strcat(names{i}, 'should be Yg'));
            yg = yg + X(:, i);
            continue;
        end
        if strncmpi(names{i}, 'a', 1)
            disp(strcat(names{i}, 'should be a'));
            a = a + X(:, i);
            continue;
        end
        if strncmpi(names{i}, 'b', 1)
            disp(strcat(names{i}, 'should be b'));
            b = b + X(:, i);
            continue;
        end
        if strncmpi(names{i}, 'y', 1)
            disp(strcat(names{i}, 'should be y'));
            y = y + X(:, i);
            continue;
        end
    end
    
    atotal = (a + ag);
    btotal = (b + bg);
    ytotal = (y + yg);
    
    % Set figure
    figure(figNo);
    box on; hold on;
%     plot(t/3600, a/1.e-9, 'LineWidth', 2.0, 'Color', 'blue');
%     plot(t/3600, b/1.e-9, 'LineWidth', 2.0, 'Color', 'red');
%     plot(t/3600, y/1.e-9, ':', 'LineWidth', 2.0, 'Color', 'green');
%     plot(t/3600, ag/1.e-9, 'LineWidth', 2.0, 'Color', 'blue');
%     plot(t/3600, bg/1.e-9, 'LineWidth', 2.0, 'Color', 'red');
%     plot(t/3600, yg/1.e-9, ':', 'LineWidth', 2.0, 'Color', 'green');
    plot(t/3600, atotal/1.e-9, 'LineWidth', 2.0, 'Color', 'blue');
    plot(t/3600, btotal/1.e-9, 'LineWidth', 2.0, 'Color', 'red');
    plot(t/3600, ytotal/1.e-9, ':', 'LineWidth', 2.0, 'Color', 'green');
    ylabel('concentration'); 
    xlabel('time');



end
