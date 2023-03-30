stopTime = 200;

% demonstrates consensus reactions reach majority
% ideal_rps(2, 10, 13, stopTime, 1)
dna_cat_rps(2, 10, 13, 3, stopTime, 1, 1, 0);
% dna_cat_rps(2, 10, 13, 3, stopTime, 2, 1, 1);

function ideal_rps(initA, initB, initC, stopTime, figNo)
    
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
    r = cell(1, 3);
    k = cell(1, length(r));
    p = cell(1, length(r));

    % Add reactions
    r{1} = addreaction(model, 'B + A -> 2 B');
    r{2} = addreaction(model, 'C + B -> 2 C');
    r{3} = addreaction(model, 'A + C -> 2 A');

    % Add kinetic laws
    s = rng;
    rands = rand(1, 3);
    disp(rands(1, 1));
    for i = 1:3
        k{i} = addkineticlaw(r{i}, 'MassAction');
        k{i}.ParameterVariableNames = {strcat('c', int2str(i))};
%         rate_i = rate + (rate/10)*rands(1, i);
        p{i} = addparameter(k{i}, strcat('c', int2str(i)), 'Value', rate);
    end

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

function ideal_cat_rps(initA, initB, initC, N, stopTime, figNo)
    
    % Rate
    rate = 1/60;

    % Base concentration
    base = 1.0;
    
    % Set tolerance levels
    absTol = 1e-12;
    relTol = 1e-12;
    % Define the model
    model = sbiomodel('Catalytic RPS');

    % Cells
    r = cell(1, 3*N);
    k = cell(1, length(r));
    p = cell(1, length(r));

    for i = 1:N
        istr = int2str(i);
        inxt = int2str(max(1, rem(i+1, N+1)));

        r{i} = addreaction(model, strcat('b', istr, ' + a', istr, ' -> b', istr, ' + b', inxt));
        r{N + i} = addreaction(model, strcat('c', istr, ' + b', istr, ' -> c', istr, ' + c', inxt));
        r{2*N + i} = addreaction(model, strcat('a', istr, ' + c', istr, ' -> a', istr, ' + a', inxt));
        
    end

    % Add kinetic laws
    for i = 1:length(r)
        k{i} = addkineticlaw(r{i}, 'MassAction');
        k{i}.ParameterVariableNames = {strcat('l', int2str(i))};
        p{i} = addparameter(k{i}, strcat('l', int2str(i)), 'Value', rate);
    end

    % Initial concentrations
    for i = 1:N
        r{i}.Reactants(1).InitialAmount = initB*base;
        r{N+i}.Reactants(1).InitialAmount = initC*base;
        r{2*N + i}.Reactants(1).InitialAmount = initA*base;
    end

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
    
    function [val] = get_plot_values(X, names, s, n)
        val = zeros(length(X));
        for iter = 1:length(names)
            if strncmp(names{iter}, s, n)
                disp(strcat(names{iter}, ' should be ', s));
                val = val + X(:, iter);
            end
        end
    end

    % Plots
    figure(figNo);
    hold on; box on;
    a = get_plot_values(X, names, 'a', 1);
    b = get_plot_values(X, names, 'b', 1);
    c = get_plot_values(X, names, 'c', 1);
    
    a = a/N; 
    b = b/N;
    c = c/N;

    plot(t, a, 'LineWidth', 2.0);
    plot(t, b, ':', 'LineWidth', 2.0);
    plot(t, c, ':', 'LineWidth', 2.0);
    ylabel('concentration'); 
    xlabel('time');
    set(gca, 'LineWidth', 2.0); legend('A', 'B', 'C')

end

function dna_cat_rps(initA, initB, initC, N, stopTime, figNo, leak, shadow)
    
    % Rate
    rate = 1/60;

    % Base concentration
    base = 1.0;
    
    % Set tolerance levels
    absTol = 1e-12;
    relTol = 1e-12;
    % Define the model
    model = sbiomodel('Catalytic DNA RPS');
    
    % Infinite concentrations and rates
    infRate = 1e5*(rate);
    infConc = 1e5*base;

 
    
    % leak rate
    leak_rate = 1e-5;

    % Cells
    r = cell(1, 45*N);
    k = cell(1, length(r));
    p = cell(1, length(r));

    for i = 1:N
        istr = int2str(i);
        inxt = int2str(max(1, rem(i+1, N+1)));
        
        % B + A -> 2 B
        r{2*i-1} = addreaction(model, strcat('b', istr, ' + Ag', istr, ' -> Ib', istr));
        r{2*i} = addreaction(model, strcat('Ib', istr, ' + Gb', istr, ' -> b', istr, ' + b', inxt));

        % C + B -> 2 C
        r{2*N + 2*i-1} = addreaction(model, strcat('c', istr, ' + Bg', istr, ' -> Ic', istr));
        r{2*N + 2*i} = addreaction(model, strcat('Ic', istr, ' + Gc', istr, ' -> c', istr, ' + c', inxt));

        % A + C -> 2 A
        r{4*N + 2*i-1} = addreaction(model, strcat('a', istr, ' + Cg', istr, ' -> Ia', istr));
        r{4*N + 2*i} = addreaction(model, strcat('Ia', istr, ' + Ga', istr, ' -> a', istr, ' + a', inxt));

        % Linker reactions
        r{6*N + i} = addreaction(model, strcat('a', istr, ' -> Ag', istr));
        r{7*N + i} = addreaction(model, strcat('b', istr, ' -> Bg', istr));
        r{8*N + i} = addreaction(model, strcat('c', istr, ' -> Cg', istr));
        r{9*N + i} = addreaction(model, strcat('Ag', istr, ' -> a', istr));
        r{10*N + i} = addreaction(model, strcat('Bg', istr, ' -> b', istr));
        r{11*N + i} = addreaction(model, strcat('Cg', istr, ' -> c', istr));

        % Shadow B + A -> 2 B
        r{12*N + 2*i-1} = addreaction(model, strcat('sb', istr, ' + sAg', istr, ' -> sIb', istr));
        r{12*N + 2*i} = addreaction(model, strcat('sIb', istr, ' + sGb', istr, ' -> sb', istr, ' + sb', inxt));

        % Shadow C + B -> 2 C
        r{14*N + 2*i-1} = addreaction(model, strcat('sc', istr, ' + sBg', istr, ' -> sIc', istr));
        r{14*N + 2*i} = addreaction(model, strcat('sIc', istr, ' + sGc', istr, ' -> sc', istr, ' + sc', inxt));

        % Shadow A + C -> 2 A
        r{16*N + 2*i-1} = addreaction(model, strcat('sa', istr, ' + sCg', istr, ' -> sIa', istr));
        r{16*N + 2*i} = addreaction(model, strcat('sIa', istr, ' + sGa', istr, ' -> sa', istr, ' + sa', inxt));

        % Shadow Linker reactions
        r{18*N + i} = addreaction(model, strcat('sa', istr, ' -> sAg', istr));
        r{19*N + i} = addreaction(model, strcat('sb', istr, ' -> sBg', istr));
        r{20*N + i} = addreaction(model, strcat('sc', istr, ' -> sCg', istr));
        r{21*N + i} = addreaction(model, strcat('sAg', istr, ' -> sa', istr));
        r{22*N + i} = addreaction(model, strcat('sBg', istr, ' -> sb', istr));
        r{23*N + i} = addreaction(model, strcat('sCg', istr, ' -> sc', istr));

        % Regular leaks
        r{24*N + 2*i-1} = addreaction(model, strcat('Ag', istr, ' -> Ib', istr));
        r{24*N + 2*i} = addreaction(model, strcat('Gb', istr, ' -> b', istr, ' + b', inxt));
        r{26*N + 2*i-1} = addreaction(model, strcat('Bg', istr, ' -> Ic', istr));
        r{26*N + 2*i} = addreaction(model, strcat('Gc', istr, ' -> c', istr, ' + c', inxt));
        r{28*N + 2*i-1} = addreaction(model, strcat('Cg', istr, ' -> Ia', istr)); 
        r{28*N + 2*i} = addreaction(model, strcat('Ga', istr, ' -> a', istr, ' + a', inxt));

        % Shadow leaks
        r{30*N + 2*i-1} = addreaction(model, strcat('sAg', istr, ' -> sIb', istr));
        r{30*N + 2*i} = addreaction(model, strcat('sGb', istr, ' -> sb', istr, ' + sb', inxt));
        r{32*N + 2*i-1} = addreaction(model, strcat('sBg', istr, ' -> sIc', istr));
        r{32*N + 2*i} = addreaction(model, strcat('sGc', istr, ' -> sc', istr, ' + sc', inxt));
        r{34*N + 2*i-1} = addreaction(model, strcat('sCg', istr, ' -> sIa', istr)); 
        r{34*N + 2*i} = addreaction(model, strcat('sGa', istr, ' -> sa', istr, ' + sa', inxt));

        % Annihilation reactions
        r{36*N + i} = addreaction(model, strcat('a', istr, ' + sa', istr, ' -> w'));
        r{37*N + i} = addreaction(model, strcat('b', istr, ' + sb', istr, ' -> w'));
        r{38*N + i} = addreaction(model, strcat('c', istr, ' + sc', istr, ' -> w'));
        r{39*N + i} = addreaction(model, strcat('Ia', istr, ' + sIa', istr, ' -> w'));
        r{40*N + i} = addreaction(model, strcat('Ib', istr, ' + sIb', istr, ' -> w'));
        r{41*N + i} = addreaction(model, strcat('Ic', istr, ' + sIc', istr, ' -> w'));
        r{42*N + i} = addreaction(model, strcat('Ag', istr, ' + sAg', istr, ' -> w'));
        r{43*N + i} = addreaction(model, strcat('Bg', istr, ' + sBg', istr, ' -> w'));
        r{44*N + i} = addreaction(model, strcat('Cg', istr, ' + sCg', istr, ' -> w'));
        
        
    end



    % Add kinetic laws
    for i = 1:length(r)
        k{i} = addkineticlaw(r{i}, 'MassAction');
        k{i}.ParameterVariableNames = {strcat('l', int2str(i))}; 
    end

    % Add reaction rates
    for i = 1:N
        % B + A -> 2C
        p{2*i-1} = addparameter(k{2*i-1}, strcat('l', int2str(2*i-1)), 'Value', rate);
        p{2*i} = addparameter(k{2*i}, strcat('l', int2str(2*i)), 'Value', infRate);
        % C + B -> 2B
        p{2*N + 2*i-1} = addparameter(k{2*N + 2*i-1}, strcat('l', int2str(2*N + 2*i-1)), 'Value', rate);
        p{2*N + 2*i} = addparameter(k{2*N + 2*i}, strcat('l', int2str(2*N + 2*i)), 'Value', infRate);
        % A + C -> 2A 
        p{4*N + 2*i-1} = addparameter(k{4*N + 2*i-1}, strcat('l', int2str(4*N + 2*i-1)), 'Value', rate);
        p{4*N + 2*i} = addparameter(k{4*N + 2*i}, strcat('l', int2str(4*N + 2*i)), 'Value', infRate);
        % Linker reactions
        p{6*N + i} = addparameter(k{6*N + i}, strcat('l', int2str(6*N + i)), 'Value', infRate);
        p{7*N + i} = addparameter(k{7*N + i}, strcat('l', int2str(7*N + i)), 'Value', infRate);
        p{8*N + i} = addparameter(k{8*N + i}, strcat('l', int2str(8*N + i)), 'Value', infRate);
        p{9*N + i} = addparameter(k{9*N + i}, strcat('l', int2str(9*N + i)), 'Value', infRate);
        p{10*N + i} = addparameter(k{10*N + i}, strcat('l', int2str(10*N + i)), 'Value', infRate);
        p{11*N + i} = addparameter(k{11*N + i}, strcat('l', int2str(11*N + i)), 'Value', infRate);

        % Shadow B + A -> 2C
        p{12*N + 2*i-1} = addparameter(k{12*N + 2*i-1}, strcat('l', int2str(12*N + 2*i-1)), 'Value', shadow*rate);
        p{12*N + 2*i} = addparameter(k{12*N + 2*i}, strcat('l', int2str(12*N + 2*i)), 'Value', shadow*infRate);
        % Shadow C + B -> 2B
        p{14*N + 2*i-1} = addparameter(k{14*N + 2*i-1}, strcat('l', int2str(14*N + 2*i-1)), 'Value', shadow*rate);
        p{14*N + 2*i} = addparameter(k{14*N + 2*i}, strcat('l', int2str(14*N + 2*i)), 'Value', shadow*infRate);
        % Shadow A + C -> 2A 
        p{16*N + 2*i-1} = addparameter(k{16*N + 2*i-1}, strcat('l', int2str(16*N + 2*i-1)), 'Value', shadow*rate);
        p{16*N + 2*i} = addparameter(k{16*N + 2*i}, strcat('l', int2str(16*N + 2*i)), 'Value', shadow*infRate);
    end
    

    % Shadow linker reactions
    for i = 18*N + 1 : 24*N
        p{i} = addparameter(k{i}, strcat('l', int2str(i)), 'Value', shadow*infRate);
    end
    % Regular leak rates
    for i = 24*N + 1 : 30*N
        p{i} = addparameter(k{i}, strcat('l', int2str(i)), 'Value', leak*leak_rate);
    end

    % Shadow leak rates
    for i = 30*N + 1 : 36*N
        p{i} = addparameter(k{i}, strcat('l', int2str(i)), 'Value', shadow*leak*leak_rate);
    end

    % Annihilation reactions
    for i = 36*N + 1:45*N
        p{i} = addparameter(k{i}, strcat('l', int2str(i)), 'Value', shadow*leak*infRate);
    end



    % Initial concentrations
    for i = 1:N
        % Substrates
        r{2*i-1}.Reactants(1).InitialAmount = initB*base;
        r{2*N + 2*i - 1}.Reactants(1).InitialAmount = initC*base;
        r{4*N + 2*i - 1}.Reactants(1).InitialAmount = initA*base;
        
        % Gate complexes
        r{2*i}.Reactants(2).InitialAmount = infConc;
        r{2*N + 2*i}.Reactants(2).InitialAmount = infConc;
        r{4*N + 2*i}.Reactants(2).InitialAmount = infConc;

        % Shadow gate complexes 
        r{12*N + 2*i}.Reactants(2).InitialAmount = infConc;
        r{14*N + 2*i}.Reactants(2).InitialAmount = infConc;
        r{16*N + 2*i}.Reactants(2).InitialAmount = infConc;
    end

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
    
    function [val] = get_plot_values(X, names, s, n)
        val = zeros(length(X));
        for iter = 1:length(names)
            if strncmp(names{iter}, s, n)
                disp(strcat(names{iter}, ' should be ', s));
                val = val + X(:, iter);
            end
        end
    end

    % Plots
    figure(figNo);
    hold on; box on;
    a = get_plot_values(X, names, 'a', 1);
    b = get_plot_values(X, names, 'b', 1);
    c = get_plot_values(X, names, 'c', 1);
    
    a = a/N; 
    b = b/N;
    c = c/N;

    plot(t, a, 'LineWidth', 2.0);
    plot(t, b, ':', 'LineWidth', 2.0);
    plot(t, c, ':', 'LineWidth', 2.0);
    ylabel('concentration'); 
    xlabel('time');
    set(gca, 'LineWidth', 2.0); legend('A', 'B', 'C')

end

