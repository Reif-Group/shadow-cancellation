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
    shadow = 1;
    annihilation = 1;
    
    % Rate
    rateConstant = 1.0;
    rate = rateConstant*alpha;

    % Set infConc and infRate
    infConc = 1e5*base;
    infRate = 1e5*rate;

    % Set leak rate
    leak_rate = 1e-6;

    % Define the model
    model = sbiomodel('Ideal Consensus');
    
    % Cells
    r = cell(1, 24*N);
    k = cell(1, length(r));
    p = cell(1, length(r));
    
    % Regular circuit
    for i = 1:N
        istr = int2str(i);
        inxt = int2str(max(1, rem(i+1, N+1)));
        r{2*i-1} = addreaction(model, strcat('a', istr, ' + Bg', istr, ' -> Iy', istr));
        r{2*i} = addreaction(model, strcat('Iy', istr, ' + Gy', istr, ' -> y', istr, ' + y', inxt));
        r{2*N + 2*i-1} = addreaction(model, strcat('b', istr, ' + Yg', istr, ' -> Ib',istr));
        r{2*N + 2*i} = addreaction(model, strcat('Ib', istr, ' + Gb', istr, ' -> b', istr, ' + b', inxt));
        r{4*N + 2*i-1} = addreaction(model, strcat('y', istr, ' + Ag', istr, ' -> Ia', istr));
        r{4*N + 2*i} = addreaction(model, strcat('Ia', istr, ' + Ga', istr, ' -> a', istr, ' + a', inxt));
        % Linker reactions
        r{6*N + 2*i-1} = addreaction(model, strcat('a', istr, ' -> Ag', istr));
        r{6*N + 2*i} = addreaction(model, strcat('Ag', istr, ' -> a', istr));
        r{8*N + 2*i-1} = addreaction(model, strcat('b', istr, ' -> Bg', istr));
        r{8*N + 2*i} = addreaction(model, strcat('Bg', istr, ' -> b', istr));
        r{10*N + 2*i-1} = addreaction(model, strcat('y', istr, ' -> Yg', istr));
        r{10*N + 2*i} = addreaction(model, strcat('Yg', istr, ' -> y', istr));
        % Shadow circuit
        r{12*N + 2*i-1} = addreaction(model, strcat('sa', istr, ' + sBg', istr, ' -> sIy', istr));
        r{12*N + 2*i} = addreaction(model, strcat('sIy', istr, ' + sGy', istr, ' -> sy', istr, ' + sy', inxt));
        r{14*N + 2*i-1} = addreaction(model, strcat('sb', istr, ' + sYg', istr, ' -> sIb',istr));
        r{14*N + 2*i} = addreaction(model, strcat('sIb', istr, ' + sGb', istr, ' -> sb', istr, ' + sb', inxt));
        r{16*N + 2*i-1} = addreaction(model, strcat('sy', istr, ' + sAg', istr, ' -> sIa', istr));
        r{16*N + 2*i} = addreaction(model, strcat('sIa', istr, ' + sGa', istr, ' -> sa', istr, ' + sa', inxt));
        % Shadow Linker
        r{18*N + 2*i-1} = addreaction(model, strcat('sa', istr, ' -> sAg', istr));
        r{18*N + 2*i} = addreaction(model, strcat('sAg', istr, ' -> sa', istr));
        r{20*N + 2*i-1} = addreaction(model, strcat('sb', istr, ' -> sBg', istr));
        r{20*N + 2*i} = addreaction(model, strcat('sBg', istr, ' -> sb', istr));
        r{22*N + 2*i-1} = addreaction(model, strcat('sy', istr, ' -> sYg', istr));
        r{22*N + 2*i} = addreaction(model, strcat('sYg', istr, ' -> sy', istr));
    end

    % Leak and Shadow
    for i = 1:N
        istr = int2str(i);
        inxt = int2str(max(1, rem(i+1, N+1)));
        % Gy1 -> y1 + y2
        r{24*N + i} = addreaction(model, strcat('Gy', istr, ' -> y', istr, ' + y', inxt));
        r{25*N + i} = addreaction(model, strcat('sGy', istr, ' -> sy', istr, ' + sy', inxt));
        r{26*N + i} = addreaction(model, strcat('y', istr, ' + sy', istr, ' -> w'));
        % Ga1 -> a1 + a2
        r{27*N + i} = addreaction(model, strcat('Ga', istr, ' -> a', istr, ' + a', inxt));
        r{28*N + i} = addreaction(model, strcat('sGa', istr, ' -> sa', istr, ' + sa', inxt));
        r{29*N + i} = addreaction(model, strcat('a', istr, ' + sa', istr, ' -> w'));
        % Gb1 -> b1 + b2
        r{30*N + i} = addreaction(model, strcat('Gb', istr, ' -> b', istr, ' + b', inxt));
        r{31*N + i} = addreaction(model, strcat('sGb', istr, ' -> sb', istr, ' + sb', inxt));
        r{32*N + i} = addreaction(model, strcat('b', istr, ' + sb', istr, ' -> w'));

    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% Leak Profile %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Add rate kinetics
    for i = 1:length(r)
        k{i} = addkineticlaw(r{i}, 'MassAction');
        k{i}.ParameterVariableNames = {strcat('c', int2str(i))};
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% Rate constants %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    gamma = 2;
    bimol_scaled_rate = N*rate; % When xi = x/N, yi = y/N
    bimol = gamma*bimol_scaled_rate; % Accounting for the presence of y strand in both gate form and strand form
    
    for i = 1:N
        % Regular reactions
        p{2*i-1} = addparameter(k{2*i-1}, strcat('c', int2str(2*i-1)), 'Value', bimol);
        p{2*i} = addparameter(k{2*i}, strcat('c', int2str(2*i)), 'Value', infRate);
        p{2*N + 2*i-1} = addparameter(k{2*N + 2*i-1}, strcat('c', int2str(2*N + 2*i-1)), 'Value', bimol);
        p{2*N + 2*i} = addparameter(k{2*N + 2*i}, strcat('c', int2str(2*N + 2*i)), 'Value', infRate);
        p{4*N + 2*i-1} = addparameter(k{4*N + 2*i-1}, strcat('c', int2str(4*N + 2*i-1)), 'Value', bimol);
        p{4*N + 2*i} = addparameter(k{4*N + 2*i}, strcat('c', int2str(4*N + 2*i)), 'Value', infRate);
        % Linker reactions
        p{6*N + 2*i-1} = addparameter(k{6*N + 2*i-1}, strcat('c', int2str(6*N + 2*i-1)), 'Value', infRate);
        p{6*N + 2*i} = addparameter(k{6*N + 2*i}, strcat('c', int2str(6*N + 2*i)), 'Value', infRate);
        p{8*N + 2*i-1} = addparameter(k{8*N + 2*i-1}, strcat('c', int2str(8*N + 2*i-1)), 'Value', infRate);
        p{8*N + 2*i} = addparameter(k{8*N + 2*i}, strcat('c', int2str(8*N + 2*i)), 'Value', infRate);
        p{10*N + 2*i-1} = addparameter(k{10*N + 2*i-1}, strcat('c', int2str(10*N + 2*i-1)), 'Value', infRate);
        p{10*N + 2*i} = addparameter(k{10*N + 2*i}, strcat('c', int2str(10*N + 2*i)), 'Value', infRate);

        % Shadow Regular reactions
        p{12*N + 2*i-1} = addparameter(k{12*N + 2*i-1}, strcat('c', int2str(12*N + 2*i-1)), 'Value', shadow*bimol);
        p{12*N + 2*i} = addparameter(k{12*N + 2*i}, strcat('c', int2str(12*N + 2*i)), 'Value', shadow*infRate);
        p{14*N + 2*i-1} = addparameter(k{14*N + 2*i-1}, strcat('c', int2str(14*N + 2*i-1)), 'Value', shadow*bimol);
        p{14*N + 2*i} = addparameter(k{14*N + 2*i}, strcat('c', int2str(14*N + 2*i)), 'Value', shadow*infRate);
        p{16*N + 2*i-1} = addparameter(k{16*N + 2*i-1}, strcat('c', int2str(16*N + 2*i-1)), 'Value', shadow*bimol);
        p{16*N + 2*i} = addparameter(k{16*N + 2*i}, strcat('c', int2str(16*N + 2*i)), 'Value', shadow*infRate);
        % Shadow Linker reactions
        p{18*N + 2*i-1} = addparameter(k{18*N + 2*i-1}, strcat('c', int2str(18*N + 2*i-1)), 'Value', shadow*infRate);
        p{18*N + 2*i} = addparameter(k{18*N + 2*i}, strcat('c', int2str(18*N + 2*i)), 'Value', shadow*infRate);
        p{20*N + 2*i-1} = addparameter(k{20*N + 2*i-1}, strcat('c', int2str(20*N + 2*i-1)), 'Value', shadow*infRate);
        p{20*N + 2*i} = addparameter(k{20*N + 2*i}, strcat('c', int2str(20*N + 2*i)), 'Value', shadow*infRate);
        p{22*N + 2*i-1} = addparameter(k{22*N + 2*i-1}, strcat('c', int2str(22*N + 2*i-1)), 'Value', shadow*infRate);
        p{22*N + 2*i} = addparameter(k{22*N + 2*i}, strcat('c', int2str(22*N + 2*i)), 'Value', shadow*infRate);

        % Leaks, shadow and Annihilation
        % Gy1 -> y1 + y2
        p{24*N + i} = addparameter(k{24*N + i}, strcat('c', int2str(24*N + i)), 'Value', leak*leak_rate);
        p{25*N + i} = addparameter(k{25*N + i}, strcat('c', int2str(25*N + i)), 'Value', leak*shadow*leak_rate);
        p{26*N + i} = addparameter(k{26*N + i}, strcat('c', int2str(26*N + i)), 'Value', leak*shadow*infRate*1e3);
        % Ga1 -> a1 + a2
        p{27*N + i} = addparameter(k{27*N + i}, strcat('c', int2str(27*N + i)), 'Value', leak*leak_rate);
        p{28*N + i} = addparameter(k{28*N + i}, strcat('c', int2str(28*N + i)), 'Value', leak*shadow*leak_rate);
        p{29*N + i} = addparameter(k{29*N + i}, strcat('c', int2str(29*N + i)), 'Value', leak*shadow*infRate*1e3);
        % Gb1 -> b1 + b2
        p{30*N + i} = addparameter(k{30*N + i}, strcat('c', int2str(30*N + i)), 'Value', leak*leak_rate);
        p{31*N + i} = addparameter(k{31*N + i}, strcat('c', int2str(31*N + i)), 'Value', leak*shadow*leak_rate);
        p{32*N + i} = addparameter(k{32*N + i}, strcat('c', int2str(32*N + i)), 'Value', leak*shadow*infRate*1e3);
        
        
    end

    % Add initial amounts
    for i = 1:N
        % Substrates
        r{2*i-1}.Reactants(1).InitialAmount = initA/N;
        r{2*N + 2*i-1}.Reactants(1).InitialAmount = initB/N;
        r{4*N + 2*i-1}.Reactants(1).InitialAmount = initY/N;

        % Regular gates
        r{2*i}.Reactants(2).InitialAmount = infConc;
        r{2*N + 2*i}.Reactants(2).InitialAmount = infConc;
        r{4*N + 2*i}.Reactants(2).InitialAmount = infConc;
        
        % Shadow gates
        r{12*N + 2*i}.Reactants(2).InitialAmount = infConc;
        r{14*N + 2*i}.Reactants(2).InitialAmount = infConc;
        r{16*N + 2*i}.Reactants(2).InitialAmount = infConc;
       
       
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
    

    function [val] = get_plot_values(X, names, s, n)
        val = zeros(length(X));
        for iter = 1:length(names)
            if strncmp(names{iter}, s, n)
                disp(strcat(names{iter}, ' should be ', s));
                val = val + X(:, iter);
            end
        end
    end
    

    a = get_plot_values(X, names, 'a', 1);
    ag = get_plot_values(X, names, 'Ag', 2);
    b = get_plot_values(X, names, 'b', 1);
    bg = get_plot_values(X, names, 'Bg', 2);
    y = get_plot_values(X, names, 'y', 1);
    yg = get_plot_values(X, names, 'Yg', 2);
    
    atotal = (a + ag);
    btotal = (b + bg);
    ytotal = (y + yg);

    % Shadow variables
    sa = get_plot_values(X, names, 'sa', 2);
    sag = get_plot_values(X, names, 'sAg', 3);
    sb = get_plot_values(X, names, 'sb', 2);
    sbg = get_plot_values(X, names, 'sBg', 3);
    sy = get_plot_values(X, names, 'sy', 2);
    syg = get_plot_values(X, names, 'sYg', 3);
    
    satotal = (sa + sag);
    sbtotal = (sb + sbg);
    sytotal = (sy + syg);

    % Set figure
    figure(figNo);
    box on; hold on;

    plot(t/3600, atotal/1.e-9, 'LineWidth', 2.0);
    plot(t/3600, btotal/1.e-9, 'LineWidth', 2.0, 'Color', 'red');
    plot(t/3600, ytotal/1.e-9, ':', 'LineWidth', 2.0, 'Color', 'green');
%     plot(t/3600, satotal/1.e-9, 'LineWidth', 2.0, 'Color', 'blue');
%     plot(t/3600, sbtotal/1.e-9, 'LineWidth', 2.0, 'Color', 'red');
%     plot(t/3600, sytotal/1.e-9, ':', 'LineWidth', 2.0, 'Color', 'green');
%     plot(t/3600, ag/1.e-9, 'LineWidth', 2.0);
%     plot(t/3600, bg/1.e-9, 'LineWidth', 2.0, 'Color', 'red');
%     plot(t/3600, yg/1.e-9, ':', 'LineWidth', 2.0, 'Color', 'green');
%     plot(t/3600, atotal/1.e-9, 'LineWidth', 2.0, 'Color', 'blue');
%     plot(t/3600, btotal/1.e-9, 'LineWidth', 2.0, 'Color', 'red');
%     plot(t/3600, ytotal/1.e-9, ':', 'LineWidth', 2.0, 'Color', 'green');
%     plot(t/3600, sa/1e-9, ':', 'LineWidth', 2.0);
%     plot(t/3600, sag/1e-9, ':', 'LineWidth', 2.0);
    ylabel('concentration'); 
    xlabel('time');



end
