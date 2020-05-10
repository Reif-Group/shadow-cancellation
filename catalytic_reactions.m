clf;
crnidealrpsoscillator(10, 1);


function crnidealrpsoscillator(stopTime, figNo)
    rate = 3;
    model = sbiomodel('Testing Catalytic Reactions');
    
    % reaction set
    r1 = addreaction(model,'x0 -> x0 + x1'); 
    r2 = addreaction(model,'x1 -> x1 + x2'); 
    r3 = addreaction(model,'x2 -> x2 + x3'); 
    r4 = addreaction(model,'x3 -> x3 + x4'); 
    r5 = addreaction(model,'x4 -> x4 + x5');
    r6 = addreaction(model,'x5 -> x5 + x6'); 
    r7 = addreaction(model,'x6 -> x6 + x7'); 
    r8 = addreaction(model,'x7 -> x7 + x8'); 
    r9 = addreaction(model,'x8 -> x8 + x9'); 
    r10 = addreaction(model,'x9 -> x9 + x10'); 
   
    
    % kinetic laws
    k1 = addkineticlaw(r1, 'MassAction'); 
    k1.ParameterVariableNames = {'c1'};
    k2 = addkineticlaw(r2, 'MassAction'); 
    k2.ParameterVariableNames = {'c2'};
    k3 = addkineticlaw(r3, 'MassAction'); 
    k3.ParameterVariableNames = {'c3'};
    k4 = addkineticlaw(r4, 'MassAction'); 
    k4.ParameterVariableNames = {'c4'};
    k5 = addkineticlaw(r5, 'MassAction'); 
    k5.ParameterVariableNames = {'c5'};
    k6 = addkineticlaw(r6, 'MassAction'); 
    k6.ParameterVariableNames = {'c6'};
    
    k7 = addkineticlaw(r7, 'MassAction'); 
    k7.ParameterVariableNames = {'c7'};
    k8 = addkineticlaw(r8, 'MassAction'); 
    k8.ParameterVariableNames = {'c8'};
    k9 = addkineticlaw(r9, 'MassAction'); 
    k9.ParameterVariableNames = {'c9'};
    k10 = addkineticlaw(r10, 'MassAction'); 
    k10.ParameterVariableNames = {'c10'};
    
    p1 = addparameter(k1, 'c1', 'Value',  rate);
    p2 = addparameter(k2, 'c2', 'Value', rate);
    p3 = addparameter(k3, 'c3', 'Value', rate);
    p4 = addparameter(k4, 'c4', 'Value', rate);
    p5 = addparameter(k5, 'c5', 'Value', rate);
    p6 = addparameter(k6, 'c6', 'Value', rate);
    p7 = addparameter(k7, 'c7', 'Value', rate);
    p8 = addparameter(k8, 'c8', 'Value', rate);
    p9 = addparameter(k9, 'c9', 'Value', rate);
    p10 = addparameter(k10, 'c10', 'Value', rate);

    % initial amounts
    base = 1;
    r1.Reactants(1).InitialAmount = base;
    r2.Reactants(1).InitialAmount = base;
    r3.Reactants(1).InitialAmount = base;
    r4.Reactants(1).InitialAmount = base;
    r5.Reactants(1).InitialAmount = base;
    r6.Reactants(1).InitialAmount = base;
    r7.Reactants(1).InitialAmount = base;
    r8.Reactants(1).InitialAmount = base;
    r9.Reactants(1).InitialAmount = base;
    r10.Reactants(1).InitialAmount = base;
    r10.Products(2).InitialAmount = base;

    

    model;
    model.Reactions
    model.Species

    cs = getconfigset(model,'active'); 
    odeset('RelTol',2e+20);
    cs.SolverType = 'ode15s';
    cs.StopTime = stopTime;
    [t,X, names] = sbiosimulate(model);
    
    % plot
    figure(figNo);
    box on; hold on;
    
    plot(t, X, 'LineWidth', 2.0);
    ylabel('Concentration'); xlabel('Time');
    set(gca, 'LineWidth', 2.0); 
    legend(names)
    title('oscillation')
   
    plot([0:0.01:10], 1 + [0:0.01:10], ':', 'LineWidth', 2.0)
    x = 1 + [0:0.01:10]; y = x.^2;
    plot([0:0.01:10], y, ':', 'LineWidth', 2.0)
    x = 1 + [0:0.01:10]; y = x.^3;
    plot([0:0.01:10], y, ':', 'LineWidth', 2.0)
    x = 1+[0:0.01:10]; y = exp(x);
    plot([0:0.01:10], y, ':', 'LineWidth', 2.0)
    
    ylim([0 100]); xlim([0, 10])
    
    
    
end