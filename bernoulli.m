n = 1:1:20;
hold on; box on;
xlim([0 20]); ylim([0 1]); set(gca, 'LineWidth', 2.0); 

p = (0.1).^n;
plot(n, p, 'LineWidth', 2)
p = (0.3).^n;
plot(n, p, 'LineWidth', 2)
p = (0.5).^n;
plot(n, p, 'LineWidth', 2)
p = (0.8).^n;
plot(n, p, 'LineWidth', 2)

ylabel('Failure probability of the system')
xlabel('Number of components')