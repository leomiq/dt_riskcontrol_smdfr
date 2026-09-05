figure
hold on
for i = 50:10:500
    scatter(i,CalculateTs(i,543), 'blue')
end

% Add labels
xlabel('Opening Time (s.)', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('Required temperature (fuel interface) [K]', 'Interpreter', 'Latex', 'FontSize', 15);

%%
figure
hold on
for i = 55:1:80
    scatter(i,CurrentToMassFlow(i), 'blue')
end

% Add labels
xlabel('current [A]', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('$\dot{m}_{\kappa, \rm ACS}$, [kg/s]', 'Interpreter', 'Latex', 'FontSize', 15);