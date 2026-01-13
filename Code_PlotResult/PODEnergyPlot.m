%% Function: plot energy spertrum of POD modes (in MATLAB)
function PODEnergyPlot(energy)
    figure('Position', [100, 100, 1200, 500]);
    
    % figure a: plot POD mode energy spectrum
    subplot(1,2,1)
    loglog(1:length(energy), energy, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 3, 'MarkerFaceColor', 'blue');
    grid on;
    xlabel('\bf{Mode Order,} $\mathbf{m}$', 'Interpreter', 'latex');
    ylabel('\bf{Normalized Energy,} $\mathbf{ \lambda_i / \sum_{j=1}^m \lambda_j }$', 'Interpreter', 'latex');
    title('(a) POD Energy Spectrum', 'FontSize', 14, 'FontWeight', 'bold');
    set(gca, 'FontSize', 12, 'XMinorGrid', 'on', 'YMinorGrid', 'on');
    hold on
    plot([10 10], ylim, 'k--', 'LineWidth', 1, 'Color', [0, 0, 0]);
    plot([20 20], ylim, 'k--', 'LineWidth', 1, 'Color', [0, 0, 0]);
    plot([100 100], ylim, 'k--', 'LineWidth', 1, 'Color', [0, 0, 0]);
    loglog([10 20 100], energy([10 20 100]), ...
         'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'yellow', 'LineWidth', 1.5);
    hold off
    ylim([1e-15, 1e0]);
    set(gca, 'FontSize', 16, 'XMinorGrid', 'on');
    
    % figure b: plot cumulative energy spectrum
    subplot(1,2,2)
    cumulative_energy = cumsum(energy);
    semilogx(1:length(cumulative_energy), cumulative_energy, 'r-', 'LineWidth', 2.5);
    hold on;
    plot([10 10], [0 1], 'k--', 'LineWidth', 1, 'Color', [0, 0, 0]);
    plot([20 20], [0 1], 'k--', 'LineWidth', 1, 'Color', [0, 0, 0]);
    plot([100 100], [0 1], 'k--', 'LineWidth', 1, 'Color', [0, 0, 0]);
    plot([10 20 100], cumulative_energy([10 20 100]), ...
         'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'yellow', 'LineWidth', 1.5);
    grid on;
    xlabel('\bf{Mode Order,} $\mathbf{m}$', 'Interpreter', 'latex');
    ylabel('\bf{Cumulative Energy,} $\mathbf{ \sum_{k=1}^i \lambda_k / \sum_{j=1}^m \lambda_j }$', 'Interpreter', 'latex');
    title('(b) Cumulative Energy', 'FontSize', 14, 'FontWeight', 'bold'); 
    % text
    text(11, 0.2, sprintf('10 modes', cumulative_energy(10)*100), 'FontSize', 18, 'FontName', 'Times');
    text(22, 0.32, sprintf('20 modes', cumulative_energy(10)*100), 'FontSize', 18, 'FontName', 'Times');
    text(110, 0.63, sprintf('100 modes', cumulative_energy(10)*100), 'FontSize', 18, 'FontName', 'Times');
    set(gca, 'FontSize', 16, 'XMinorGrid', 'on');
    ylim([0, 1]);
    
    % figure c: POD energy spectrum for mode 1-10
    figure('Position', [100, 100, 300, 200]);
    plot(1:10, energy(1:10), 'b-o', 'LineWidth', 2, 'MarkerSize', 5, 'MarkerFaceColor', 'blue');
    xlim([1 10])
    ylim([0,0.2]);
    grid on
    set(gca, 'FontSize', 16);
    
end