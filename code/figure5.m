%% Figure (5)
% Figure showing the computation time of the kernel based approach as a
% function of the sample size.

%% Run performance tests.
results = runperf('comp_time_funcofM');
save('../results/comp_funcofM_alg.mat', 'results');

%% Plot results.

mean_times = zeros(1, length(results));

max_times = zeros(1, length(results));
min_times = zeros(1, length(results));

for k = 1:length(results)
    mean_times(k) = mean(results(k).Samples.MeasuredTime);
    
    max_times(k) = max(results(k).Samples.MeasuredTime) - mean_times(k);
    min_times(k) = mean_times(k) - min(results(k).Samples.MeasuredTime);
end

figure;
ax = axes;
ax.NextPlot = 'add';
ax.Units = 'points';

ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = 'Sample Size $$M$$';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = 'Computation Time [$$s$$]';
set(ax, 'FontSize', 8);

errorbar(100:500:8100, mean_times, min_times, max_times);

%% Save the figure as 'figure5'.
saveas(gcf, '../results/figure5.png')
savefig('../results/figure5.fig')
