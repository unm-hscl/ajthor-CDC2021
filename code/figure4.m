%% Figure (4)
% Figure showing the computation time of the kernel based approach as a
% function of the number of evaluation points.

%% Run performance tests.
results = runperf('comp_time_funcofR/testOptimalAlgorithmVaryingTestPoints');
save('../results/comp_funcofR_cvx.mat', 'results');

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
ax.XLabel.String = 'Number Of Eval. Points';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = 'Computation Time [$$s$$]';
set(ax, 'FontSize', 8);

errorbar(1:5:51, mean_times, min_times, max_times);

%% Run performance tests.
results = runperf('comp_time_funcofR/testAlgorithmVaryingTestPoints');
save('../results/comp_funcofR_alg.mat', 'results');

%% Plot results.

mean_times = zeros(1, length(results));

max_times = zeros(1, length(results));
min_times = zeros(1, length(results));

for k = 1:length(results)
    mean_times(k) = mean(results(k).Samples.MeasuredTime);
    
    max_times(k) = max(results(k).Samples.MeasuredTime) - mean_times(k);
    min_times(k) = mean_times(k) - min(results(k).Samples.MeasuredTime);
end

errorbar(1:5:51, mean_times, min_times, max_times);

%% Save the figure as 'figure4'.
saveas(gcf, '../results/figure4.png')
savefig('../results/figure4.fig')
