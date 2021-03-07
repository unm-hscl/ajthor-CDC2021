%% Run performance tests.
results = runperf('comp_time_funcofR/testOptimalAlgorithmVaryingTestPoints');

%% Plot results.

mean_times = zeros(1, length(results));

max_times = zeros(1, length(results));
min_times = zeros(1, length(results));

for k = 1:length(results)
    mean_times(k) = mean(results(k).Samples.MeasuredTime);
    
    max_times(k) = max(results(k).Samples.MeasuredTime) - mean_times(k);
    min_times(k) = mean_times(k) - min(results(k).Samples.MeasuredTime);
end

figure('Units', 'points', ...
       'Position', [0, 0, 120, 120]);
ax = axes;
ax.NextPlot = 'add';
ax.Units = 'points';

% ax.Position = [30, 30, 200, 200];

ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = 'Number Of Eval. Points';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = 'Computation Time [$$s$$]';
set(ax, 'FontSize', 8);

errorbar(1:5:51, mean_times, min_times, max_times);

%% Run performance tests.
results = runperf('comp_time_funcofR/testAlgorithmVaryingTestPoints');

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

%% Save the figure as 'figure3a'.
saveas(gcf, '../results/figure3a.png')
savefig('../results/figure3a.fig')

%% Run performance tests.
results = runperf('comp_time_funcofM');

%% Plot results.

mean_times = zeros(1, length(results));

max_times = zeros(1, length(results));
min_times = zeros(1, length(results));

for k = 1:length(results)
    mean_times(k) = mean(results(k).Samples.MeasuredTime);
    
    max_times(k) = max(results(k).Samples.MeasuredTime) - mean_times(k);
    min_times(k) = mean_times(k) - min(results(k).Samples.MeasuredTime);
end

figure('Units', 'points', ...
       'Position', [0, 0, 120, 120]);
ax = axes;
ax.NextPlot = 'add';
ax.Units = 'points';

% ax.Position = [30, 30, 200, 200];

ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = 'Sample Size $$M$$';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = 'Computation Time [$$s$$]';
set(ax, 'FontSize', 8);

errorbar(100:500:8100, mean_times, min_times, max_times);

%% Save the figure as 'figure3b'.
saveas(gcf, '../results/figure3b.png')
savefig('../results/figure3b.fig')