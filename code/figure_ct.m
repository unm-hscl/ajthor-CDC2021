%% Run performance tests.
results = runperf('comp_time');

%% Plot results.

mean_times = zeros(1, length(results));

max_times = zeros(1, length(results));
min_times = zeros(1, length(results));

for k = 1:length(results)
    mean_times(k) = mean(results(k).Samples.MeasuredTime);
    
    max_times(k) = max(results(k).Samples.MeasuredTime);
    min_times(k) = min(results(k).Samples.MeasuredTime);
end

figure('Units', 'points', ...
       'Position', [0, 0, 240, 120]);
ax = axes;
ax.NextPlot = 'add';
ax.Units = 'points';

% ax.Position = [30, 30, 200, 200];

ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = 'Number Of Test Points';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = 'Computation Time [$$s$$]';
set(ax, 'FontSize', 8);

errorbar(1:50, mean_times, min_times, max_times);