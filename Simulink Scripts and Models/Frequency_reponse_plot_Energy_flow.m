% Assuming 'logsout' is your Simulink logging object
data = logsout{1};
data = data.Values;
data = data.Data;

% Plotting the data
% Creating a time vector from 0 to 1 second with a step of 2e-6 seconds
time_vector = 0:2e-6:1;

% Ensure that the time_vector and data have the same length
if length(time_vector) > length(data)
    time_vector = time_vector(1:length(data));
elseif length(time_vector) < length(data)
    data = data(1:length(time_vector));
end

figure; % Create a new figure
plot(time_vector, data); % Plot the data against the time vector

% Log-log plot
% Creating a frequency vector from 1 to 10000 Hz
frequency_vector = 1:9999/5e5:10000;

% Ensure that the frequency_vector and data have the same length
if length(frequency_vector) > length(data)
    frequency_vector = frequency_vector(1:length(data));
elseif length(frequency_vector) < length(data)
    data = data(1:length(frequency_vector));
end

figure; % Create a new figure for the log-log plot
loglog(frequency_vector, data); % Create a log-log plot
xlim([100, 5000]); % Set the x-axis limits
ylabel('Response');
xlabel('Frequency (Hz)');

% Add warning if negative data is present
if any(data < 0)
    warning('Negative data ignored');
end
