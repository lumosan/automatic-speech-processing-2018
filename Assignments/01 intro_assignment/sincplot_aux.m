% SINCPLOT Time and frequency plot of a sinc signal
%    SINCPLOT() creates a new figure window with two subplots, which show
%    time and frequency plots for a sinc function.

clear all;

% Define constants
F_SAMPLE = 5000;
DURATION = 0.05;
B = 1000;

% Set time and frequency scales
t = linspace(0, DURATION, DURATION*F_SAMPLE + 1);

% Define the signal
m = sinc(B*(t - max(t)/2));

% Plot
tfplot(m,F_SAMPLE,'m','A sinc signal')
