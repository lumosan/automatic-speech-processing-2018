function tfplot(s, fs, name, plottitle)
% TFPLOT Time and frequency plot
%    TFPLOT(S, FS, NAME, TITLE) displays a figure window with two
%    subplots. Above, the signal S is plotted in time domain; below,
%    the signal is plotted in frequency domain. NAME is the "name" of the
%    signal, e.g., if NAME is 's', then the labels on the y-axes will be
%    's(t)' and '|s_F(f)|', respectively.  TITLE is the title that will
%    appear above the two plots.

	% Obtain time vector
	t = 0:1/fs:1/fs*(length(s)-1);
    
    % Obtain frequency vector
    f = linspace(-fs/2, fs/2, length(s));
    
    % First plot
    figure
	subplot(2,1,1);
	plot(t, s)
    title(plottitle)
	xlabel('t [s]')
	ylabel(strcat(name, '(t)'))
    
    % Second plot
    subplot(2,1,2);
    spectrum = fftshift(abs(fft(s)));
    % TODO: modify the spectrum somehow?
	plot(f, spectrum/fs)
	xlabel('f [Hz]')
	ylabel(strcat('|', name, '_F(f)|'))
end