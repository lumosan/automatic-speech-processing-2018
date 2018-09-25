function tfplot(s, fs, name, plottitle)
% TFPLOT Time and frequency plot
%    TFPLOT(S, FS, NAME, TITLE) displays a figure window with two
%    subplots.  Above, the signal S is plotted in time domain; below,
%    the signal is plotted in frequency domain. NAME is the "name" of the
%    signal, e.g., if NAME is 's', then the labels on the y-axes will be
%    's(t)' and '|s_F(f)|', respectively.  TITLE is the title that will
%    appear above the two plots.