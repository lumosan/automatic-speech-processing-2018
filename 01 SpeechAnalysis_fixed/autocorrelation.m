function [corr] = autocorrelation(data, order, fignu, plottitle);

% this routine computes the autocorrelation of the given signal and
% plots it, see below the usage of this routine
%
% autocorrelation(data, order, fignu, 'plottitle');
%
% data is the array containing the region of signal for which
% autocorrelation has to be computed. npoint is the length of the 
% data array. fignu is the figure number in which the autocorrelation
% is plotted, with a title plottitle.

corr = real(xcorr(data, order));
figure(fignu)
     subplot(2, 1, 1)
     plot(corr);
     xlabel('sequence number');
     ylabel('amplitude');
     subplot(2, 1, 2)
     plot(corr(order+1:(order+order+ 1)), 'r');
     grid on;
     title(plottitle);
     xlabel('sequence number');
     ylabel('amplitude');
     zoom on;

