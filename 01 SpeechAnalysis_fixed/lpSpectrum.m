function [lpspec] = lpSpectrum(data, lporder, win, order, sf, fignu, plottitle)

% this routine computes the linear prediction coefficients and plots
% the lp specturm in the figure with number fignu and plottitle as the
% title, along with the fourier spectrum. The usage is as given below
%
% lpspec = lpSpectrum(data, lporder, win, order, sf, fignu, 'plottitle');
%
% where data is the short-time signal, lporder is the order of the linear
% predicition, win is 0 for rectangular window and 1 for Hamming window.
% order is the FFT order used during the computation of lp spectrum,
% sf is the sampling frequency the computed lp spectrum is plotted in 
% figure fignu with plottitle as the title.

if win == 1
 data = data .* hamming(length(data));
end

% preemphasize to remove the spectral tilt due to glottal pulse spectrum.
difdata(1) = data(1);
difdata(2:length(data)) = diff(data);
if size(difdata, 2)>1
   difdata = difdata';
end

for n=1:order freq(n) = ((n-1) * sf)/(order * 1000.0); end

% Computation of lp Spectrum
a = real(lpc(difdata, lporder))
     lpspec = -20 * log10(abs(fft(a, order)));
     %lpspec = -20 * log10(abs(freqz(a, (order/2))));

% Computation of fourier Spectrum
fourSpec = 20 * log10(abs(fft(data, order)));

figure(fignu)
     subplot(2, 1, 1)
     plot(freq(1:(order/2)), fourSpec(1:(order/2)));
     xlabel('Frequency (in kHz)');
     ylabel('Log Amplitude');
     title('Fourier Spectrum');
     subplot(2, 1, 2)
     plot(freq(1:(order/2)), lpspec(1:(order/2)));
     xlabel('Frequency (in kHz)');
     ylabel('Log Amplitude');
     title(plottitle);
     zoom on;





