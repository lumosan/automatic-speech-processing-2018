function [residual] = lpResidual(stdata, nsample, lporder, fignu, plottitle)

% This routine takes as input the speech file name, the short-time speech 
% signal, the duration of the short-time
% signal expressed in number of samples, the linear prediction order,
% figure and the title of the plot. See the usage below
%
%  residual = lpResidual(stdata, nsample, lporder, fignu, 'plottitle')
%
% This routine returns the linear prediction residual

%     data = load(fname);
%     data = data(17:length(data));
%     stdata = data(begin:end);
     

%sp = stdata .* hamming(nsample);
%a = lpc(sp, lporder);

% Preemphasize to remove the spectral tilt due to glottal pulse spectrum.
difsp(1) = stdata(1);
difsp(2:nsample) = diff(stdata);
if size(difsp,2)>1
   difsp = difsp';   % Convert to a column vector.
end

%find the autocorrelation coefficients
sp = difsp .* hamming(nsample);
arcoef = real(lpc(sp, lporder));



prevsig = zeros(lporder, 1);
augsignal = prevsig;
augsignal(lporder+1:lporder+nsample) = stdata;
for n=lporder+1:lporder+nsample,
   predict = 0;
   for j=1:lporder,
      predict = predict + arcoef(j+1)*augsignal(n-j);
   end
   residual(n-lporder) = augsignal(n) + predict;
end
length(residual)
figure(fignu)
subplot(2, 1, 1)
     plot(stdata, 'b');
     xlabel('time');
     ylabel('amplitude');
     title('short time speech signal');
     grid on;
subplot(2, 1, 2)
     plot(residual, 'r');
     xlabel('time');
     ylabel('amplitude');
     title(plottitle);
     grid on;
