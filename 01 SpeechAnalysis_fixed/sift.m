function [pitch] = sift(fname, frmsize, frmshift, lporder, sf, fignu, plottitle)
% This routine computes the pitch contour of a given speech signal. It takes
% as input the speech file name (fname), frame size (frmsize), frame shift(frmshift), linear prediction
% order (lporder), sampling frequency (sf), figure number (fignu) and the title of the plot (plottitle). The
% usage is as given below
%
% pitch = sift('fname', frmsize, frmshift, lporder, sf, fignu, 'plottitle')
%
% this routine returns the pitch contour.

sfby2 = sf/2;
wn = 800/sfby2;
%[b, a] = butter(6, wn);
 b(1) = 0.0357081667;
  b(2) = -0.0069956244;
  b(3) = b(2);
  b(4) = b(1);
  
  a(1) = 1.0;
  a(2) = -2.34036589;
  a(3) = 2.01190019;
  a(4) = -0.61419218;
%sp = filter(b, a, stdata);
tmpdata = load(fname);
data = tmpdata(17:length(tmpdata));
nsample = length(data);
maxframes=floor((nsample - frmsize + frmshift)/frmshift);
%pitch = zeros(maxframes, 1);
divisor = 1.7;
for n = 1:maxframes
loc = (n - 1)*frmshift;
sp = data(loc+1:loc+frmsize);
sp = filter(b, a, sp);
residual = lpResidual_SIFT(sp, frmsize, lporder, fignu, 'lp Residual');
residual = real(xcorr(residual, 256));

maxres = residual(256+1)/divisor;
maxloc = 256+1;
for k = 256+1+(sf/400)+1:256+1+(sf/80)
if(residual(k) > maxres)
     maxres = residual(k);
     maxloc = k;
end
end

     maxloc = maxloc - 256;
if(maxloc > 1)
     pitch(n) = 16000/maxloc;
else
pitch(n) = 0.0;
end

end
     maxpitch = max(pitch);
     figure(fignu+1)
     subplot(2, 1, 1);
     plot(data);
     title('Speech Signal');
     xlabel('time');
     ylabel('amplitude');
     axis([0 nsample -30000 30000])

     subplot(2, 1, 2)
     plot(pitch, 'r.');
     title(plottitle);
     xlabel('time');
     ylabel('Pitch Frequency in Hz');
     axis([0 maxframes 0 maxpitch])
