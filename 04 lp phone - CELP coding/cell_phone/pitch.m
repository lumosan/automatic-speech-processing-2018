function T0=pitch(speech_frame);
% T0=pitch(speech_frame) 
%
% This function estimates the fundamental period (in samples) of a 30 ms
% speech frame sampled at 8 kHz. T0=0 if the frame is detexted as unvoiced.
% T0 is computed from the max. of the autocorrelation of the LPC residual.
% Voiced/unvoiced decision is based on the the ratio of this maximum by the
% variance of the residual. Pitch is searched for in the range
% [60Hz,300Hz], i.e. for autocorrelation indices in [26,133].
%
% 12:10 21/02/2007 T. Dutoit

ai=lpc(speech_frame,10);
lpc_residual=filter(ai,1,speech_frame);
C=xcorr(lpc_residual,133); % autocorrelation
Cxx=C(133+1:2*133+1)/C(134); % relative to its variance
Cxx(1:26)=0;
[Amax,Imax]=max(Cxx); % position of the maximum
if (Amax>0.20) % *very* rough V/UV condition.
   T0=Imax-1;
else 
   T0=0;
end;
