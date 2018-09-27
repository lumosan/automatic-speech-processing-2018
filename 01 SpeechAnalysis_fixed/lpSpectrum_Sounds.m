function [] = lpSpectrum_Sounds(fignu)

% this routine plots the linear prediction spectrum of different speech
% sounds. The usage is as given below
%
% lpSpectrum_Sounds(fignu)
%
% the variable to be passed to this routine is just the figure number
% fignu

spdata = load('m_s2_t1_a');
stspdata = spdata(9001:9480);
lpSpectrum(stspdata, 16, 1, 512, 16000, fignu,'LP spectrum of sound /a/');

fignu = fignu + 1;
spdata = load('m_s2_t1_e');
stspdata = spdata(8001:8480);
lpSpectrum(stspdata, 16, 1, 512, 16000, fignu,'LP spectrum of sound /e/');

fignu = fignu + 1;
spdata = load('m_s2_t1_i');
stspdata = spdata(14001:14480);
lpSpectrum(stspdata, 14, 1, 512, 16000, fignu,'LP spectrum of sound /i/');

fignu = fignu + 1;
spdata = load('m_s2_t1_o');
stspdata = spdata(12001:12480);
lpSpectrum(stspdata, 18, 1, 512, 16000, fignu,'LP spectrum of sound /o/');

fignu = fignu + 1;
spdata = load('m_s2_t1_u');
stspdata = spdata(17001:17480);
lpSpectrum(stspdata, 18, 1, 512, 16000, fignu,'LP spectrum of sound /u/');

