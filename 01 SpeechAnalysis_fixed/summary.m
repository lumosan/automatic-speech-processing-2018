function [] = summary()
     close all
     data = load('f_s1_t1_a');
     data = data(17:length(data));
     stdata = data(15001:15960);
     glottal(stdata, 1);
