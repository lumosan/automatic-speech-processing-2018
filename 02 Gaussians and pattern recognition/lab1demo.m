% LAB4DEMO Demo for session 4 of Speech Processing labs
%
%    This script executes the exercises described in
%    the lab package available from:
%        ftp.idiap.ch/pub/sacha/labs/lab4.tgz
%

clear all; close all;
colordef none;
clc; echo on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Samples from a Gaussian density
N = 10000;

%%%%%%%%
%% Simulate the samples for each process

mu = [730 1090];

sigma_1 = [8000 0; 0 8000];
X1 = randn(N,2) * sqrtm(sigma_1) + repmat(mu,N,1);

sigma_2 = [8000 0; 0 18500];
X2 = randn(N,2) * sqrtm(sigma_2) + repmat(mu,N,1);

sigma_3 = [8000 8400; 8400 18500];
X3 = randn(N,2) * sqrtm(sigma_3) + repmat(mu,N,1);

%%%%%%%%
%% Plot these samples and the corresponding densities

hf1 = gausview(X1,mu,sigma_1,'Sample X1');
hf2 = gausview(X2,mu,sigma_2,'Sample X1');
hf3 = gausview(X3,mu,sigma_3,'Sample X1');

pause; % Press a key...    

clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Gaussian modeling: estimation of the mean and variance of a sample

%% Mean and variance of the whole sample (10000 points)
X = X3; N = size(X,1);
mu_10000 = sum(X)/N
sigma_10000 = (X - repmat(mu_10000,N,1))' * (X - repmat(mu_10000,N,1)) / (N-1)

%% Estimation error:
e_mu =  sqrt( (mu_10000 - mu) * (mu_10000 - mu)' )
% (This is the Euclidean distance between mu_10000 and mu_3)

e_sigma = norm( sigma_10000 - sigma_3 )
% (This is the 2-norm of the difference between sigma_10000 and sigma_3)


pause; % Press a key...    

clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Gaussian modeling: estimation of the mean and variance of a sample

%% Mean and variance with 1000 points
X = X3(1:1000,:); N = size(X,1);
mu_1000 = sum(X)/N
sigma_1000 = (X - repmat(mu_1000,N,1))' * (X - repmat(mu_1000,N,1)) / (N-1)

%% Estimation error
e_mu =  sqrt( (mu_1000 - mu) * (mu_1000 - mu)' )
% (This is the Euclidean distance between mu_1000 and mu_3)

e_sigma = norm( sigma_1000 - sigma_3 )
% (This is the 2-norm of the difference between sigma_1000 and sigma_3)


pause; % Press a key...    

clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Gaussian modeling: estimation of the mean and variance of a sample

%% Mean and variance with 100 points
X = X3(1:100,:); N = size(X,1);
mu_100 = sum(X)/N
sigma_100 = (X - repmat(mu_100,N,1))' * (X - repmat(mu_100,N,1)) / (N-1)

%% Estimation error
e_mu =  sqrt( (mu_100 - mu) * (mu_100 - mu)' )
% (This is the Euclidean distance between mu_100 and mu_3)

e_sigma = norm( sigma_100 - sigma_3 )
% (This is the 2-norm of the difference between sigma_100 and sigma_3)


pause; % Press a key...    


clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Likelihood of a sample with respect to a Gaussian model

% N = size(X3,1)
% logLike = 0;
% for i = 1:N;
%   logLike = logLike + (X3(i,:) - mu)*inv(sigma)*(X3(i,:) - mu)';
% end;
% logLike =  - 0.5 * ( logLike + N*log(det(sigma)) + 2*N*log(2*pi) )

%% Model N1
mu_1 = [730 1090]; sigma_1 = [8000 0; 0 8000];
echo off;
logLike1 = gloglike(X3,mu_1,sigma_1);
disp(' ');
disp(['log p(X|N_1) = ' num2str(logLike1)]);
disp(' ');
echo on;

%% Model N2
mu_2 = [730 1090]; sigma_2 = [8000 0; 0 18500];
echo off;
logLike2 = gloglike(X3,mu_2,sigma_2);
disp(' ');
disp(['log p(X|N_2) = ' num2str(logLike2)]);
disp(' ');
echo on;

%% Model N3
mu_3 = [730 1090]; sigma_3 = [8000 8400; 8400 18500];
echo off;
logLike3 = gloglike(X3,mu_3,sigma_3);
disp(' ');
disp(['log p(X|N_3) = ' num2str(logLike3)]);
disp(' ');
echo on;

%% Model N4
mu_4 = [270 1690]; sigma_4= [8000 8400; 8400 18500];
echo off;
logLike4 = gloglike(X3,mu_4,sigma_4);
disp(' ');
disp(['log p(X|N_4) = ' num2str(logLike4)]);
disp(' ');
echo on;


pause; % Press a key...    

clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% A-priori class probabilities
clear all;
load vowels.mat; whos
Na = size(a,1); Ne = size(e,1); Ni = size(i,1); No = size(o,1); Ny = size(y,1);
N = Na + Ne + Ni + No + Ny;

Pa = Na/N

Pe = Ne/N

Pi = Ni/N

Po = No/N

Py = Ny/N


pause; % Press a key...

clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Gaussian modeling of classes
hf = figure;
plot(a(:,1),a(:,2),'y+');
hold on;
set(gca,'xlim',[0 3000],'ylim',[0 3000],'dataaspectratio',[1 1 1e-6]);
xlabel('F1'); ylabel('F2'); grid on;
plot(e(:,1),e(:,2),'r+');
plot(i(:,1),i(:,2),'g+');
plot(o(:,1),o(:,2),'b+');
plot(y(:,1),y(:,2),'m+');
legend('/a/','/e/','/i/','/o/','/y/');
title('Simulated F1-F2 data for various vowels');

% (DO NOT CLOSE the obtained figure.)

pause; % Press a key...

mu_a = mean(a)
sigma_a = cov(a)
plotgaus(mu_a,sigma_a,[0 1 1]);


mu_e = mean(e)
sigma_e = cov(e)
plotgaus(mu_e,sigma_e,[0 1 1]);


mu_i = mean(i)
sigma_i = cov(i)
plotgaus(mu_i,sigma_i,[0 1 1]);


mu_o = mean(o)
sigma_o = cov(o)
plotgaus(mu_o,sigma_o,[0 1 1]);


mu_y = mean(y)
sigma_y = cov(y)
plotgaus(mu_y,sigma_y,[0 1 1]);


shg; pause; % Press a key...    

clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Bayesian classification

point = [400,1800];

echo off;
plot(point(1), point(2), 'c*', 'markersize',10); shg;
disp(['log P(/a/|x) = ' num2str( gloglike(point,mu_a,sigma_a) + log(Pa) )]);
disp(['log P(/e/|x) = ' num2str( gloglike(point,mu_e,sigma_e) + log(Pe) )]);
disp(['log P(/i/|x) = ' num2str( gloglike(point,mu_i,sigma_i) + log(Pi) )]);
disp(['log P(/o/|x) = ' num2str( gloglike(point,mu_o,sigma_o) + log(Po) )]);
disp(['log P(/y/|x) = ' num2str( gloglike(point,mu_y,sigma_y) + log(Py) )]);

echo on;

pause; % Press a key... 

clc;

point = [400,1000];

echo off;
plot(point(1), point(2), 'c*', 'markersize',10); shg;
disp(['log P(/a/|x) = ' num2str( gloglike(point,mu_a,sigma_a) + log(Pa) )]);
disp(['log P(/e/|x) = ' num2str( gloglike(point,mu_e,sigma_e) + log(Pe) )]);
disp(['log P(/i/|x) = ' num2str( gloglike(point,mu_i,sigma_i) + log(Pi) )]);
disp(['log P(/o/|x) = ' num2str( gloglike(point,mu_o,sigma_o) + log(Po) )]);
disp(['log P(/y/|x) = ' num2str( gloglike(point,mu_y,sigma_y) + log(Py) )]);

echo on;

pause; % Press a key... 

clc;

point = [530,1000];

echo off;
plot(point(1), point(2), 'c*', 'markersize',10); shg;
disp(['log P(/a/|x) = ' num2str( gloglike(point,mu_a,sigma_a) + log(Pa) )]);
disp(['log P(/e/|x) = ' num2str( gloglike(point,mu_e,sigma_e) + log(Pe) )]);
disp(['log P(/i/|x) = ' num2str( gloglike(point,mu_i,sigma_i) + log(Pi) )]);
disp(['log P(/o/|x) = ' num2str( gloglike(point,mu_o,sigma_o) + log(Po) )]);
disp(['log P(/y/|x) = ' num2str( gloglike(point,mu_y,sigma_y) + log(Py) )]);

echo on;

pause; % Press a key... 

clc;

point = [600,1300];

echo off;
plot(point(1), point(2), 'c*', 'markersize',10); shg;
disp(['log P(/a/|x) = ' num2str( gloglike(point,mu_a,sigma_a) + log(Pa) )]);
disp(['log P(/e/|x) = ' num2str( gloglike(point,mu_e,sigma_e) + log(Pe) )]);
disp(['log P(/i/|x) = ' num2str( gloglike(point,mu_i,sigma_i) + log(Pi) )]);
disp(['log P(/o/|x) = ' num2str( gloglike(point,mu_o,sigma_o) + log(Po) )]);
disp(['log P(/y/|x) = ' num2str( gloglike(point,mu_y,sigma_y) + log(Py) )]);

echo on;

pause; % Press a key... 

clc;

point = [670,1300];

echo off;
plot(point(1), point(2), 'c*', 'markersize',10); shg;
disp(['log P(/a/|x) = ' num2str( gloglike(point,mu_a,sigma_a) + log(Pa) )]);
disp(['log P(/e/|x) = ' num2str( gloglike(point,mu_e,sigma_e) + log(Pe) )]);
disp(['log P(/i/|x) = ' num2str( gloglike(point,mu_i,sigma_i) + log(Pi) )]);
disp(['log P(/o/|x) = ' num2str( gloglike(point,mu_o,sigma_o) + log(Po) )]);
disp(['log P(/y/|x) = ' num2str( gloglike(point,mu_y,sigma_y) + log(Py) )]);

echo on;

pause; % Press a key... 

clc;

point = [420,2500];

echo off;
plot(point(1), point(2), 'c*', 'markersize',10); shg;
disp(['log P(/a/|x) = ' num2str( gloglike(point,mu_a,sigma_a) + log(Pa) )]);
disp(['log P(/e/|x) = ' num2str( gloglike(point,mu_e,sigma_e) + log(Pe) )]);
disp(['log P(/i/|x) = ' num2str( gloglike(point,mu_i,sigma_i) + log(Pi) )]);
disp(['log P(/o/|x) = ' num2str( gloglike(point,mu_o,sigma_o) + log(Po) )]);
disp(['log P(/y/|x) = ' num2str( gloglike(point,mu_y,sigma_y) + log(Py) )]);

echo on;

pause; % Press a key... 

clc;

echo off;

