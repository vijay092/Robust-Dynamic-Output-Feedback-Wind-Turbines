
clear all;
close all;

% The following matrices are generated from FAST V7
% for different azimuth angle (:,:,az,:) and different
% wind speed (:,:,:,wnd)
load LinMat

% What I want to examine is the stability of the system 
% when the wind speed varies, i.e.,
% what if Ax is perturbed as (A + delta A)x 
% This is essentially what happens when the wind speed changes.
% So, can we model this an uncertainty block? 

% Define the matrices;
az = 13; % select the azimuth angle you want to choose for the nom sys
wnd = 20; % select the azimuth angle you want to choose for the nom sys

% I have used the following notation wrt Prof. Seiler's notes
% Lecture 18
gamma_d = 1e-10;

% xdot
A = Ag(:,:,az,wnd); 
n = length(A);
B1 = Bg(:,:,az,wnd);
nu = size(B1,2);
B3 = eye(n); 

% z
C1 = Cg(1,:,az,wnd);
D13 = Dg(1,:,az,wnd); 
D11 = zeros(1,n); 

% v
C3 = eye(n); 
D33 = zeros(n,nu);
D31 =  zeros(n,n);


[P, lam, gamma_p_sq] = dynUncertaintyLMI(A,B1,B3,C1,C3,D11,D13,D31,D33,gamma_d)
