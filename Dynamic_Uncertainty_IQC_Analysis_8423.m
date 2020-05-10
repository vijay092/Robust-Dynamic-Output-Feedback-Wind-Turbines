
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
gamma_d = 1e-20;

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
del = 1e-3;

cvx_solver sedumi
cvx_precision best
cvx_begin sdp

    variable P(n,n) symmetric
    variable gamma_p_sq 
    variable lam 
    
     %LMI
      P >= del*eye(n)
      lam >= del
      gamma_p_sq >= del
      [P*A + A'*P,  P*B1, P*B3;...
      B1'*P,      -gamma_p_sq*eye(nu) ,    zeros(nu,n);...
      B3'*P ,     zeros(n,nu),   zeros(n)] +...
      lam*[C3 D31 D33; zeros(n)  zeros(n,nu) eye(n)]'*...
      [gamma_d*eye(n) zeros(n); zeros(n) -eye(n)] *...
      [C3 D31 D33; zeros(n) zeros(n,nu) eye(n)] + ...
      [C1'; D11'; D13']*[C1 D11 D13] <= -del*eye(2*n +nu)
  
cvx_end

