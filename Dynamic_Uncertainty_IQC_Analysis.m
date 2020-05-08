
clear all;
close all;

% The following matrices are generated from FAST V7
% for different azimuth angle (:,:,az,:) and different
% wind speed (:,:,:,wnd)
load LinMat

% What I wanted to examine is the stability of the system 
% when the wind speed varies, i.e.,
% what if Ax is perturbed as (A + delta A)x 
% This is essentially what happens when the wind speed changes.
% So, can we model this an uncertainty block? 

% Define the matrices;
az = 1; % select the azimuth angle you want to choose for the nom sys
wnd = 20; % select the azimuth angle you want to choose for the nom sys

% I have used the following notation wrt Prof. Seiler's notes
% Lecture 18

% xdot
A = Ag(:,:,az,wnd); 
n = length(A);
B2 = Bg(:,:,az,wnd);
nu = size(B2,2);
B1 = eye(n); 

% e
C2 = Cg(1,:,az,wnd);
D22 = Dg(1,:,az,wnd); 
D21 = zeros(1,n); 

% v
C1 = eye(n); 
D12 = zeros(n,nu); 
D11 =  zeros(n,n);
del = 1e-100;

cvx_solver sdpt3
cvx_precision best
cvx_begin sdp

    variable P(n,n) symmetric
    variable gamma_p_sq 
    variable lam 
    
     %LMI
      P >= del*eye(n)
      lam >= del
      
      [P*A + A'*P,  P*B1, P*B2;...
      B1'*P,      zeros(n) ,    zeros(n,nu);...
      B2'*P ,     zeros(nu,n),   -gamma_p_sq*eye(nu)] +...
      lam*[C1 D11 D12; zeros(n) eye(n) zeros(n,nu)]'*...
      [eye(n) zeros(n); zeros(n) -eye(n)] *...
      [C1 D11 D12; zeros(n) eye(n) zeros(n,nu)] + ...
      [C2'; D21'; D22']*[C2 D21 D22] <= -del*eye(2*n +nu)
  
cvx_end

