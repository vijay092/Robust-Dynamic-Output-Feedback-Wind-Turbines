
clear all;
close all;

% What I want to examine is the stability of the system 
% when the wind speed varies, i.e.,
% what if Ax is perturbed as (A + delta A)x 
% This is essentially what happens when the wind speed changes.
% So, can we model this an uncertainty block? 
sys =  rss(4);

% I have used the following notation wrt Prof. Seiler's notes
% Lecture 18
gamma_d = 1;

% xdot
A = sys.A; 
n = length(A);
B1 = sys.B;
nu = size(B1,2);
B3 = eye(n); 

% z
C1 = sys.C;
D13 = sys.D; 
D11 = ones(1,n); 

% v
C3 = eye(n); 
D33 = ones(n,nu);
D31 =  ones(n,n);
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