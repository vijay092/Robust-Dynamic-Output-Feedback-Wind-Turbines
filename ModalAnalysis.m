clear all;
close all;

load CLMats
% The following matrices are generated from FAST V7
% for different azimuth angle (:,:,az,:) and different
% wind speed (:,:,:,wnd)

% What I wanted to examine is the stability of the system
% when the wind speed varies, i.e.,
% what if Ax is perturbed as (A + delta A)x
% This is essentially what happens when the wind speed changes.
% So, can we model this an uncertainty block?

wnd = 5;
idx = 16:31;

% I have used the following notation wrt Prof. Seiler's notes
% Lecture 18
gamma_vec = 0;
gamma_d = 0.5;
wnd = 1;
k = 21;

% xdot
A = A_CL(k,k,wnd);
n = length(A);
B2 = B_CL(k,:,wnd);
nu = size(B2,2);
B1 = eye(n);

% e
C2 = C_CL(k,k,wnd);
ny = size(C2,1);
D22 = D_CL(k,:,wnd);
D21 = zeros(ny,n);

% v
C1 = zeros(1,n);
D12 = 0;
D11 =  1;
del = 1e-3;

cvx_solver sdpt3
cvx_precision high
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
    [gamma_d*eye(n) zeros(n); zeros(n) -eye(n)] *...
    [C1 D11 D12; zeros(n) eye(n) zeros(n,nu)] + ...
    [C2'; D21'; D22']*[C2 D21 D22] <= -del*eye(2*n +nu)

cvx_end

gam = sqrt(gamma_p_sq)

% figure(1);
% plot([1:15],gam,'*-');
% xlabel('Modes $(k : 16)$','interpreter', 'latex')
% ylabel('$\gamma$','interpreter', 'latex')
% garyfyFigure

% hold on


% legend('$v = 6 m/s$','$v = 7 m/s$','$v = 8 m/s$','$v = 9 m/s$','$v = 10 m/s$','$v = 11 m/s$','interpreter', 'latex')
% set(gca, 'FontName', 'Times New Roman')
% box off;