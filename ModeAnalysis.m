clear all;
close all;

load LMIMats
% The following matrices are generated from FAST V7
% for different azimuth angle (:,:,az,:) and different
% wind speed (:,:,:,wnd)

% What I wanted to examine is the stability of the system
% when the wind speed varies, i.e.,
% what if Ax is perturbed as (A + delta A)x
% This is essentially what happens when the wind speed changes.
% So, can we model this an uncertainty block?

wnd = 5;
idx = 1:31;


Ac = A;
B1c = B1;
B2c = B2;
C1c = C1;
C2c = C2;
D11c = D11;
D12c = D12;
D21c = D21;
D22c = D22;
% I have used the following notation wrt Prof. Seiler's notes
% Lecture 18
gamma_vec = 0;
gamma_d = 0.01;
for wnd = 3
    i = 1;
    for k = 1:30
        
        
        % xdot
        A = Ac(k:31,k:31);
        n = length(A);
        B2 = B2c(k:31,:);
        nd = size(B2,2);
        B1 = B1c(k:31,k:31);
        nw = size(B1,2);
        %B1(10) =1;
        
        % e
        C2 = C2c(k:31,k:31);
        ne = size(C2,1);
        D22 = D22c(k:31,:);
        D21 = D21c(k:31,k:31);
        
        % v
        nv = size(C1,1);
        C1 = C1c(:,k:31);
        D11 = D11c(:,k:31);
        D12 = D12c;
        del = 1e-2;
        
        cvx_solver sedumi
        cvx_precision low
        cvx_begin sdp
        
        variable P(n,n) symmetric
        variable gamma_p_sq
        variable lam
        
        %LMI
        P >= del*eye(n)
        lam >= del
        
        minimize(gamma_p_sq)
        [P*A + A'*P,  P*B1, P*B2;...
         B1'*P,      zeros(nw,nw) ,    zeros(nw,nd);...
         B2'*P ,     zeros(nd,nw),   -gamma_p_sq*eye(nd)] +...
         lam*[C1 D11 D12; zeros(nw,n) eye(nw) zeros(nw,nd)]'*...
         [gamma_d*eye(nv) zeros(nv,nw); zeros(nw,nv) -eye(nw)] *...
         [C1 D11 D12; zeros(nw,n) eye(nw) zeros(nw,nd)] + ...
         [C2'; D21'; D22']*[C2 D21 D22] <= -del*eye(n + nw + nd)
        
        cvx_end
        
        gam(i) = sqrt(gamma_p_sq)
        i = i+1;
    end
    figure(1);
    plot([1:30],gam,'m*-');
    
    xlabel('Modes $(k : 31)$','interpreter', 'latex')
    ylabel('$\gamma$','interpreter', 'latex')
    garyfyFigure
    hold on
end

legend('$v = 8 m/s$','interpreter', 'latex')
set(gca, 'FontName', 'Times New Roman')
box off;