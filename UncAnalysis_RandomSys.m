clear all;
close all;


% The following matrices are generated from FAST V7
% for different azimuth angle (:,:,az,:) and different
% wind speed (:,:,:,wnd)

% What I wanted to examine is the stability of the system
% when the wind speed varies, i.e.,
% what if Ax is perturbed as (A + delta A)x
% This is essentially what happens when the wind speed changes.
% So, can we model this an uncertainty block?


sys = rss(2);
gamma_vec = linspace(0, 0.4, 20);

for wnd = 3
    i = 1; gam = 0;
    for gamma_d = gamma_vec
        
        
        % xdot
        A = sys.A;
        n = length(A);
        B2 = sys.B;
        nd = size(B2,2);
        nw = 1;
        B1 = ones(n,nw);

        
        % e
        C2 = sys.C;
        ne = size(C2,1);
        D22 = sys.D;
        D21 = ones(ne,nw);
        
        % v
        nv = n;
        C1 = 0*eye(nv);
        D11 = zeros(nv,nw);
        D12 = ones(nv,nd);
        del = 1e-2;
        
        cvx_solver sedumi
        cvx_precision low
        cvx_begin sdp
        
        variable P(n,n) symmetric
        variable gamma_p_sq
        variable lam
        
        %LMI
        P >= del*eye(n);
        lam >= del;
        
        
        
        [P*A + A'*P,  P*B1, P*B2;...
            B1'*P,      zeros(nw,nw) ,    zeros(nw,nd);...
            B2'*P ,     zeros(nd,nw),   -gamma_p_sq*eye(nd)] +...
            lam*[C1 D11 D12;
            zeros(nw,n) eye(nw) zeros(nw,nd)]'*...
            [gamma_d^2*eye(nv) zeros(nv,nw); zeros(nw,nv) -eye(nw)] *...
            [C1 D11 D12; zeros(nw,n) eye(nw) zeros(nw,nd)] + ...
            [C2'; D21'; D22']*[C2 D21 D22] <= -del*eye(n + nw + nd);
        
        minimize(gamma_p_sq)
        cvx_end
        
        
         gamma(i) = sqrt(gamma_p_sq);
         gam = gamma(i);
         gam_plot(i) = gamma_d;
         i = i+1;

    end
    figure(1);
    plot(gam_plot,gamma,'m-');
    xlabel('$\gamma_\Delta$','interpreter', 'latex')
    ylabel('$\gamma$','interpreter', 'latex')
    garyfyFigure
end

legend('$8 m/s$','interpreter', 'latex')
set(gca, 'FontName', 'Times New Roman')
box off;