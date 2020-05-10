function [P, lam, gamma_p_sq] = dynUncertaintyLMI(A,B1,B3,C1,C3,D11,D13,D31,D33,gamma_d)

n = length(A);
nu = size(B,2);
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


end