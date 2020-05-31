% Closed-loop matrices for the system:

clear all;
load NREL5MW_Data

% Consider only Region 2. The control law we use here
% is the standard kw^2 law for generator torque.


%% Linearized control input equation.
GBRatio = 97;

% Trim conditions from the .lin file
% This is rotor speed trim expressed in rad/s.
% These are for wind speeds [6,7....11]m/s
kopt = 2.3323;
wt_trim = [0.8378; 0.89012;0.942476;1.09956;1.20428;1.25664];
wg_trim = wt_trim*GBRatio;

% Linearized control
n = 31; % no. of states
K = zeros(1,n,6);
A_CL = zeros(n,n,6);
B_CL = zeros(n,1,6);
C_CL = zeros(n,n,6);
D_CL = zeros(n,1,6);

for i = 1:6
    % Controller
    K(:,21,i) = 2*kopt*wg_trim(i)*GBRatio;
    
    
    % Construct (A+BK)
    A_CL(:,:,i) = Am(:,:,i)+Bm(:,1,i)*K(:,:,i);
    B_CL(:,:,i) = Bdm(:,:,i);
    C_CL(:,:,i) = eye(n);
    D_CL(:,:,i) = zeros(n,1);
    
    
end

%save('CLMats','A_CL','B_CL','C_CL','D_CL')






