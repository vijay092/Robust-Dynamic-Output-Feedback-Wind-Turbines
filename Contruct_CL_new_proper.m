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
K = zeros(1,n);
nw = n;
nd = 1;
nv = 1;
ne = n;


% Controller
K(:,21) = 2*kopt*wg_trim(3)*GBRatio;



A = Am(:,:,3) + Bm(:,1,3)*K;
B1 = zeros(n,nw);
B2 = Bdm(:,1,3);
C1 = K;
D11 = zeros(nv,nw);
D12 =  zeros(nv,nd);
C2 = eye(n);
D21 = eye(n);
D22 = zeros(ne,nd);

save('LMIMats','A','B1','B2','C1','C2','D11','D12','D21','D22')






