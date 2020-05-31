% Construct a stable system of the 
clear all;

load NREL5MW_Data
% Am(6,:,:) =[];
% Am(:,6,:) =[];
% Bm(6,:,:) =[];
% Bdm(6,:,:) =[];
% Cm(:,6,:) =[];

n = size(Am,1);

B = zeros(n,2,20); D = zeros(2,2,20);

B(:,1,1:6) = Bm(:,:,1:6);
B(:,2,7:end) = Bm(:,:,7:end);
D(:,1,1:6) = Dm(:,:,1:6);   
D(:,2,7:end) = Dm(:,:,7:end);

Bm = B;
Dm = D;


% 
% Ag(1:15,:,:,:) =[];
% Ag(:,1:15,:,:)=[];
% Bg(1:15,:,:,:)=[];
% Cg(:,1:15,:,:)=[];
max(real(eig(Am(:,:,20))))


