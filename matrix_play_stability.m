% Construct a stable system of the 
clear all;

load LinMat
Ag(6,:,:,:) =[];
Ag(:,6,:,:) =[];
Bg(6,:,:,:) =[];
Cg(:,6,:,:) =[];

n = size(Ag,1);

B = zeros(n,2,36,20); D = zeros(2,2, 36, 20);
B(:,1,:,1:6) = Bg(:,:,:,1:6);
B(:,2,:,7:end) = Bg(:,:,:,7:end);
D(:,1,:,1:6) = Dg(:,:,:,1:6);   
D(:,2,:,7:end) = Dg(:,:,:,7:end);

Bg = B;
Dg = D;

max(real(eig(Ag(:,:,1,20))))


