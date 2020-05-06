clear all;
close all;


load LinMat



for i= 1:19
    g(i) = norm(Bg(:,:,1,1) - Bg(:,:,1,i+1),'fro')
end

plot(g)