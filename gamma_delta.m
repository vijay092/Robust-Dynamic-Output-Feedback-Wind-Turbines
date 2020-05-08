load LinMat



for i= 1:19
    g(i) = norm(Ag(:,:,1,1) - Ag(:,:,1,i+1),inf);
end

% plot(g);
% xlabel('Wind Speed (m/s)','interpreter', 'latex')
% ylabel('$H_{\infty}$ norm of $\Delta A$', 'interpreter', 'latex')
% garyfyFigure;