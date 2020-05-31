clear all;
close all
load NREL5MW_Data
idx1 = {[1:31]};

Cm = eye(31);
Dm = zeros(31,1);
for i = 1
    
    idx = idx1{i};
    nom = 3;
    Gnom = ss(Am(:,:,nom),Bm(:,:,nom),Cm...
        ,Dm);
    for i = 1:6
        g(i) = norm(Gnom -...
            ss(Am(:,:,i),Bm(:,:,i),...
            Cm,Dm),inf);
    end
    
    wnd = 5 +[1:6];
    plot(wnd, g,'*-');
    hold on
    xlabel('Wind Speed (m/s)','interpreter', 'latex')
    ylabel(' $\|G_{nom} - G_i\|_\infty$', 'interpreter', 'latex')
    garyfyFigure;
    set(gca, 'FontName', 'Times New Roman')
    box off;
    
end


%legend('Model 1','Model 2', 'interpreter', 'latex')