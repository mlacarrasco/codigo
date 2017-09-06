
% Proyecto Elecmetal
% Miguel A. Carrasco. (mlacarrasco@gmail.com)
% v0.01. 30-08-2017  %buscamos una combinación con el menor error)
% v0.02. 06-09-2017  %busqueda de polinomio que intersecta

function [SEL, indices,POL_SEL]=ransac_full(pts, id_cont, id_sel, POL, pant)

SEL      = [];
best_cmb = [];
indices  = [];
backup   = pts;
POL_SEL  = [];

%id_sel: indices de puntos correlacionados
pts       = pts(:,id_sel);
id_frame  = id_cont(id_sel);

if (size(pts,2)>=4)
    pts = pts';
   [pts_s, ci] = unique(pts,'rows');
    pts_s = pts_s';
    nodes =size(pts_s,2);
    cmb =  nchoosek(1:nodes,4);
    indices=id_frame(ci);
    backup_ind= indices;
    pts= pts_s;
end
err   = zeros(size(cmb,1),1);
pol_s = zeros(size(cmb,1),1);

for i=1:size(cmb,1)
    
    cols= cmb(i,:);
    %coordenadas seleccionadas (combinatorial)
    SEL = pts(:, cols); 
    
    %P= [polyfit(SEL(1,:), SEL(2,:), 1) 1];
    tmp=[];
    %para cada punto
    for j=1:size(SEL,2)
        
        pt=[SEL(:,j)' 1];
        G=[];
        
        %calculamos la distancia del punto a todas las rectas.
        for k=1:size(POL,1)
            G(k) = (abs(POL(k,1)*pt(1) - 1*pt(2) + POL(k,2))) / (sqrt(POL(k,1)^2 + 1));
        end
        
        %buscamos la mejor recta.
        [val, id]=min(G);
        tmp = [tmp; val id];
        
    end
    s= sum(tmp,1);
    err(i) = s(1);
    pol_s(i) = mode(tmp(:,2));
    
    
%     %ploteamos una recta que intersecta los puntos agrupados en el
%     %subcluster.
%     if (pant)
%         ax= axis(); xa= ax(1):ax(2); yfit=P(1)*xa'+P(2);
%         hold on;plot(xa, yfit,'y-.');drawnow;
%     end
%     
%     distance = zeros(size(pts,2),1);
%     for p=1:size(pts,2)
%         %calculamos la distancia de un punto a la recta (polinomio)
%         %forma del polinomio y= -a/b*x - c/b;
%         distance(p) = (abs(-P(1)*pts(1,p) + 1*pts(2,p) - P(3))) / (sqrt(P(1)^2 + P(2)^2));
%         %plot(SEL(1,p),SEL(2,p),'ys','MarkerSize',10);
%     end
%     err(i) = mean(distance);
    
end

try
[~, id_sort]=sort(err);
best_cmb = cmb(id_sort(1),:)
SEL= pts(:,best_cmb);
indices = indices(best_cmb);
thickness = size(SEL,2);
POL_SEL = POL(pol_s(id_sort(1)),:);
catch
    SEL=[];
    best_cmb=[];
    thickness =0;
end

end



