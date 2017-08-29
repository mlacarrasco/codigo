% Proyecto Elecmetal
% Miguel A. Carrasco. (mlacarrasco@gmail.com)
% v.0.1.% Determina la distancia y clustering segun correlacion
% v.0.2% agregado a GitHub
% v.0.3% analisis con kmeans iterativo.
% v.0.4% analisis vectorial con kmeans aglomerativo.


function SEL=mapa_distancia(pos, idx_i, idx_j, res, delta)

pts   = size(pos,2);
nodes =  size(idx_i,2);
M     = zeros(pts);
MC    = zeros(pts);

pos_sel = ones(2,pts)*inf;

for p=1:nodes
   
    x=idx_i(p);
    y=idx_j(p);
    
    if (x>y)
        A=  pos(:,x);
        B=  pos(:,y);
        pos_sel(:,x) = A;
        pos_sel(:,y) = B;
        %M(x,y)  = sqrt( (A(1)-B(1))^2+(A(2)-B(2))^2);
        %MC(x,y) = M(x,y)/res(x,y);
    end
        
   
end

X=[];
D=[];

for cluster=2:nodes
    
    %obtenemos los cluster variando el numero de nodos
    K = kmeans(pos_sel',cluster);
    
    %para cada subgrupo determinamos el error de proyeccion
    for j=1:max(K)
        SEL = pos_sel(:, K==j);
        idx= find(j==K);
        P= [polyfit(SEL(1,:), SEL(2,:), 1) 1];
      
       %ploteamos una recta que intersecta los puntos agrupados en el
       %subcluster.
        ax= axis(); xa= ax(1):ax(2); yfit=P(1)*xa'+P(2);
        hold on;plot(xa, yfit,'m-.');drawnow;
        
        distance=[];
        hold on;
        for p=1:size(SEL,2)
            %calculamos la distancia de un punto a la recta (polinomio)
              %forma del polinomio y= -a/b*x - c/b;
                distance(p) = (abs(-P(1)*SEL(1,p) + 1*SEL(2,p) - P(3))) / (sqrt(P(1)^2 + P(2)^2));
                plot(SEL(1,p),SEL(2,p),'ys');
        end
        %determinamos el error promedio de la distancia a la recta.
        %esto verifica cuan bien predice el cluster la distancia del punto
        %modelado con el polinomio
        err(j) = mean(distance);
        D(idx,cluster-1)=err(j);
    end
    
    %almacenamos los índices de los cluster
    X=[X, K];
    
    
    %determinamos el error
    
end

% OBJETIVO: Buscar el numero de cluster optimo que minimice la distancia
% hacia la recta.
pts_sort= sort(idx_i);

sel_D= D(pts_sort,:);
sel_X= X(pts_sort,:);


if(size(D,2)>2)
   
    %determinamos el histograma de nodos
    tmp= hist(sel_X, nodes);
  
    %calculamos los maximos por columnas (cada columna es un incremento en
    %el cluster). Clculamos desde la segunda columna en adelante.
   [val_max, row_tmp] =max(tmp(:,2:end));
   [~,col]=max(val_max);
   
   val_search=row_tmp(col); %maximo valor de histograma
   
   idx_ids=find(val_search ==sel_X(:,col+1))
   SEL = pos(:, pts_sort(idx_ids));
  
   
    %puntos seleccionados
end

%https://stackoverflow.com/questions/37460067/cluster-distances-like-the-matlab-cluster-function-in-opencv
%Z= linkage(MC, 'ward');  %
%K= cluster(Z, 'maxclust',2);
%K= clusterdata(Z,2);
%K= kmeans(Z,2);
%cl=hist(K, unique(K))
%[val, id_min]= min(cl);
%SEL = pos(:,K==id_min)

end
