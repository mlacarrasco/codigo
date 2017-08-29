% Proyecto Elecmetal
% Miguel A. Carrasco. (mlacarrasco@gmail.com)
% v.0.1.% Determina la distancia y clustering segun correlacion
% v.0.2% agregado a GitHub

function SEL=mapa_distancia(pos, idx_i, idx_j, res, delta)

pts = size(pos,2);
M=zeros(pts);
MC=zeros(pts);
SEL=[];
for i=idx_i
    for j=idx_j
        %if (j>i)
        A=  pos(:,i);
        B=  pos(:,j);
        M(i,j)  = sqrt( (A(1)-B(1))^2+(A(2)-B(2))^2);
        MC(i,j) = M(i,j)/res(i,j);
        %end
        
    end
end

X=[];
D=[];

for cluster=1:pts
    
    %obtenemos los cluster
    K = kmeans(M,cluster);
    
    %para cada subgrupo determinamos el error de proyeccion
    for j=1:max(K)
        SEL = pos(:, K==j);
        idx= find(j==K);
        P= [polyfit(SEL(1,:), SEL(2,:), 1) 1];
        %forma del polinomio y= -a/b*x - c/b;
        ax= axis(); xa= ax(1):ax(2); yfit=P(1)*xa'+P(2);
        hold on;plot(xa, yfit,'m-.');drawnow;
        
        distance=[];
        hold on;
        for p=1:size(SEL,2)
          distance(p) = (abs(-P(1)*SEL(1,p) + 1*SEL(2,p) - P(3))) / (sqrt(P(1)^2 + P(2)^2));
          plot(SEL(1,p),SEL(2,p),'ys');
        end
        err(j) = mean(distance);
        D(idx,cluster)=err(j);
    end
    
    X=[X, K];
    
    
    %determinamos el error
    
end
%https://stackoverflow.com/questions/37460067/cluster-distances-like-the-matlab-cluster-function-in-opencv
Z= linkage(MC, 'ward');  %
%K= cluster(Z, 'maxclust',2);
K= clusterdata(Z,2);
K= kmeans(Z,2);
cl=hist(K, unique(K))
[val, id_min]= min(cl);
SEL = pos(:,K==id_min)

end
