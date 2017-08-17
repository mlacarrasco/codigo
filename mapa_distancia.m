% Proyecto Elecmetal
% Miguel A. Carrasco. (mlacarrasco@gmail.com)
% v.0.1.% Determina la distancia y clustering segun correlacion
% v.0.2% agregado a GIT

function SEL=mapa_distancia(pos, idx_i, idx_j, delta)


M=zeros(size(pos,2));
SEL=[];
for i=idx_i
    for j=idx_j
        if (j>i)
        A=  pos(:,i);
        B=  pos(:,j);
        M(i,j)= sqrt( (A(1)-B(1))^2+(A(2)-B(2))^2);
        end
        
    end
end

%https://stackoverflow.com/questions/37460067/cluster-distances-like-the-matlab-cluster-function-in-opencv
Z= linkage(M, 'ward');  %
K= cluster(Z, 'maxclust',2);
cl=hist(K, unique(K))
[val, id_min]= min(cl);
SEL = pos(:,K==id_min)

end
