%codigo descartado. 

function [M, POSS]= malla_analisis(T, D_TRACK, PL, PCUT,xa,gr, radio_px)

if exist('mat_files/net.mat')
    load('mat_files/net.mat');    %carga el modelo de red neuronal
    %load('mat_files/vSted.mat');  %carga los vectores propios de PCA
end

M=[];
POSS=[];
range_i= -10:4:10;
range_j= 1:4:20;

REG= D_TRACK.data(1:end-3,:);
ci=1; cj=1;
for t=range_i
    yfit_pol=PL(1)*xa'+PL(2)+t;
    PLL = [PL(1) PL(2)+t 1];
    
    for k=range_j
        yfit_corte=PCUT(1)*xa'+PCUT(2)+k;
        hold on; plot(xa, yfit_pol,'g-', 'lineWidth',0.5);drawnow;
        hold on; plot(xa, yfit_corte,'r-', 'lineWidth',0.5);drawnow;
        
        PLC = [PCUT(1) PCUT(2)+k 1];
        PI =cross(PLL,PLC);
        point_inter =[ PI(1)/PI(2) -PI(3)/PI(2)];
        plot(point_inter(1),point_inter(2), 'rs','lineWidth',0.5);
        
        center=ceil(point_inter);
        %patch
      
        
        D= xy_ray(gr, center, radio_px, 500 , 0 , 0);
            
        input=D(:)';%*vSted;
        res= corr([REG, input']);
        res>0.95
        %si salida es 1. red neuronal clasifica como bola.
        M(cj+(ci-1)*length(range_i),:)  = net(input');
        POSS(cj+(ci-1)*length(range_i),:) = center;
        
        sw=1;
        
        cj=cj+1;
    end
    ci=ci+1;
    
end

id=find(M>0.96);
hold on;
plot(POSS(id,1), POSS(id,2), 'bs');
