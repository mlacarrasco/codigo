% Proyecto Elecmetal
% Miguel A. Carrasco. (mlacarrasco@gmail.com)
%v.0.01. analisis de regiones candidatas (putativas)

function [class sw] = search_putative(T, D_TRACK, POL, frame, gr, scale, PCUT, Template)

if exist('mat_files/net.mat')
    load('mat_files/net.mat');    %carga el modelo de red neuronal
    %load('mat_files/vSted.mat');  %carga los vectores propios de PCA
end

sw=0;
pant=0;
class= [];
dim =size(POL.P,2);
ax= axis();xa= ax(1):ax(2);
radio_px = 6;


if dim>0
    id_projections = POL.P(6,:); %todos las proyecciones
    frame_median= floor(median(id_projections));
    check=frame==frame_median;
    output.class=[];
    
    if check
        
        for i=1:dim
            try
                if (id_projections(i)==frame_median)
                    tmp = POL.P(:,i);
                    PL = tmp(1:2);  %polinomio
                    
                    
                    
                    xy_cross=tmp(end-1:end);
                    yfit_corte=PL(1)*xa'+PL(2);
                    hold on; plot(xa, yfit_corte,'r-', 'lineWidth',0.5);drawnow;
                    
                    
                    %fprintf('Bola detected.\y');
                    %codigo de busqueda de bola-->
                    
                    B= gr(xy_cross(2)-20:xy_cross(2)+20, xy_cross(1)-20:xy_cross(1)+20);
                    %figure, imshow(B,[]); colormap(jet);
                    
                    
                    c= normxcorr2(T, B);
                    
                    [values id] = sort(c(:), 'descend');
                    id_sel= find(values(1)> 0.9);
                    values_sel = values(id_sel);
                    [ypeak, xpeak] = find(c==values_sel(1));
                    pxy=[ xy_cross(1)-20+xpeak-size(T,1)/2 xy_cross(2)-20+ypeak-size(T,2)/2]
                    % plot(pxy(1), pxy(2),'gs','markerSize',15, 'lineWidth',3);
                    center= pxy;
                    
                    plot(center(1), center(2),'gs','markerSize',15, 'lineWidth',3);
                    D= xy_ray(gr, center, radio_px, 500 , 0 , 0);
                    input=D(:)';%*vSted;
                    
                    %si salida es 1. red neuronal clasifica como bola.
                    output.class(i)  = net(input');
                    output.poss(i,:) = center;
                    ouput.frame= frame;
                    sw=1;
                    
                    if (pant)
                        
                        figure,
                        imshow(g);
                        ax= axis();xa= ax(1):ax(2);
                        yfit_scaled = PL(1)*xa'+PL(2)/scale;
                        yfit_cut    = PCUT(1)*xa'+PCUT(2)/scale;
                        hold on; plot(xa, yfit_scaled,'r-', 'lineWidth',0.5);drawnow;
                        hold on; plot(xa, yfit_cut,'r-', 'lineWidth',0.5);drawnow;
                        
                        
                        % imagen original
                        PLC = [PCUT(1) PCUT(2)/scale 1];
                        PLL = [PL(1) PL(2)/scale 1];
                        PI =cross(PLL,PLC);
                        %coordenada de interseccion (imagen original)
                        point_inter =[ PI(1)/PI(2) -PI(3)/PI(2)];
                        plot(point_inter(1),point_inter(2), 'ys','lineWidth',0.5);
                    end
                    
                end
            catch
                class=[];
            end
            
        end
        
        class=output;
    end
end



end