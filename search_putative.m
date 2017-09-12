% Proyecto Elecmetal
% Miguel A. Carrasco. (mlacarrasco@gmail.com)
%v.0.01. analisis de regiones candidatas (putativas)

function class = search_putative(POL, frame, g, scale, PCUT, Template)

class= [];
dim =size(POL.P,2);
ax= axis();xa= ax(1):ax(2);
radio_px = 6;


if dim>0
    id_projections = POL.P(6,:); %todos las proyecciones
    id = id_projections==frame;
    
    
    for i=1:dim
        try
            if (id(i)==1)
                tmp = POL.P(:,i);
                PL = tmp(1:2);  %polinomio
                frame_id = PL(4);
                
                xy_cross=tmp(end-1:end);
                yfit_corte=PL(1)*xa'+PL(2);
                hold on; plot(xa, yfit_corte,'r-', 'lineWidth',0.5);drawnow;
                plot(xy_cross(1), xy_cross(2),'gs','markerSize',15, 'lineWidth',3);
                %fprintf('Bola detected.\y');
                %codigo de busqueda de bola-->
                center= xy_cross;
                
                  D= xy_ray(roi_area, center, radio_px, ~ , ~, 0);
                
                if (pant)
                    
                    figure,
                    imshow(g);
                    ax= axis();xa= ax(1):ax(2);
                    yfit_scaled = PL(1)*xa'+PL(2)/scale;
                    yfit_cut    = PCUT(1)*xa'+PCUT(2)/scale;
                    hold on; plot(xa, yfit_scaled,'r-', 'lineWidth',0.5);drawnow;
                    hold on; plot(xa, yfit_cut,'r-', 'lineWidth',0.5);drawnow;
                end
                
                % imagen original
                PLC = [PCUT(1) PCUT(2)/scale 1];
                PLL = [PL(1) PL(2)/scale 1];
                PI =cross(PLL,PLC);
                %coordenada de interseccion (imagen original)
                point_inter =[ PI(1)/PI(2) -PI(3)/PI(2)];
                plot(point_inter(1),point_inter(2), 'ys','lineWidth',0.5);
                
                
            end
        catch
            class=[];
        end
        
    end
    
    
end



end