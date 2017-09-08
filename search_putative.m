% Proyecto Elecmetal
% Miguel A. Carrasco. (mlacarrasco@gmail.com)
%v.0.01. analisis de regiones candidatas (putativas)

function class = search_putative(POL, frame, g)

class= [];
dim =size(POL.P,2);
ax= axis();xa= ax(1):ax(2);


if dim>0
    id_projections = POL.P(6,:); %todos las proyecciones
    id = id_projections==frame;
    
    
    for i=1:dim
        try
        if (id(i)==1)
            tmp = POL.P(:,i);
            PL = tmp(1:2);
            xy_cross=tmp(end-1:end);
            yfit_corte=PL(1)*xa'+PL(2);
            hold on; plot(xa, yfit_corte,'r-', 'lineWidth',0.5);drawnow;
            plot(xy_cross(1), xy_cross(2),'gs','markerSize',15, 'lineWidth',3);
            fprintf('Bola detected.\y');
        end
        catch 
            class=[];
        end
        
    end
    
    
end



end