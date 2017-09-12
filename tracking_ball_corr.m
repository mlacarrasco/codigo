
% Proyecto Elecmetal
% Miguel A. Carrasco. (mlacarrasco@gmail.com)
% v0.01. 12-04-2017  %busqueda de regiones (correlación- template)
% v0.02. 02-05-2017  %analisis de multiples regiones
% v0.03. 05-06-2017  %Filtrado con clasificador neural
% v0.04. 08-06-2017  %analisis vectorial de movimiento. (Window)
% v0.05. 14-06-2017  %tracking de vectores por ángulo y distancia
% v0.06. 23-06-2017  %tracking en una zona y classificación en otra
% v0.07. 30-06-2017  %guardamos las regiones del tracking
% v0.08 10-07-2017   %analizamos las trayectorias de un mismo punto. Bloqueamos la que sean del mismo frame
% v0.09 14-07-2017   %empleamos la correlación entre dos señales para buscar similitud
% v0.10 14-08-2017   %borrado de codigo extra
% v0.11 17-08-2017   %analisis por bandas
% v0.12 17-08-2017 % codigo subido a GitHub https://github.com/mlacarrasco/codigo/
% v0.13 18-08-2017 % agregado la captura de vectores a una estructura
% v0.14 24-08-2017 % analisis de los vectores temporales
% v0.15 25-08-2017 % modificacion seccion de analisis de red neuronal
% v0.16 31-08-2017 % frame_idrol de Polinomio con RANSAC de 4 puntos (minimza el error del punto a la recta)
% v0.17 31-08-2017 % analizamos la trayectoria futura de un punto..para el clasificador
% v0.18 06-09-2017 % Interseccion con familia de polinomios
% v0.19 08-09-2017 % Buscamos en la posición proyectada


function track_ball()
close all;
clear all;
clc
warning off;


ratio         = 0.85;   % indice de correlacion de una region
factor        = 0.6;    % factor de escala del video

radio_px      = 6;    % radio bola (red neuronal entrenada con dicho tamaño)
factor_corr   = 0.95; % factor de correlacion entre imagenes
best_putatives= 5;    % numero de bolas candidatas maximo
draw          = 1;    %1: grafica
plot_lines    = 0;
max_frames    = 10;        % maximo numero de regiones
max_lines     = 200;  %máximo numero de lineas que cubre la zona

%parametros para grabar video
myVideo = VideoWriter('video_output/myfile_0.14.avi');
myVideo.FrameRate=10;
open(myVideo);
v= VideoReader('video_input/section_2.mov');
%%cargamos el archivo con el modelo de red neuroal
if exist('mat_files/net.mat')
    load('mat_files/net.mat');    %carga el modelo de red neuronal
    load('mat_files/vSted.mat');  %carga los vectores propios de PCA
end


T = imread('templates/template_2.png');
TR = imresize(T,factor);
BW= imresize(imread('templates/roi_area_detection_1280_banda.tif'),factor);

BW=BW>0;
tmp=regionprops(edge(BW),'PixelList');

%determinamos los perimetros de la mascara
perBW=[];
for i=1:size(tmp,1)
    perBW= [perBW; tmp(i).PixelList];
end

frame_id=1;

D_TRACK.data=[];     % almacena las regiones candidatas
D_POLY.P =[];

sw=0;

%POLINOMIOS
PL =[ 0.77701      -88.567 1]; %polinomio linea corte
PN_UP =  [-0.14883  176.15  1]; %polinomio movimiento superior
PN_DOWN= [-0.40257   456.4 1]; %polinomio movimiento inferior
POL(:,1)=linspace(PN_UP(1),PN_DOWN(1),max_lines);
POL(:,2)=linspace(PN_UP(2),PN_DOWN(2),max_lines);


while hasFrame(v)
    video = readFrame(v);
    g= video(:,:,2);
    gr = imresize(g, factor);
    roi_area= uint8(BW).*gr;
    %roi_area = g;

    msge= sprintf('frame id: %i\n',frame_id);
    imshow(gr);  hold on
    plot(perBW(:,1), perBW(:,2),'g.','MarkerSize',0.5); 
    text(size(gr,2)-90,size(gr,1)-10, msge, 'FontSize',14,'Color','green');
    axis on; drawnow;
    
    ax= axis();xa= ax(1):ax(2);
    yfit_corte=PL(1)*xa'+PL(2);
    plot(xa, yfit_corte,'r-', 'lineWidth',0.5);drawnow;
  
    
    %>>despliega las lineas por pantalla.
    if (plot_lines)
        for i=1:size(POL,1)
            yfit_normal = POL(i,1)*xa'+POL(i,2);
            plot(xa, yfit_normal,'y-', 'lineWidth',0.1);
        end
    end
    
    
    %>>>> TEMPLATE-MATCHING
    c= normxcorr2(TR, roi_area);
    
    %>> obtenemos un set de candidatos..
    %   los ordenamos de mayor a menor
    [values id] = sort(c(:), 'descend');
    id_sel= find(values(1:best_putatives)> ratio);
   
    %>>> SEARCH PUTATIVE 
    class = search_putative(D_POLY, frame_id, g, factor, PL, T);
    %>> analizamos si se encuentran regiones de interés
    if (not(isempty(id_sel)))
        
        values_sel = values(id_sel);
        
        for k_point=1:length(id_sel)
            
            [ypeak, xpeak] = find(c==values_sel(k_point));
            yoffSet = ypeak-size(T,1);
            xoffSet = xpeak-size(T,2);
            
            center= [xoffSet, yoffSet]+radio_px-1;
            D= xy_ray(roi_area, center, radio_px, frame_id, k_point, draw);
            
            
            %%>> RED NEURONAL
            %    tomamos la region y la convertimos en una columna
            %    evaluamos la región con una RED NEURONAL
            input=D(:)';%*vSted;
            outputs = net(input');
            
            %Si output es mayor a 0.1 significa que la región es putativa
            if (outputs>0.96)
                
                %>> incoramos esta info a la estructura temporal
                D_TRACK= add_region(D_TRACK, input, frame_id, max_frames, center);
                
                try
                    
                    res= corr(D_TRACK.data(1:end-3,:));
                    
                    %>> obtenemos las coordenadas
                    poss_xy= D_TRACK.data(end-1:end,:);
                    %>> obtenemos los indices de dichos frames
                    id_frame_id = D_TRACK.data(end-2,:);
                    
                    %>> buscamos solo las coordenadas que sean correspondientes
                    [px,~ ] =find(res>factor_corr & res < 1);
                    
                    %>> si no esta vacio significa que existe una relación de
                    %   similitud entre un patron y otro.
                    
                    if (not(isempty(px)) && length(px)>2)
                        
                        %>> analizamos las combinaciones de puntos
                        [SEL, id_frames, POL_SEL] =ransac_full(poss_xy, id_frame_id, px,POL, 0);
                        
                        %>> si obtengo más de dos puntos en correspondencias
                        if (size(SEL,2)>=4)
                            
                            hold on; plot(SEL(1,:), SEL(2,:), 'ys', 'markerSize',5,'LineWidth',2);
                            
                            %>> INTERSECCION
                            %seleccionamos el polinomio
                            P= [POL_SEL 1];
                            PT= cross(PL,P); %PL polinomio recta de borde
                            inter= PT./PT(2); %punto de interseccion entre recta
                            xy_cross=[inter(1), -inter(3)];
                            
                            
                            %>> ploteamos la linea de tendencia.
                            yfit=P(1)*xa'+P(2);
                            hold on;plot(xa, yfit,'m-.');drawnow;
                            plot(SEL(1,:), SEL(2,:), 'yx', 'markerSize',6)
                            plot(xy_cross(1), xy_cross(2),'gs','markerSize',15, 'lineWidth',3);
                            
                            %>> prediccion del frame
                            pred_frame = dist_pts(xy_cross',SEL,id_frames, frame_id , xa, yfit, radio_px)
                            
                            %>> agregamos el polinomio a un registro
                            D_POLY = add_poly(D_POLY, P, xa, yfit, frame_id, pred_frame, xy_cross, max_frames);
                            sw=1;
                        end
                    end
                    
                    
                catch
                    fprintf('waiting for more frames... \n')
                end
            end %fin output
            %
            
            if (frame_id>124)
                
                break;
            end
            
            %imrect(hAx, [xoffSet, yoffSet, size(T,2), size(T,1)]); drawnow;
            %F=getframe();
            
            % pause
            %imshow(J,'InitialMagnification',200);
            
            %s=sprintf('data/data_move_%i.png',frame_id);
            %imwrite(g,s);
            %end
        end%
        
        fprintf('%s',msge);
     
        frame_id= frame_id+1;
        
    end
    
    if (sw==0)
        F=getframe();
        writeVideo(myVideo,F.cdata);
    else
        F=getframe();
        for i=1:40
            writeVideo(myVideo,F.cdata);
        end
        sw=0;
    end
    
    
    if (frame_id>300)
        break;
    end
    
end
close(myVideo);

end




%% >>> FUNCIONES >>
function M= add_region(M, input, frame_id, max_frames, center)

frames =size(M.data,2);

if (frames<max_frames)
    M.data(:,frames+1) = [input, frame_id, center];
    
else
    %agregamos una última columna al final y borramos la primera
    tmp = repmat(M.data,1,3);
    M.data(:,1:end-1)= tmp(:,frames+2:frames+max_frames);
    M.data(:,end) = [input, frame_id, center];
    
end

end


%%
function D= xy_ray_pintar(center, radio_px, draw)

if (draw)
    angle=0:1:359;
    radian=deg2rad(angle);
    
    M=radio_px.*[sin(radian)', cos(radian)']+repmat(center,length(angle),1);
    CO=repmat(center,length(angle),1);
    xi=[CO(:,1), M(:,1)];
    yi = [CO(:,2), M(:,2)];
    
    hold on
    plot(xi, yi, 'g.'); drawnow;
    hold off
    
    
    D=[];
end
end

%%
function D= xy_ray(ima, center, radio_px, frame_id, k, draw)

angle=0:1:359;
radian=deg2rad(angle);

M=radio_px.*[sin(radian)', cos(radian)']+repmat(center,length(angle),1);
CO=repmat(center,length(angle),1);
xi=[CO(:,1), M(:,1)];
yi = [CO(:,2), M(:,2)];
D=[];

for t=1:length(angle)
    D(t,:)=improfile(ima, xi(t,:),yi(t,:), radio_px);
end

%bw= D>0;
%[I,J]= find(bw==1);
%idx  = find(bw==1);
%C= double(D(idx))./255;
%HA= hu(I,J,C);

if (draw)
    hold on
    plot(xi, yi, 'r.'); drawnow;
    hold off
end

sw=0;
if (frame_id>350 && sw)
    % Construct a questdlg with three options
    choice = questdlg('Es bola?', ...
        'Opciones', ...
        'Yes','No', 'No');
    % Handle response
    switch choice
        case 'Yes'
            disp([choice ' Es bola.'])
            hold on
            plot(xi, yi, 'r.'); drawnow;
            hold off
            option = 1;
        case 'No'
            disp([choice ' No es bola.'])
            option = 0;
    end
    %dat = [HA, option];
    
    s=sprintf('data/region_move_%i_%i.mat',frame_id,k);
    save(s,'D', 'dat');
end
%imwrite(uint8(D),s);
%figure, surf(D);
%shading interp;
end


