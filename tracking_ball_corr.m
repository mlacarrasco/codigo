
% Proyecto Elecmetal
% Miguel A. Carrasco. (mlacarrasco@gmail.com)
% v0.01. 12-04-2017  %busqueda de regiones (correlación- template)
% v0.02. 02-05-2017  %analisis de multiples regiones
% v0.03. 05-06-2017  %Filtrado con clasificador neural
% v0.04. 08-06-2017  %analisis vectorial de movimiento. (Window)
% v0.05. 14-06-2017  %tracking de vectores por ángulo y distancia
% v0.06. 23-06-2017  %tracking en una zona y classificación en otra
% v0.07. 30-06-2017  %guardamos las regiones del tracking
% v0.08 10-07-2017   %analizamos las trayectorias de un mismo punto.
%Bloqueamos la que sean del mismo frame
% v0.09 14-07-2017   %empleamos la correlación entre dos señales para buscar similitud
% v0.10 17-08-2017   %analisis por bandas



function track_ball()
close all;
clear all;
warning off;

ratio =0.90;
myVideo = VideoWriter('video_output/myfile_0.8.avi');
myVideo.FrameRate=5;

factor=0.6;
open(myVideo);
v= VideoReader('video_input/section_2.mov');

%%cargamos el archivo con el modelo de red neuroal
if exist('mat_files/net.mat')
    load('mat_files/net.mat');    %carga el modelo de red neuronal
    load('mat_files/vSted.mat');  %carga los vectores propios de PCA
end


T = imread('templates/template_2.png');
T = imresize(T,factor);
BW= imresize(imread('templates/roi_area_detection_1280_banda.tif'),factor);

BW=BW>0;
tmp=regionprops(edge(BW),'PixelList');

%determinamos los perimetros de la mascara
perBW=[];
for i=1:size(tmp,1)
    perBW= [perBW; tmp(i).PixelList];
end

cont=1;

Data=[];
val=[];
dist=[];
track=[];

radio_px=6;          % radio bola
factor_corr = 0.95;  %factor de correlacion entre imagenes
window = 10;         % frames
best_putatives=5;    % numero de bolas candidatas maximo
draw=1;
D_TRACK.data=[];     % almacena las regiones candidatas
max_frames=40;        % maximo numero de regiones

PAIRS=[];
while hasFrame(v)
    video = readFrame(v);
    g = imresize(video(:,:,2), factor);
    roi_area= uint8(BW).*g;
    %roi_area = g;
    
    
    c= normxcorr2(T, roi_area);
    %%obtenemos un set de candidatos..
    %%los ordenamos de mayor a menor
    [values id] = sort(c(:), 'descend');
    id_sel= find(values(1:best_putatives)> ratio);
    
    
    imshow(g);
    
    hold on
    plot(perBW(:,1), perBW(:,2),'b.');
    axis on; drawnow;
    
    
    if (not(isempty(id_sel)))
        %pos_list= id(id_sel);
        values_sel = values(id_sel);
        
        for k_point=1:length(id_sel)
            
            [ypeak, xpeak] = find(c==values_sel(k_point));
            yoffSet = ypeak-size(T,1);
            xoffSet = xpeak-size(T,2);
            Data(cont,:) =[cont, yoffSet xoffSet];
            dist(cont)=inf;
            center= [xoffSet, yoffSet]+radio_px-1;
            D= xy_ray(g, center, radio_px, cont, k_point, draw);
            
            %evaluamos la región con una red neuronal
            input=D(:)';%*vSted;
            outputs = net(input');
            
            D_TRACK= add_region(D_TRACK, input, cont, max_frames, center);
            
            
            %Si output es mayor a 0.1 significa que la región es putativa
            if (outputs>0.96)
                
                %almacena una region candidata para un posterior analisis
                
                
                fprintf('output:  %g\n',outputs);
                
                %almacenamos las coordenadas de las regiones putativas
                track= [track; cont center];
                xy_ray_pintar(center, radio_px, draw);
                
                %revisamos los últimos xx estados
                %determino los ciclos que son únicos
                
                %[WARNING]: esto habría que mejorar.
                %            Potencial riego de limite de memoria
                
                
                % seleccionamos los indices unicos del listado de tracking
                ids = unique(track(:,1));
                try
                    %determinamos el id del ciclo con xx cuadros antes
                    cont_sel= ids(end- window);
                    
                    id_start = find(track(:,1)==cont_sel);
                    id_finish =  find(track(end,1)==track);
                    
                    %analizamos la distancia de recorrido entre todos los
                    %puntos
                    DST=[];
                    ANG=[];
                    PTV=[];
                    
                    MATRIX=zeros(1,3);
                    
                    res= corr(D_TRACK.data(1:end-3,:));
                    poss_xy= D_TRACK.data(end-1:end,:);
                    [px, py] =find(res>factor_corr & res < 1);
                    
                    %si no esta vacio significa que existe una relación de
                    %similitud entre un patron y otro.
                    
                    if (not(isempty(px)))
                        %El clustering primero determina la distancia entre solo
                        %los puntos que cumplen una condición de similaridad.
                        %Luego realizamos un clustering basado en similitud de
                        %distancia. (analisis de dendograma)
                        delta =1;
                        SEL= mapa_distancia(poss_xy, px', py', delta);
                        
                        hold on; plot(SEL(1,:), SEL(2,:), 'bo')
                        
                        %determinamos un polinomio de grado 1 para crear
                        %una linea de proyeccion del movimiento.
                        P=polyfit(SEL(1,:), SEL(2,:), 1)
                        
                        ax= axis();
                        xa= ax(1):ax(2);
                        yfit=P(1)*xa'+P(2)
                        
                        %ploteamos la linea de tendencia.
                        hold on;plot(xa, yfit,'g-.');drawnow;
                    end
                    
                    
                    
                    
                    for i=id_start(1):id_finish(end)
                        
                        %punto (x,y) inicial
                        xi= track(i,2);
                        yi= track(i,3);
                        
                        %buscamos la columna donde se encuentra esta
                        %seccion
                        positions= D_TRACK.data(end-1:end,:); %coordenadas X-Y
                        id_ini= find(sum(positions==repmat([xi,yi]',1,size(positions,2)))==2);
                        
                        % ¿cual es el mejor destino en el siguiente frame?
                        % criterio 1. angulo y distancia
                        % criterio 2. analisis de la forma en otras
                        % secuencias (tracking en frames previos)
                        
                        %una vez que tenemos la trayectoria de todos
                        %los puntos de interes, analizamos el matching
                        
                        
                        
                    end %end j
                    
                catch
                    fprintf('waiting for more frames... \n')
                end
                
                
            end
            %
            
            if (cont==36)
                
                fprintf('Bola');
            end
            
            %imrect(hAx, [xoffSet, yoffSet, size(T,2), size(T,1)]); drawnow;
            F=getframe();
            
            % pause
            %imshow(J,'InitialMagnification',200);
            
            %s=sprintf('data/data_move_%i.png',cont);
            %imwrite(g,s);
            %end
        end%
        fprintf('cont:%i \t val:%f\n',cont, length(values_sel));
        cont= cont+1;
        
    end
    
    %F=getframe();
    %writeVideo(myVideo,F.cdata);
    
    if (cont>300)
        break;
    end
    
end
close(myVideo);
figure,
plot(Data(:,2), Data(:,3), 'b+');

figure, plot(val);
figure, plot(dist);

end


%%
function M= add_region(M, input, cont, max_frames, center)

frames =size(M.data,2);

if (frames<max_frames)
    M.data(:,frames+1) = [input, cont, center];
    
else
    %agregamos una última columna al final y borramos la primera
    tmp = repmat(M.data,1,3);
    M.data(:,1:end-1)= tmp(:,frames+2:frames+max_frames);
    M.data(:,end) = [input, cont, center];
    
    
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
function D= xy_ray(ima, center, radio_px, cont, k, draw)

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

bw= D>0;
[I,J]= find(bw==1);
idx  = find(bw==1);
C= double(D(idx))./255;
HA= hu(I,J,C);

if (draw)
    hold on
    plot(xi, yi, 'b.'); drawnow;
    hold off
end

sw=0;
if (cont>350 && sw)
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
    dat = [HA, option];
    
    s=sprintf('data/region_move_%i_%i.mat',cont,k);
    save(s,'D', 'dat');
end
%imwrite(uint8(D),s);
%figure, surf(D);
%shading interp;
end

%%
function [V, idx, mu]=pca(x, pct, cols)

% n: filas
% m: Características o columnas

n=size(x,1);
m=size(x,1);

mu=mean(x,1);                %0o. Calcular la media por columnas

B= x-repmat(mu,n,1);         %1o.  Restar la media
C= 1/(n-1) * B'*B;           %2o.  Calcular la covarianza

[V, D]  = eig(C);            %3o. Calcular los valores y vectores propios

d       = diag(D);           %4to. Tomar los valores de la diagonal y ordernarlos
[t, idx]= sort(d,'descend'); %     en forma descendiente

e_total = sum(d);            %5to. Calcular la energia total y acumulada de los
%     valores propios
cumE= cumsum(d(idx))/e_total;

sel= find(cumE<=pct);     %6to. Buscar un numero de columnas de acuerdo
W= V(:,idx(1:cols));         %     al criterio de la energia

KLT= B*W;                    %7mo. Terminamos. El resultados es la transformacion
%     Karhunen-Loeve

% Datos transformados con perdida
LOSS= KLT*W'+repmat(mu,n,1);


end
