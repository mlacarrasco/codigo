
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
                     %Bloqueamos de sean del mismo frame
% v0.09 14-07-2017  %empleamos la correlación entre dos señales para buscar similitud





function track_ball()
close all;
clear all;
warning off;

ratio =0.90;
myVideo = VideoWriter('video_output/myfile_0.6.avi');
myVideo.FrameRate=10;

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
BW= imresize(imread('templates/roi_area_detection_1280_bandas.tif'),factor);

BW=BW>0;
tmp=regionprops(edge(BW),'PixelList');

perBW=[tmp.PixelList];

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
                        
                        for j= i:id_finish
                            
                            if (i~=j)
                                
                                %punto (x,y) de destino
                                xe= track(j,2);
                                ye= track(j,3);
                                
                                
                                id_end= find(sum(positions==repmat([xe,ye]',1,size(positions,2)))==2);
                                
                                if (draw)
                                    hold on
                                    plot(xi,yi, 'r*');
                                    plot(xe,ye, 'g*');
                                    hold off
                                end
                                
                                %analizamos la trayectoria de inicio con
                                %todos los destinos
                                DST(i,j)= sqrt((xi-xe)^2+(yi-ye)^2);
                                ANG(i,j)= atan2d(yi-ye,xi-xe);
                                
                                %diferencia entre dos vectores
                                %log(sqrt(sum( D_TRACK.data(1:end-3,id_ini)- D_TRACK.data(1:end-3,id_end)).^2));
                                DFF(i,j)= corr(D_TRACK.data(1:end-3,id_ini), D_TRACK.data(1:end-3,id_end));
                                
                                
                                %if (abs(ANG(i,j))<30 && DST(i,j)<5 && track(i)~=track(j))
                                if (abs(ANG(i,j))<190 &&  DST(i,j)<20 && DFF(i,j)>factor_corr && track(i)~=track(j))
                                    line([xi; xe],[yi; ye], 'LineWidth',2, 'Color','m'); drawnow;
                                    PTV(i,j)= 1;
                                    tmp = cross([xi yi 1], [xe, ye 1]);
                                    %tmp= normal_vector([xi, yi, xe, ye])
                                    MATRIX = MATRIX+ tmp;
                                else
                                    % marcamos un indice cero si la region 
                                    % no cumple los criterios anteriores
                                    PTV(i,j)= 0; 
                                end
                                
                            else
                                DST(i,j)=inf;
                                ANG(i,j)=inf;
                                DFF(i,j)=inf;
                            end
                            
                        end %end j
                      
                        
                        %borramos los valores con distancia cero
                        %DFF(i,find(DFF(i,:)==0))=inf;
                        
                        %dejamos en inf los puntos que no son putativos
                        DFF(i,find(PTV(i,:)==0))=inf; 
                        
                        %buscamos el mejor indice que sea tracking
                        [val_tmp, best_id] = min(DFF(i,:));
                      
                        if (not(isinf(val_tmp)))
                            PTV(i,best_id)= best_id;
                           xe_sel = track(best_id,2);
                           ye_sel = track(best_id,3);
                           PAIRS = [PAIRS; i, best_id];
                           line([xi; xe_sel],[yi; ye_sel], 'LineWidth',4, 'Color','g'); drawnow;
                           hold on; plot_line(MATRIX, xi, yi,size(BW,2), 'g'); drawnow
                        end
                       
                        
                    end
                catch
                    fprintf('waiting for more frames... \n')
                end
               
                
            end
            %
            
             if (cont==35)
                            
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
    
    F=getframe();
    writeVideo(myVideo,F.cdata);
    
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



function m=mrc(r,s,I,J,C)
% Funcion de momento R+S en color
i=I.^r;
j=J.^s;
m= sum(i.*j.*C);
end

function res=eta(r,s,I,J,C)
% Funcion eta (necesaria para momentos de Hu)
t= (r+s)/2 + 1;
a=  m_central(r,s,I,J,C);
b=  m_central(0,0,I,J,C);
res=a/(b^t);
end

function mc= m_central(r,s,I,J,C)
% Funcion de momento central
m00=mrc(0,0,I,J,C);
m10=mrc(1,0,I,J,C);
m01=mrc(0,1,I,J,C);

ci=m10/m00;
cj=m01/m00;

i= (I-ci).^r;
j= (J-cj).^s;

mc= sum(i.*j.*C);

end

function H=hu(I,J,C)
% Funcion de momentos de Hu
eta11=(eta(1,1,I,J,C));
eta12=(eta(1,2,I,J,C));
eta20=(eta(2,0,I,J,C));
eta21=(eta(2,1,I,J,C));
eta02=(eta(0,2,I,J,C));
eta03=(eta(0,3,I,J,C));
eta30=(eta(3,0,I,J,C));

H(1)=eta20+eta02;
H(2)=(eta20-eta02)^2+4*eta11^2;
H(3)=(eta30-3*eta12)^2+(3*eta21-eta03)^2;
H(4)=(eta30+eta12)^2+(eta21+eta03)^2;
H(5)=(eta30-3*eta12)*(eta30+eta12)*( (eta30+eta12)^2-3*(eta21+eta03)^2)+...
    (3*eta21-eta03)*(eta21+eta03)*(3* (eta30+eta12)^2- (eta21+eta03)^2 );
H(6)= (eta20-eta02)*...
    ((eta30+eta12)^2-(eta21+eta03)^2+ 4*eta11*(eta30+eta12)*(eta21+eta03));
H(7)= (3*eta21-eta03)*(eta30+eta12)*((eta30+eta12)^2-3*(eta21+eta03)^2)+...
    (eta30-3*eta12)*(eta21+eta03)* (3*(eta30+eta12)^2-(eta21+eta03)^2);
end


