% Proyecto Elecmetal
% Miguel A. Carrasco. (mlacarrasco@gmail.com)
%v.0.01. Transforma la imagen en rayos

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