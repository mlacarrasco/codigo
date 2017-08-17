% UNIVERSIDAD DIEGO PORTALES
% Escuela de Informática y Telecomunicaciones
% Laboratorio de Vision por Computador (LVC)

% Demostracion de Momentos de Hu en Color
% (c) 2011. Miguel Carrasco (30-08-2011)  v.1.0

function test_momentos_color()

close all;
clear;
clc;

im=[ 4  0  0  2  2  3;
     2  1 12 15  5  4;
     1 11 25 34 56  2;
     3 16 34  3  4  1;
     2  0  4 33 65  3;
     4  1  6  3  4  2];


imA= round(imresize(im, 2));

% Segmentacion
bw= imA>8;
[I,J]= find(bw==1);
idx  = find(bw==1);

% Normalizacion
C= imA(idx)./255;
HA= hu(I,J,C);

% Resultados
fprintf('\n\n');
for i=1:7
    fprintf('Hu_A(%i): %e\n',i,HA(i));
    
end

%Ploteo
figure, imshow(imA,[])
end

%% FUNCIONES


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



