% Proyecto Elecmetal
% Miguel A. Carrasco. (mlacarrasco@gmail.com)
% Analisis de regiones para la extracción de características
% v0.02. 02-05-2017  %analisis de multiples regiones

function data_analisis()
clear all;
close all;
%image directory
%bd_images= dir('data/region_move_*.png');
bd_images= dir('data/region_move_*.mat');
n = size(bd_images,1)
%for i=1:n
%    namefile = bd_images(i).name;
%    fprintf(' ID: %i -foto:%s\n', i,namefile);
%end

cont=1;
K=zeros(n,8);
IMA=[];
for id_ima=1:n
    namefile = bd_images(id_ima).name;
    str= sprintf('data/%s',namefile);
    load(str);
    %raw_image= double(imread(str));
    
    %bw= raw_image>0;
    %[I,J]= find(bw==1);
    %idx  = find(bw==1);
    %C= raw_image(idx)./255;
    %HA= hu(I,J,C);
    IMA(cont,:) = D(:)';
    K(cont,:) = dat;
    cont=cont+1;
    
end


figure, plot(K(:,1))
hold on
idx=find(K(:,end)==1);
idx_not = find(K(:,end)==0);
input= IMA;
input_k= K(:,1:end-1);
output = K(:,end);
DATA = [IMA K(:,end)];
[pcaDAT, input_PCA, ~, vSted] = pca_function(IMA, 0.99);
save('vSted.mat', 'vSted');
%save('inputs.mat', 'input_PCA');

training_nn(IMA', output');
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





