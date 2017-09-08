% Proyecto Elecmetal
% Miguel A. Carrasco. (mlacarrasco@gmail.com)
% Distancia entre un punto y un vector
% v0.01. 01-09-2017  
% v0.02. 07-09-2017  %alicamos un modelo de regresion


function regress_frame = dist_pts(point, V, id_frames, now, x,y, radio)

nodes=size(V,2);
P=repmat(point,1, nodes);

%calculamos angulo de la linea
lx = x(end)-x(1);
ly = y(1)-y(end);
angle = atan(ly/lx)*180/pi;

%distancia entre el punto proyectado y los puntos del vector
distance = sqrt(sum((P-V).^2))+radio*2;
time     = abs(id_frames-now)+1;

x=[ones(size(distance,2),1), distance'];
y = id_frames';

xfit= 0:1:max(x(:,2));
model =regress (y,x);
yfit = model(1) + model(2)*xfit;

%figure, plot(y(:), x(:,2)); hold on
%plot(yfit, xfit)

%frame de posicion proyectada
regress_frame =ceil(model(1));
end