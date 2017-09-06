% Proyecto Elecmetal
% Miguel A. Carrasco. (mlacarrasco@gmail.com)
% Distancia entre un punto y un vector
% v0.01. 01-09-2017  

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

%velocidad de pixeles por frame
vel = mean(distance./time);

next_frame = ceil(min(distance)/vel);

regress_frame =  max(id_frames)+next_frame;
end