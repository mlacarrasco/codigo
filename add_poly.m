% Proyecto Elecmetal
% Miguel A. Carrasco. (mlacarrasco@gmail.com)
%v.0.01. Estructura que almance los vectores y predicciones
function M= add_poly(M, P, x, y, cont, pred_frame, xy_cross, max_frames)

frames =size(M.P,2);

%calculamos angulo de la linea
lx = x(end)-x(1);
ly = y(1)-y(end);
angle = atan(ly/lx)*180/pi;

if (frames<max_frames)
    M.P(:,frames+1) = [P'; cont; angle; pred_frame; xy_cross'];
else
    %agregamos una última columna al final y borramos la primera
    tmp = repmat(M.P,1,3);
    M.P(:,1:end-1)= tmp(:,frames+2:frames+max_frames);
    M.P(:,end) = [P'; cont; angle; pred_frame; xy_cross'];
    
end

end
