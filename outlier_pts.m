
% Proyecto Elecmetal
% Miguel A. Carrasco. (mlacarrasco@gmail.com)
% v0.01. 30-08-2017  %outlier basado en Leys

function idx= outlier_pts(err)
% buscamos un outlier
b= 1.4826;
mm = median(err(:));
MAD =b*median(abs(err-mm));
BIT = abs((err-mm)/MAD) > 3;
idx =find(BIT==1);
end
