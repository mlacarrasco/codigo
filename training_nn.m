% Proyecto Elecmetal
% Miguel A. Carrasco. (mlacarrasco@gmail.com)
% Analisis de regiones para la extracción de características
% v0.03. 05-06-2017  %Filtrado con clasificador neural


function [net]=training_nn(inputs, targets)


% Create a Pattern Recognition Network
hiddenLayerSize = 20;
net = patternnet(hiddenLayerSize);


% Set up Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;

% Train the Network
[net,tr] = train(net,inputs,targets);


save('mat_files/net.mat','net');
% Test the Network
outputs = net(inputs);
%errors = gsubtract(targets,outputs);
performance = perform(net,targets,outputs)
%figure, plotperform(tr)
figure, plotconfusion(targets,outputs)
% View the Network
%view(net)

end