%%% Universite Pierre et Marie Curie. UPMC
%%% Institute des Systemes Intelillents et Robotique
%%% Miguel Carrasco
%%% $ R.1.0 $ $02/05/2009$

function plot_line(V, x0, x1, dim_y, col)

%# Linear equations
a  = V(1); b = V(2); c = V(3);       
lx0 =  1;
ly0 = -(a*lx0 + c)/b;
lx1 =  dim_y;
ly1 = -(a*lx1 + c)/b;
line([lx0 lx1], [ly0 ly1], 'Color', col);

return