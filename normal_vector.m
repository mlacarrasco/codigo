%%% Universite Pierre et Marie Curie. UPMC
%%% Institute des Systemes Intelillents et Robotique
%%% Miguel Carrasco  miguel.carrasco@etu.upmc.fr
%%% $ R.2.0 $ $05/05/2009$

% Computes the Normal vector to last point

function NV=normal_vector(V)

xa=V(1);ya=V(2);
xb=V(3);yb=V(4);

NV= [ xa-xb;
       ya-yb;
       xb*(xb-xa)+ yb*(yb-ya)];
   
return