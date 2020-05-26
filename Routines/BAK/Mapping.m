function [X1,X2] = Mapping(x1,x2,L,R,courbe);
% Idem dans le code FUNAMBULE
% Je le sors juste pour abreger les sorties graphiques
ddx1 = min(diff(x1))/1000;

[Xc1 ,Xc2 ,dXc1dL ,dXc2dL ] = feval(courbe,L,x1     );
[XcP1,XcP2,dXc1dLp,dXc2dLp] = feval(courbe,L,x1+ddx1);

% Calcul de ses derivees numeriques en x1
dXc1dx1 = (XcP1-Xc1)/(ddx1);
dXc2dx1 = (XcP2-Xc2)/(ddx1);
ndXcdx1 = sqrt(dXc1dx1.*dXc1dx1 + dXc2dx1.*dXc2dx1);
% if min(ndXcdx1)==0;error('Point stationnaire');end
es1 = dXc1dx1./ndXcdx1; % vecteur tangeant
es2 = dXc2dx1./ndXcdx1;
er1 = es2;              % vecteur normal
er2 =-es1;

% Derivees secondes en dlambda dx1
d2Xc1dLdx1  = (dXc1dLp-dXc1dL)/ddx1;
d2Xc2dLdx1  = (dXc2dLp-dXc2dL)/ddx1;


% LES POINTS DE CALCUL --------------------------------------------


% Mapping : construction points X1 = Xc1 + R x2 er
X1 = ones(size(x2))*Xc1 + R*x2*er1;
X2 = ones(size(x2))*Xc2 + R*x2*er2;
