function [X1 X2] = DrawVirtualImage(F,L,R,courbe,x1,x2);
% Computes and plot the vitual image

global cF cG cD trait_gro trait_moy trait_fin myfontsize
    
ddx1 = min(diff(x1))/1000;
% Curve and derivatives / parameters L
[Xc1 ,Xc2 ,dXc1dL ,dXc2dL ] = feval(courbe,L,x1     );
[XcP1,XcP2,dXc1dLp,dXc2dLp] = feval(courbe,L,x1+ddx1);
% Numerical derivatives in x1
dXc1dx1 = (XcP1-Xc1)/(ddx1);
dXc2dx1 = (XcP2-Xc2)/(ddx1);
ndXcdx1 = sqrt(dXc1dx1.*dXc1dx1 + dXc2dx1.*dXc2dx1);
% Numerical second derivatives in x1
d2Xc1dLdx1  = (dXc1dLp-dXc1dL)/ddx1;
d2Xc2dLdx1  = (dXc2dLp-dXc2dL)/ddx1;
% Tangent vector es and normal vector er to the curve
es1 = dXc1dx1./ndXcdx1;
es2 = dXc2dx1./ndXcdx1;
er1 = es2;
er2 =-es1;
% Mapping : current points X = Xc1 (of the curve) + R x2 er
X1 = ones(size(x2))*Xc1 + R*x2*er1;
X2 = ones(size(x2))*Xc2 + R*x2*er2;

plot(Xc1,Xc2,[cG,'-'],'LineWidth',trait_gro);
% Borders of the virtual image
plot([X1(1  ,:);X1(end,:)]',[X2(1  ,:);X2(end,:)]',[cG,'-'],'LineWidth',trait_moy);
% Possible supplementary graphics
feval(courbe,L);
drawnow;