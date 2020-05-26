% Figure 6 sans mouvement de corps rigide
figure(6);clf;
set(6,'Name',['Image initiale f deformee (sans mvt de corps rigide) x ',num2str(ampli)]);

% Image deformee
pcolor(X1+ampli*V1,X2+ampli*V2,F);
shading flat;colormap gray;
axis ij equal tight off;
hold on

% Maillage
if exist('XN1');
    % Champ de deplacement aux noeuds
    quiver(XNP1 ,XNP2 ,VNP1*ampli,VNP2*ampli,0,'w-','LineWidth',3)
    quiver(XNP1 ,XNP2 ,VNP1*ampli,VNP2*ampli,0,'k-','LineWidth',1)
    % Maillage initial
    plot(XNP1 ,XNP2 ,[cF,'-'],'LineWidth',trait_fin);
    plot(XNP1',XNP2',[cF,'-'],'LineWidth',trait_fin);
    % Maillage deforme
    plot(XNP1 +VNP1 *ampli,XNP2 +VNP2 *ampli,[cG,'-'],'LineWidth',trait_fin);
    plot(XNP1'+VNP1'*ampli,XNP2'+VNP2'*ampli,[cG,'-'],'LineWidth',trait_fin);
end

% Bords du domaine
plot(P1,P2,[cF,'-'],'LineWidth',trait_gro);
if exist('p1')
    plot(P1+ampli*VP1,P2+ampli*VP2,[cG,'-'],'LineWidth',trait_gro);
    quiver(P1,P2,ampli*VP1,ampli*VP2,0,'w-','LineWidth',3)
    quiver(P1,P2,ampli*VP1,ampli*VP2,0,'k-','LineWidth',1)
end
