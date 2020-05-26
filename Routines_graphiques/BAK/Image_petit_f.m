figure(5);clf;
set(5,'Name',['Image initiale f deformee x ',num2str(ampli)]);

% Image deformee
pcolor(X1+ampli*U1,X2+ampli*U2,F);
shading flat;colormap gray;
axis ij equal tight off;
hold on

% Maillage
if exist('XN1');
    % Champ de deplacement aux noeuds
    quiver(XNP1 ,XNP2 ,UNP1*ampli,UNP2*ampli,0,'w-','LineWidth',3)
    quiver(XNP1 ,XNP2 ,UNP1*ampli,UNP2*ampli,0,'k-','LineWidth',1)
    % Maillage initial
    plot(XNP1 ,XNP2 ,[cF,'-'],'LineWidth',trait_fin);
    plot(XNP1',XNP2',[cF,'-'],'LineWidth',trait_fin);
    % Maillage deforme
    plot(XNP1 +UNP1 *ampli,XNP2 +UNP2 *ampli,[cG,'-'],'LineWidth',trait_fin);
    plot(XNP1'+UNP1'*ampli,XNP2'+UNP2'*ampli,[cG,'-'],'LineWidth',trait_fin);
end

% Bords du domaine
plot(P1,P2,[cF,'-'],'LineWidth',trait_gro);
if exist('p1')
    plot(P1+ampli*UP1,P2+ampli*UP2,[cG,'-'],'LineWidth',trait_gro);
    quiver(P1,P2,ampli*UP1,ampli*UP2,0,'w-','LineWidth',3)
    quiver(P1,P2,ampli*UP1,ampli*UP2,0,'k-','LineWidth',1)
end
