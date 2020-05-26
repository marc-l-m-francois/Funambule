function DEFORMEE(F,DG,X1,X2,x1,x2,P1,P2,p1,p2,XC1,XC2,xc1,xc2)
% Trace les deformees
% F image
% DG domaine de def correlation
% X1 X2 coo des pixels dans repere de F
% x1 x2 coo des pixels dans repere de g
% P1 P2 coo des points caracteristiques (bords ou maillage) rep F
% p1 p2 coo des points caracteristiques (bords ou maillage) rep g
% XC1 XC2 coo des lignes caracteristiques rep F
% xc1 xc2 coo des lignes caracteristiques rep g

% Elimination points hors domaine
x1(~DG) = NaN;
x2(~DG) = NaN;

% Calcul des deplacements
U1 = x1-X1;
U2 = x2-X2;

% Mouvement de corps rigide moyen. voir bibliographie "Calculs specifiques"
X1m = mean(X1(DG));
X2m = mean(X2(DG));
X12m = mean(X1(DG).^2);
X22m = mean(X2(DG).^2);
U1m = mean(U1(DG));
U2m = mean(U2(DG));
X1U2m = mean(X1(DG).*U2(DG));
X2U1m = mean(X2(DG).*U1(DG));
theta_cs = ( X1U2m-X1m*U2m-X2U1m+X2m*U1m) / (X12m+X22m-X1m^2-X2m^2);
u1cs = U1m+theta_cs*X2m;
u2cs = U2m-theta_cs*X1m;
disp('====== mouvement de corps rigide moyen ======');
disp(['    translation (',num2str(u1cs),',',...
num2str(u2cs),')']);
disp(['    angle ',num2str(theta_cs),' radians soit ',...
num2str(theta_cs*180/pi),' degres']);

% Champ de deplacement sans mouvement de corps rigide moyen
V1 = U1 - (u1cs-theta_cs*X2);
V2 = U2 - (u2cs+theta_cs*X1);

% Deplacement des points caracateristiques sans mouvement de corps rigide moyen
UP1 = p1-P1;
UP2 = p2-P2;
UP1CS = u1cs-theta_cs*P2;
UP2CS = u2cs+theta_cs*P1;
VP1 = UP1-UP1CS;
VP2 = UP2-UP2CS;

if ~isempty(XC1)
    % Deplacement des lignes caracateristiques sans mouvement de corps rigide moyen
    UC1 = xc1-XC1;
    UC2 = xc2-XC2;
    UC1CS = u1cs-theta_cs*XC2;
    UC2CS = u2cs+theta_cs*XC1;
    VC1 = UC1-UC1CS;
    VC2 = UC2-UC2CS;
end



% Choix
disp('======= Trace des deformees amplifees =======');
ampli = 1;
while ~isempty(ampli)
    
    % Deformee amplifiee
    figure(5);clf;set(gca,'Position',[0 0 1 1])
    set(5,'Name',['Image initiale f deformee x ',num2str(ampli)]);
    pcolor(X1+ampli*U1,X2+ampli*U2,F);caxis([0 1]);
    shading flat;colormap gray;
    axis ij equal tight off; hold on
    % Traces bords avant (bleu) arpes (rouge)
    plot(P1 ,P2 ,'b-','LineWidth',1);
    plot(P1',P2','b-','LineWidth',1);
    plot(P1 +ampli*UP1 ,P2 +ampli*UP2 ,'r+','LineWidth',1,'MarkerSize',24);
    plot(P1'+ampli*UP1',P2'+ampli*UP2','r+','LineWidth',1,'MarkerSize',24);
    quiver(P1,P2,ampli*UP1,ampli*UP2,0,'w-','LineWidth',3)
    quiver(P1,P2,ampli*UP1,ampli*UP2,0,'k-','LineWidth',1)
    if ~isempty(XC1)
        % Lignes caracteristiques avant (bleu) apres (rouge)
        plot(XC1,XC2,'b-','LineWidth',1);
        plot(XC1+ampli*UC1,XC2+ampli*UC2,'r-','LineWidth',1,'MarkerSize',24);
    end
    % Deformee amplifiee sans mouvement de corps rigide
    figure(6);clf;set(gca,'Position',[0 0 1 1])
    set(6,'Name',['Image initiale f deformee (sans mvt de corps rigide) x ',num2str(ampli)]);
    pcolor(X1+ampli*V1,X2+ampli*V2,F);caxis([0 1]);
    shading flat;colormap gray;
    axis ij equal tight off;hold on
    % Traces bords avant (bleu) apres (rouge)
    plot(P1 ,P2 ,'b-','LineWidth',1);
    plot(P1',P2','b-','LineWidth',1);
    plot(P1 +ampli*VP1 ,P2 +ampli*VP2 ,'r+','LineWidth',1,'MarkerSize',24);
    plot(P1'+ampli*VP1',P2'+ampli*VP2','r+','LineWidth',1,'MarkerSize',24);
    quiver(P1,P2,ampli*VP1,ampli*VP2,0,'w-','LineWidth',3)
    quiver(P1,P2,ampli*VP1,ampli*VP2,0,'k-','LineWidth',1)
    if ~isempty(XC1)
        % Lignes caracteristiques avant (bleu) apres (rouge)
        plot(XC1,XC2,'b-','LineWidth',1);
        plot(XC1+ampli*VC1,XC2+ampli*VC2,'r-','LineWidth',1,'MarkerSize',24);
    end
    % Autre choix amplification
    disp(['    Coefficient actuel d''amplification dela deformee ',num2str(ampli)]);
    ampli = input('        Autre valeur ? ou [] pour sortir ');
end
