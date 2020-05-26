function CARTES(F,DG,X1,X2,x1,x2,XN1,XN2)

% Elimination des points courants hors domaine
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
% disp('    Mouvement de corps rigide moyen :'); 
% disp(['    translation (',num2str(u1cs),',',...
% num2str(u2cs),')']);
% disp(['    angle ',num2str(theta_cs),' radians soit ',...
% num2str(theta_cs*180/pi),' degres']);
% Champ de deplacement sans mouvement de corps rigide moyen
V1 = U1 - (u1cs-theta_cs*X2);
V2 = U2 - (u2cs+theta_cs*X1);

% Calcul du tenseur des deformations de Green-Lagrange
F11 = gradient(x1 ) ;
F12 = gradient(x1')'; 
F21 = gradient(x2 ) ;
F22 = gradient(x2')';

% Green-Lagrange E = (F'*F - I)/2
E11 = (F11.^2 + F21.^2 - 1)/2;
E12 = (F11.*F12 + F21.*F22)/2;
%E21 = (F12.*F11 + F22.*F21)/2;
E22 = (F12.^2 + F22.^2 - 1)/2;

% Tenseur des petites deformations H=F-I, EPS = (H+H')/2
% Almansi-Euler A = (eye(2)-inv(FF*FF'))/2;
% Rotations
% FF = [F11 F12;F21 F22];
% R = FF*(FF'*FF)^(-0.5);
% R11 = R(1,1);R12=R(1,2):R21=R(2,1)
% Angle de rotation
% theta = angle(R(1,1)+i*R(1,2));

disp('======= Trace des cartes =======');
disp('    Deplacements : [U1], [U2], norme [NU] (pixels)')
disp('    Deplacement - deplacement de corps rigide : [V1], [V2], norme [NV] (pixels)');
%disp('    Gradient de la transformation F : [F11], [F12], [F21], [F22]');
%disp('    Tenseur des petites deformations EPS : [EPS11], [EPS12], [EPS21], [EPS22]');
disp('    Tenseur de Green-Lagrange E : [E11], [E12], [E22], norme [NE], trace [TE]');
disp('    Deviateur du tenseur de Green-Lagrange ED : [ED11], [ED12], [ED22], norme [NED]');
disp('    --------------------------------');
disp('    Axe 1 horizontal vers la droite');
disp('    Axe 2 vertical vers le bas.');

rep = 1;
while ~isempty(rep)
    commandwindow
    rep = input('    Composante ou [] pour sortir : ','s');
    if ~isempty(rep);     
        % Calcul de certaines composantes
        if strcmp(rep,'NU')
            carte = sqrt(U1.^2+U2.^2);
        elseif strcmp(rep,'NV')
            carte = sqrt(V1.^2+V2.^2);
        elseif strcmp(rep,'NE')
            carte = sqrt(E11.^2+2*E12.^2+E22.^2);
        elseif strcmp(rep,'TE')
            carte = E11+E22;
        elseif strcmp(rep,'ED11')
            carte = (E11-E22)/2;
        elseif strcmp(rep,'ED12')
            carte = E12;
        elseif strcmp(rep,'ED22')
            carte = (E22-E11)/2;
        elseif strcmp(rep,'NED')
            carte = sqrt((E11-E22).^2/2+2*E12.^2);
        else
            if ~exist(rep,'var')
                warning('This field does not exist');
                rep = 'U1';
            end
            carte = eval(rep);
        end
        % Trace de la carte
        figure(1);clf
        set(1,'Name',['Carte de ',rep,' sur image etat initial F']);
        h=imagesc(carte);
        colormap jet
        % Symetric colorscale - thus includes 0
        cax = caxis;
        caxis([-max(abs(cax)),max(abs(cax))]);
        masque = 1-F;
        masque(~DG) = F(~DG);       
        masque(isnan(carte))=0;
        set(h, 'AlphaData',masque);
        hold on;axis equal tight off;colorbar;set(gca,'FontSize',12);
        % Trace grille
        plot(XN1 ,XN2 ,'b-','LineWidth',1)
        plot(XN1',XN2','b-','LineWidth',1)
        %set(gca,'Position',[0.1 0.1 0.8 0.8]);
        commandwindow
    end
end
