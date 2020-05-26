disp(' ');
disp('       -----------------------------------------------------');
disp('       |2016                                           2019|')
disp('       |                                                   |')
disp('       |                   DigImCo_Q4                      |')
disp('       |    DIGital IMage COrrelation elements finis Q4    |');
disp('       |        Identification par Images Numeriques       |');
disp('       |                                                   |')
disp('       |Marc Francois                Labo. GeM Univ. Nantes|')
disp('       |marc.francois@univ-nantes.fr                       |');
disp('       -----------------------------------------------------');
disp(' ');


% INITS MATLAB ------------------------------------------------------------


clear all;
close all;
path(path,'Routines');
path(path,'Routines_graphiques');
path(path,'Champs');
% Detection de la toolbox image
VersionMatlab = ver;
ToolBoxImages = false;
for p=1:size(VersionMatlab,2)
    if strcmp(VersionMatlab(p).Name,'Image Processing Toolbox')
        ToolBoxImages = true;
    end
end


% INITS LOGICIEL ----------------------------------------------------------


% Type d'interpollation retenu
type_interp	= 'linear'; 
% Activer pour verifier les gradients des champs (debugage)
verif_grad = false;
% Nombre maxi d'iterations
cptmax     = 100;
% Limite de variation des li
precision  = 1E-3;


% INITS GRAPHIQUES --------------------------------------------------------


cg = 'r';   % Couleur Couleur bord du domaine
cF = 'b';   % Couleur maillage
myfontsize = 18;


% CHARGEMENT IMAGES -------------------------------------------------------


[nom_F,ImagePath] = uigetfile('*.*','Choisir l''image de l''etat initial');
IMA = imread(fullfile(ImagePath,nom_F));
F = PrepImage(IMA);
[nom_g,ImagePath] = uigetfile('*.*','Choisir l''image de l''etat final');
IMA = imread(fullfile(ImagePath,nom_g));
g = PrepImage(IMA);
clear IMA;
% Sauvegarde des images (pour les reprendre apres le flou)
save Temporaire/F0 F
save Temporaire/g0 g
% Coordonnees des pixels en repere de F
LX1 = [1:size(F,2)];
LX2 = [1:size(F,1)];
[X1,X2] = meshgrid(LX1,LX2);


% INITS SORTIES GRAPHIQUES ------------------------------------------------


% Image F
figure(1);clf;Range_mon_trace(1,1,1,2,4);
set(1,'Name','Image etat initial F','MenuBar','none',...
    'GraphicsSmoothing','off','ToolBar','none');
PlotImage(F,1);

% Image g
figure(2);Range_mon_trace(2,1,2,2,4);
set(2,'Name','Image etat final g','MenuBar','none',...
    'GraphicsSmoothing','off','ToolBar','none');
PlotImage(g,2);

% Image |G-F|
figure(3);Range_mon_trace(3,2,1,2,4);

% Convergence
figure(4);Range_mon_trace(4,2,2,2,4);clf;
set(4,'Name','Convergence','MenuBar','none',...
    'GraphicsSmoothing','off','ToolBar','none');

% Deformee
figure(5);Range_mon_trace(5,1,3,2,4);clf;
% set(5,'MenuBar','none','GraphicsSmoothing','off','ToolBar','none');
% 
% Deformee sans corps rigide
figure(6);Range_mon_trace(6,2,3,2,4);clf;
% set(6, 'MenuBar','none','GraphicsSmoothing','off','ToolBar','none');


% DOMAINE DE CORRELATION --------------------------------------------------


[P1,P2,DF] = ChoixDomaineCorrelation(F,X1,X2);


% ETAPES INITIALES AVEC FLOU -----------------------------------------------


if ToolBoxImages
    disp('    Un calcul initial sur une image floutee aide la convergence en cas de grands deplacements');
    rayon_de_flou = input('    Entrer le ou les niveaux de flou souhaite par, ex. [12 6 0], ou [0] si non : ');
else
    rayon_de_flou = 0;
end
if rayon_de_flou(end)~=0;rayon_de_flou=[rayon_de_flou,0];end


% PHASE 1 CALCUL DEFORMATION HOMOGENE -------------------------------------


disp('======= Calcul de la deformation homogene');

% Init des parametres
nom_champ = 'DEFHOM';
li = [0 0 0 0 0 0];

for rg = rayon_de_flou
    % Chargement images brutes
    load Temporaire/F0;
    load Temporaire/g0;
    % Application du flou gaussien
    if rg>0
        disp(['Etape intiale avec un flou de rayon ',num2str(rg)]);
        H = fspecial('gaussian',6*ceil(rg)+1,rg/2);
        F = filter2(H,F);
        g = filter2(H,g);
    else
        disp('Etape finale avec les images brutes');
    end

    
	% DISPLAYS IMAGES EN COURS DE TRAITEMENT ------------------------------


    PlotImage(F,1,P1,P2);
    PlotImage(g,2);

    
    % INITS CORRELATION ---------------------------------------------------


    % Gradients de l'image F
    [gradF10,gradF20]  = gradient(F);
    % Init compteur de boucles et ecart
    cpt = 1;
    % Champ actuel
    [x1,x2] = feval(nom_champ,'champ',li,X1,X2);
    % Image G = image g interpolee
    G = interp2(g,x1,x2,type_interp,NaN);
    % Son domaine de definition intersection avec celui de F
    DG = isfinite(G) & DF;
    % Surface = nombre de pixels correles
    S = sum(sum(DG));
    % Ecart GG-F
    G_F = G-F;
    % On annule l'effet des points hors domaine
    G_F(~DG) = 0;
    % Erreur actuelle
    Psi = sum(sum(G_F.^2))/S;
    % Variables internes au Newton
    % gradF_dxdl = zeros(size(F,1),size(F,2),size(li,1));
    gradF_dxdlp = zeros(size(F,1),size(F,2));
    gradF_dxdlq = zeros(size(F,1),size(F,2));
    M = zeros(max(size(li)));
    V = zeros(max(size(li)),1);


    % NEWTON --------------------------------------------------------------


    % Init correcteur pour rentrer dans la boucle
    dli = inf*ones(size(li));

    while max(abs(dli))>precision && cpt<cptmax

        cpt = cpt+1;
        % Elimination des points hors domaine
        gradF1 = gradF10;
        gradF2 = gradF20;
        gradF1(~DG) = 0;
        gradF2(~DG) = 0;
        % Calcul matrice et vecteur
        for p=1:max(size(li))
             % Calcul de la p-ieme derivee en li
            [dx1dL,dx2dL] = feval(nom_champ,'derivees',li,X1,X2,p);
            gradF_dxdlp  =  gradF1.*dx1dL +  gradF2.*dx2dL;
            % Vecteur
            V(p) = -sum(sum(gradF_dxdlp.*G_F));
            % Matrice
            for q=1:p
                if q~=p
                    [dx1dL,dx2dL] = feval(nom_champ,'derivees',li,X1,X2,q);
                    gradF_dxdlq  =  gradF1.*dx1dL +  gradF2.*dx2dL;
                    M(p,q) = sum(sum( gradF_dxdlp.*gradF_dxdlq ));
                    M(q,p) = M(p,q);
                else
                    M(p,q) = sum(sum( gradF_dxdlp.*gradF_dxdlp));
                end
            end
        end

        % Verification de la determination de la matrice
        detM = abs(det(M));
        if detM<1E-6
            warning(['      La matrice a determinant ',num2str(detM),' presque nul)']);
        end       

        % Correcteur
        dli = (M\V)';
        li = li + dli;

        % Champ actuel
        [x1,x2] = feval(nom_champ,'champ',li,X1,X2);
        % Image G = image g interpolee
        G = interp2(g,x1,x2,type_interp,NaN);
        % Son domaine de definition intersection avec celui de F
        DG = isfinite(G) & DF;
        % Surface = nombre de pixels correles
        S = sum(sum(DG));
        % Ecart GG-F
        G_F = G-F;
        % On annule l'effet des points hors domaine
        G_F(~DG) = 0;
        % Erreur actuelle
        Psi(cpt) = sum(sum(G_F.^2))/S;
        % Affichage
        disp(['   iteration ',num2str(cpt),' Psi=',num2str(Psi(cpt))]);

            
        % TESTS -----------------------------------------------------------
        
        % Test oscillations
        if cpt>=2
            if Psi(cpt)>Psi(cpt-1);% | % Raffiner
                %if norm(dli + dli_0)<=10*epsilon(etape) % Psi(cpt)>Psi(cpt-1)% | % Raffiner
                warning('   Solution oscillante');
                % Bricolage : diminution du pas
%                li = li - dli/10;
                % Bricolage : on sort avec les valeurs precedentes
                li = li - dli;
                Psi(cpt) =  Psi(cpt-1); 
                dli = 0;
                % Champ actuel
                [x1,x2] = feval(nom_champ,'champ',li,X1,X2);
                % Image G = image g interpolee
                G = interp2(g,x1,x2,type_interp,NaN);
                % Son domaine de definition intersection avec celui de F
                DG = isfinite(G) & DF;
                % Surface = nombre de pixels correles
                S = sum(sum(DG));
                % Ecart GG-F
                G_F = G-F;
                % On annule l'effet des points hors domaine
                G_F(~DG) = 0;
                % Erreur actuelle
                Psi(cpt) = sum(sum(G_F.^2))/S;
             end
        end
        % Test nombre d'iteration
        if cpt >=cptmax;warning ('    Nombre max d''iteration atteint');end

    end

    
    % DISPLAYS --------------------------------------------------------

    
%     [p1,p2] = [x1,x2] = feval(nom_champ,'champ',li,P1,P2);
%     PlotImage(g,2,p1,p2);
%     PlotEcard(G_F,3,P1,P2);
    figure(4); 
    semilogy(Psi','o-','MarkerSize',4,'LineWidth',2,'MarkerFaceColor','w');
    set(gca,'FontSize',myfontsize,'CLim',[0 0.25]);
    xlabel iteration;ylabel \psi;

    
end


% MAILLAGE EF ET INIT CHAMP DEPLACEMENT -----------------------------------


disp('======= Calcul avec champ elements finis Q4');
disp(['   Info : taille de l''image ',num2str(size(F,1)),' x ',num2str(size(F,2)),' pixels']);
TailleEF = input('   Taille des elements finis : ');

% Chargement images brutes (securite)
load Temporaire/F0;
load Temporaire/g0;

% Maillage et domaines de definition des noeuds DN et de l'image DF (reduit)
[XN1,XN2,DN,DF] = Mailleur_EF(F,TailleEF,P1,P2);
XNP1 = XN1;XNP1(~DN) = NaN;
XNP2 = XN2;XNP2(~DN) = NaN;

% Projection de champ a la barbare
disp('   Initialisation sur le champ precemment obtenu');
UN1 = interp2(X1,X2,x1-X1,XN1,XN2);
UN2 = interp2(X1,X2,x2-X2,XN1,XN2);

% Init des parametres lambda_i (li)
li = 0;
l  = 0;
for p=1:size(XN1,1)
    for q=1:size(XN2,2)
        if DN(p,q);
            l = l+2;
            li(l-1) = UN1(p,q);
            li(l  ) = UN2(p,q);
        end
    end
end


% DISPLAYS IMAGES EN COURS DE TRAITEMENT ----------------------------------


PlotImage(F,1,P1,P2,XNP1,XNP2);
PlotImage(g,2);


% INIT NEWTON EF ----------------------------------------------------------


% Champ actuel
x1 = X1;
x2 = X2;
% Init position des noeuds dans l'image g
xn1 = XN1;
xn2 = XN2;
% Calcul de la nouvelle position des pixels depuis le champ EF
l = 0;
disp('   Premier calcul du mapping');
for p=1:size(XN1,1) 
    for q=1:size(XN2,2)
        if DN(p,q)
            % Les noeuds hors domaine ne bougent pas
            l = l+2;
            % Carte aux noeuds
            xn1(p,q) = XN1(p,q)+li(l-1);
            xn2(p,q) = XN2(p,q)+li(l  );
            % Fonction de base centree au noeud
            LB1 = max(0,1 - abs(LX1-XN1(1,q))/TailleEF);
            LB2 = max(0,1 - abs(LX2-XN2(p,1))/TailleEF);
            FB   = LB2'*LB1;
            % Carte - les points hors domaine se deplacement de 0
            x1 = x1 + li(l-1)*FB;
            x2 = x2 + li(l  )*FB;
        end
    end
end
% Image G = image g interpolee
G = interp2(g,x1,x2,type_interp,NaN);
% Son domaine de definition intersection avec celui de F
DG = isfinite(G) & DF;
% Surface = nombre de pixels correles
S = sum(sum(DG));
% Ecart GG-F
G_F = G-F;
% On annule l'effet des points hors domaine
G_F(~DG) = 0;
% Erreur actuelle
Psi = sum(sum(G_F.^2))/S;

disp('   Calcul de la matrice'); % ne change pas en fonction des li
l = 0;  % indice de li
M = zeros(max(size(li)));
% Gradients
[gradF10,gradF20]  = gradient(F);
GF1GF1 = gradF10.*gradF10;
GF1GF2 = gradF10.*gradF20;
GF2GF2 = gradF20.*gradF20;

for p=1:size(XN1,1)
    disp(['   ',num2str(round(100*p/size(XN1,1))),'%']);
    for q=1:size(XN2,2)
        if DN(p,q)
            % Vecteur
            l = l+2;
            % Fonction de base centree au noeud
            LB1 = max(0,1 - abs(LX1-XN1(1,q))/TailleEF);
            LB2 = max(0,1 - abs(LX2-XN2(p,1))/TailleEF);
            m = 0;  % second indice de li
            for r=1:size(XN1,1)
                for s=1:size(XN2,2)
                    if DN(r,s)
                        m = m+2;
                        if m<=l % Symetrie OK
                            if abs(p-r)<=1 & abs(q-s)<=1 % Seult. noeuds voisins  
                                % Fonction de base centree au noeud
                                LC1 = max(0,1 - abs(LX1-XN1(1,s))/TailleEF);
                                LC2 = max(0,1 - abs(LX2-XN2(r,1))/TailleEF);
                                M(l-1,m-1) = (LB2.*LC2)*GF1GF1*(LB1.*LC1)';
                                M(l-1,m  ) = (LB2.*LC2)*GF1GF2*(LB1.*LC1)';
                                M(l  ,m  ) = (LB2.*LC2)*GF2GF2*(LB1.*LC1)';
                                M(l  ,m-1) = M(l-1,m);
                                M(m-1,l-1) = M(l-1,m-1);
                                M(m  ,l-1) = M(l-1,m  );
                                M(m-1,l  ) = M(l  ,m-1);
                                M(m  ,l  ) = M(l  ,m  );
                            end
                        end
                    end
                end
            end
        end
    end
end

% Verification de la determination de la matrice 
% tres lent sur gros maillage
% detM = abs(det(M));
% if detM<1E-6
%     warning(['   La matrice a determinant ',num2str(detM),' presque nul)']);
% end


% NEWTON EF ---------------------------------------------------------------


% Inits
dli = inf*ones(size(li));
% Init compteur de boucles et ecart
cpt = 0;

disp('   Newton');
while max(abs(dli))>precision && cpt<cptmax
    
    cpt = cpt+1;
    % Elimination des points hors domaine
    gradF1 = gradF10;
    gradF2 = gradF20;
    gradF1(~DG) = 0;
    gradF2(~DG) = 0;
    % Calcul vecteur
    V = zeros(max(size(li)),1);
    l = 0;
    G_FG1 = G_F.*gradF1;
    G_FG2 = G_F.*gradF2;
    for p=1:size(XN1,1)
        for q=1:size(XN2,2)
            if DN(p,q)
                % Vecteur
                l = l+2;
                % Fonction de base centree au noeud
                LB1 = max(0,1 - abs(LX1-XN1(1,q))/TailleEF);
                LB2 = max(0,1 - abs(LX2-XN2(p,1))/TailleEF);
                V(l-1) = -LB2*G_FG1*LB1';
                V(l  ) = -LB2*G_FG2*LB1';
            end
        end
    end
    
    % Correcteurs
    dli = (M\V)';
    li = li+dli;
    
    % Coo des noeuds dans repere de g
    xnp1 = xn1;xnp1(~DN)=NaN;
    xnp2 = xn2;xnp2(~DN)=NaN; 
    
    % Champ actuel
    x1 = X1;
    x2 = X2;
    % Init position des noeuds dans l'image g
    xn1 = XN1;
    xn2 = XN2;
    % Calcul de la nouvelle position des pixels depuis le champ EF
    l = 0;
    for p=1:size(XN1,1) 
        for q=1:size(XN2,2)
            if DN(p,q)
                % Les noeuds hors domaine ne bougent pas
                l = l+2;
                % Carte aux noeuds
                xn1(p,q) = XN1(p,q)+li(l-1);
                xn2(p,q) = XN2(p,q)+li(l  );
                % Fonction de base centree au noeud
                LB1 = max(0,1 - abs(LX1-XN1(1,q))/TailleEF);
                LB2 = max(0,1 - abs(LX2-XN2(p,1))/TailleEF);
                FB   = LB2'*LB1;
                % Carte - les points hors domaine se deplacement de 0
                x1 = x1 + li(l-1)*FB;
                x2 = x2 + li(l  )*FB;
            end
        end
    end
    % Image G = image g interpolee
    G = interp2(g,x1,x2,type_interp,NaN);
    % Son domaine de definition intersection avec celui de F
    DG = isfinite(G) & DF;
    % Surface = nombre de pixels correles
    S = sum(sum(DG));
    % Ecart GG-F
    G_F = G-F;
    % On annule l'effet des points hors domaine
    G_F(~DG) = 0;
    % Erreur actuelle
    Psi(cpt) = sum(sum(G_F.^2))/S;
    % Informations
    disp(['   iteration ',num2str(cpt),' Psi=',num2str(Psi(cpt)),' max(|dli|)=',num2str(max(abs(dli)))]);
    % rep = input('continuer []');if ~isempty(rep);error;end
    

    % TESTS -----------------------------------------------------------


    % Test oscillations
    if cpt>=2
        if Psi(cpt)>Psi(cpt-1);% | % Raffiner
            %if norm(dli + dli_0)<=10*epsilon(etape) % Psi(cpt)>Psi(cpt-1)% | % Raffiner
            warning('   Solution oscillante');
            % Bricolage : on sort avec les valeurs precedentes
            %li = li - dli;
            %Psi(cpt) =  Psi(cpt-1); 
            %dli = 0;
        end
    end
    % Test nombre d'iteration
    if cpt >=cptmax;warning ('    Nombre max d''iteration atteint');end

end


% SORTIES GRAPHIQUES ------------------------------------------------------


figure(4); 
semilogy(Psi','o-','MarkerSize',4,'LineWidth',2,'MarkerFaceColor','w');
set(gca,'FontSize',myfontsize,'CLim',[0 0.25]);
xlabel iteration;ylabel \psi;

PlotImage(g,2,xnp1,xnp2);
PlotEcard(G_F,3,XNP1,XNP2);
DEFORMEE(F,DG,X1,X2,x1,x2,XNP1,XNP2,xnp1,xnp2,[],[],[],[]);
CARTES(F,DG,X1,X2,x1,x2,XN1,XN2);
PICKVAL(F,DG,X1,X2,x1,x2,XNP1,XNP2);


% SORTIES SPECIFIQUE DIC CHALLENGE ----------------------------------------
% 
% U2 = x2-X2;
% figure(10);plot(U2(251,:))
% 
% figure(11)
% pcolor(X1,X2,U2)
% shading interp
% colormap hot
% colorbar
% caxis([-0.5 0.5])
% axis equal tight off
% 
