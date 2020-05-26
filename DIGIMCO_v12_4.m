disp(' ');
disp('       -----------------------------------------------------');
disp('       |2008                                           2019|')
disp('       |                                                   |')
disp('       |                      DigImCo                      |')
disp('       |            DIGital IMage COrrelation              |');
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
verif_grad = true;
% Nombre maxi d'iterations
cptmax     = 200;
% Limite de variation des li
precision  = 1E-3;
% Limite de variation du correcteur des parametres li
limdli = 5;


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
% Taille de F en global car utilise dans certain champs
global sizeF
sizeF = size(F);


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
set(5,'MenuBar','none','GraphicsSmoothing','off','ToolBar','none');

% Deformee sans corps rigide
figure(6);Range_mon_trace(6,2,3,2,4);clf;
set(6, 'MenuBar','none','GraphicsSmoothing','off','ToolBar','none');


% ETAPES INITIALES AVEC FLOU -----------------------------------------------


if ToolBoxImages
    disp('    Un calcul initial sur une image floutee aide la convergence en cas de grands deplacements');
    rayon_de_flou = input('    Entrer le ou les niveaux de flou souhaite par, ex. [12 6 0], ou [0] si non : ');
else
    rayon_de_flou = 0;
end
if rayon_de_flou(end)~=0;rayon_de_flou=[rayon_de_flou,0];end


% DOMAINE DE CORRELATION --------------------------------------------------


[P1,P2,DF] = ChoixDomaineCorrelation(F,X1,X2);


% CHOIX DU CHAMP ----------------------------------------------------------


disp('======= Choix du champ =======');
cd Champs/;a = dir;cpt=0;
for p=1:size(a,1)
    if strcmp(a(p).name(end),'m') && ~a(p).isdir && size(a(p).name,2)>2 && ~strcmp(a(p).name,'Q4.m')
        if strcmp(a(p).name(end-1),'.')
            cpt = cpt+1;disp([num2str(cpt),' ',a(p).name(1:end-2)]);num(cpt) = p;
        end
    end
end
rep = input('Le numero de votre choix ? ');
nom_champ = a(num(rep)).name(1:end-2);cd ..


% INIT DES LI -------------------------------------------------------------


li = zeros(1,feval(nom_champ,'nombre_param'));


% VERIFICATION DES DERIVEES ------------------------


if verif_grad
    li = (rand(size(li))-0.5)/100;
    [x1,x2] = feval(nom_champ,'champ',li,X1,X2);
    for p=1:size(li,2)
        % Derivee analytique
        [dx1dL,dx2dL] = feval(nom_champ,'derivees',li,X1,X2,p);
        % Derivee numerique
        dli = zeros(size(li)); dli(p) = precision;
        [x1p,x2p] = feval(nom_champ,'champ',li+dli,X1,X2);
        dx1dLnum = (x1p-x1)/precision;  dx2dLnum = (x2p-x2)/precision;
        % Comparaison
        if norm([dx1dLnum;dx2dLnum]-[dx1dL;dx2dL],'fro')/...
           numel([dx1dLnum;dx2dLnum])>precision
           disp(['   derivee numero ',num2str(p)]);
           warning(' fausse');
           disp(['   norme de l''ecart ',num2str(norm([dx1dLnum;dx2dLnum]-[dx1dL;dx2dL],'fro')/...
           numel([dx1dLnum;dx2dLnum])),' pour max autorise ',num2str(precision)]);
%            % activer pour debug
%            disp('dx1dL , dx1dLnum');
%            disp(dx1dL);
%            disp(dx1dLnum);
%            disp('dx2dL , dx2dLnum');
%            disp(dx2dL);
%            disp(dx2dLnum);
            error
        else
            disp(['   derivee numero ',num2str(p),' verifiee']);
        end
    end
end


% CALCUL DE CORRELATION AVEC PUIS SANS FLOU -------------------------------


disptitle('Calcul de correlation');

% Init des parametres
% li = zeros(size(li));
li = (rand(size(li))-0.5)/100; % Certains champ (EBHOMHPP) n'apprecient pas l'init a 0...

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

    
	% DISPLAY IMAGES EN COURS DE TRAITEMENT -------------------------------


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
    disp(['   iteration ',num2str(1),' Psi=',num2str(Psi(cpt))]);
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

        % Verification de la determination de la matrice (peu utile)
        % detM = abs(det(M));
        % if detM<1E-6
        % 	warning(['      La matrice a determinant ',num2str(detM),' presque nul)']);
        % end       

        % Correcteur
        dli = (M\V)';
        % Possible limitation en cas de probleme partiellement mal defini
        if all(~isfinite(dli))
            error('Systeme entierement mal conditionne');
        elseif any(~isfinite(dli))
            warning('Systeme partiellement mal conditionne');
            dli(~isfinite(dli)) = 0;
        end
        % Possible limitation du correcteur ICI
        if max(abs(dli))>limdli
            disp('   Limitation du correcteur');
            dli = dli*limdli/max(abs(dli));
        end
        li = li + dli;
        
        disp(li);

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
%                 li = li - dli;
%                 Psi(cpt) =  Psi(cpt-1); 
%                 dli = 0;
%                 % Champ actuel
%                 [x1,x2] = feval(nom_champ,'champ',li,X1,X2);
%                 % Image G = image g interpolee
%                 G = interp2(g,x1,x2,type_interp,NaN);
%                 % Son domaine de definition intersection avec celui de F
%                 DG = isfinite(G) & DF;
%                 % Surface = nombre de pixels correles
%                 S = sum(sum(DG));
%                 % Ecart GG-F
%                 G_F = G-F;
%                 % On annule l'effet des points hors domaine
%                 G_F(~DG) = 0;
%                 % Erreur actuelle
%                 Psi(cpt) = sum(sum(G_F.^2))/S;
%                 % Affichage
%                 disp(['   iteration ',num2str(cpt),' Psi=',num2str(Psi(cpt))]);
            end
        end
        % Test nombre d'iteration
        if cpt >=cptmax;warning ('    Nombre max d''iteration atteint');end
        
        
        % SORTIES GRAPHIQUES EN COURS DE NEWTON ---------------------------

        
    end

    
    % SORTIES GRAPHIQUES EN FIN DE NIV DE FLOU ----------------------------
    
    
    % Fig 2 Maille deforme sur g
    [p1,p2] = feval(nom_champ,'champ',li,P1,P2);
    PlotImage(g,2,p1,p2);
    % Fig 3 Ecart
    PlotEcard(G_F,3,P1,P2);
    % Fig 4 Convergence
    figure(4);semilogy(Psi','o-','MarkerSize',4,'LineWidth',2,'MarkerFaceColor','w');
    set(gca,'FontSize',myfontsize,'CLim',[0 0.25]); xlabel iteration;ylabel \psi;

    
end
        

% SORTIES TEXTE EN FIN DE CALCUL ------------------------------------------


disptitle('resultats')
[Xx1,Xx2] = feval(nom_champ,'infos',li);
% Lignes caracteristique
if ~isempty(Xx1)
    XC1=Xx1(1,:);XC2=Xx2(1,:);xc1 = Xx1(2,:);xc2=Xx2(2,:);
else
    XC1 = [];XC2=[];xc1=[];xc2=[];
end
figure(1),plot(XC1,XC2,'r-');figure(2),plot(xc1,xc2,'r-');


% SORTIES GRAPHIQUES EN FIN DE CALCUL -------------------------------------


disptitle('sorties graphiques')

% Fig 5 et 6 deformee amplifiee
DEFORMEE(F,DG,X1,X2,x1,x2,P1,P2,p1,p2,XC1,XC2,xc1,xc2);
% Fig 7 cartes
CARTES(F,DG,X1,X2,x1,x2,P1,P2);
% Fig 1 picking des valeurs des champs
PICKVAL(F,DG,X1,X2,x1,x2,P1,P2);

