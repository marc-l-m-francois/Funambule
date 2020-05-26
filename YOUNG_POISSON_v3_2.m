% Porgamme MesuresYoungPoisson 
% Mesure Young et Poisson sur une suite d'iages
% Developpe pour les essais sur le chanvre avec Y et M Lecieux des 2017



% ------------------------------ ACCUEIL ----------------------------------


disp(' ');
disp('       -----------------------------------------------------');
disp('       |2017                                           2020|')
disp('       |                                                   |')
disp('       |                   Young Poisson                   |')
disp('       |            DIGital IMage COrrelation              |');
disp('       |        Identification par Images Numeriques       |');
disp('       |                                                   |')
disp('       |Marc Francois                Labo. GeM Univ. Nantes|')
disp('       |marc.francois@univ-nantes.fr                       |');
disp('       -----------------------------------------------------');
disp(' ');

disp('------------------------------- NOTES -------------------------------');
disp('- Les images doivent etre nommees sous le format NN_FXXXX ');
disp('- ou NN est le numero d''ordre et XXXX la force appliquee en N');
disp('- Le chargement est supose vertical');
disp('- Le repertoire ne doit contenir que les images');
disp('- Toutes les images doivent avoir meme nombre de pixels')
disp('---------------------------------------------------------------------');


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


% INITS DIGIMCO ----------------------------------------------------------


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


% INITS SPECIFIQUES YOUNG POISSON -----------------------------------------


% Rayon de flou pour etape initiale defaut 6, (on peut mettre un vecteur)
% Doit finir par 0 pour la precision max
rayon_de_flou = [6 0];
% Type de champ = deformation homogene
nom_champ = 'DEFHOM';


% INITS GRAPHIQUES --------------------------------------------------------


cg = 'r';   % Couleur Couleur bord du domaine
cF = 'b';   % Couleur maillage
myfontsize = 18;
trait_gro = 2;


% CHOIX DE LA LISTE D'IMAGES ----------------------------------------------


ActPath = cd;	% Actuel
%cd Images
string = 'Choisir une image de la serie ';
[nom_F,ImagePath] = uigetfile('*.*',string);
Liste_images0 = dir(ImagePath);
Liste_images = [];
for p=1:size(Liste_images0,1)
    if Liste_images0(p).name(1)~='.';
        Liste_images = [Liste_images;Liste_images0(p).name];
    end
end
clear Liste_images0;


% INITS SORTIES GRAPHIQUES ------------------------------------------------


% Image g
figure(2);Range_mon_trace(2,1,2,2,4);
set(2,'Name','Image etat final g','MenuBar','none',...
    'GraphicsSmoothing','off','ToolBar','none');

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


% DOMAINE DE CORRELATION SUR IMAGE 1 --------------------------------------


% Chargement premiere image
IMAGE = imread(fullfile(ImagePath,Liste_images(1,:)));
F = PrepImage(IMAGE);

% Coordonnees des pixels en repere de F
LX1 = [1:size(F,2)];
LX2 = [1:size(F,1)];
[X1,X2] = meshgrid(LX1,LX2);
% Taille de F en global car utilise dans certain champs
global sizeF
sizeF = size(F);

% Image F
figure(1);clf;Range_mon_trace(1,1,1,2,4);
set(1,'Name','Image etat initial F','MenuBar','none',...
    'GraphicsSmoothing','off','ToolBar','none');
PlotImage(F,1);

% Domaine de correlation
[P1,P2,DF] = ChoixDomaineCorrelation(F,X1,X2);


% CALCUL DE CORRELATION ---------------------------------------------------


% Init des increments de deformation
DE11 = zeros(1,size(Liste_images,1));
DE12 = zeros(1,size(Liste_images,1));
DE22 = zeros(1,size(Liste_images,1));


for N_image=2:size(Liste_images,1)
    disp(['================== Traitement image ',num2str(N_image),'/',...
        num2str(size(Liste_images,1))]);
    % Chargement image de reference
    IMAGE = imread(fullfile(ImagePath,Liste_images(N_image-1,:)));
    F = PrepImage(IMAGE);
    % Chargement image deformee
    IMAGE = imread(fullfile(ImagePath,Liste_images(N_image,:)));
    g = PrepImage(IMAGE);
    % Sauvegarde des images (pour les reprendre apres le flou)
    save Temporaire/F0 F
    save Temporaire/g0 g
    % Init des parametres a 0
    li = zeros(1,feval(nom_champ,'nombre_param'));
    
    
    % FLOUTAGE POUR AIDER LA CONVERGENCE ----------------------------------
    
    
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


        % DISPLAY IMAGES EN COURS DE TRAITEMENT ---------------------------


        PlotImage(F,1,P1,P2);
        PlotImage(g,2);


        % INITS CORRELATION -----------------------------------------------


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


            % SORTIES GRAPHIQUES EN COURS DE NEWTON -----------------------


        end


        % SORTIES GRAPHIQUES EN FIN DE NIV DE FLOU ------------------------


        % Fig 2 Maille deforme sur g
        [p1,p2] = feval(nom_champ,'champ',li,P1,P2);
        PlotImage(g,2,p1,p2);
        % Fig 3 Ecart
        PlotEcard(G_F,3,P1,P2);
        % Fig 4 Convergence
        figure(4);semilogy(Psi','o-','MarkerSize',4,'LineWidth',2,'MarkerFaceColor','w');
        set(gca,'FontSize',myfontsize,'CLim',[0 0.25]); xlabel iteration;ylabel \psi;

    end
    
    
    % TRAITEMENT ET STOCKAGE DES MESURES
    
    % Tenseur gradient (init a identite)
    nbp = 1;
    F11  = 1+li(3)/nbp;
    F12  =   li(4)/nbp;
    F21  =   li(5)/nbp;
    F22  = 1+li(6)/nbp;
    % Tenseur des deformations de Green-Lagrange(incremental)
    DE11(N_image) = (F11^2+F21^2-1)/2;
    DE12(N_image) = (F11*F12+F21*F22)/2;
    DE22(N_image) = (F12^2+F22^2-1)/2;
    
    % Stockage du domaine deforme pour futur calcul
    P1 = p1;
    P2 = p2;
    
    % Stockage erreur de correlation finale
    psis(N_image) = Psi(end);
    
    pause(0.1);
end


% POST-TRAITMENT ----------------------------------------------------------


disp(' ')
disp('====== RESULTATS ========');
disp(' ')


Surface = input('Surface d''appui initiale en m2 ? ');

% Force depuis le nom des images
Force = -str2num(Liste_images(:,5:9))';

% Tenseur des deformation de Green Lagrance
E11 = cumsum(DE11);
E12 = cumsum(DE12);
E22 = cumsum(DE22);

% Tenseur de Piola Kirchhoff 2
for N_image=2:size(Liste_images,1)
    % Contrainte de traction (selon 2)
    SIGL(N_image) = (Force(N_image)-Force(1))/Surface;
    % Young secant
    YoungSec(N_image) =  SIGL(N_image)/E22(N_image);
    % Poisson secant
    PoissonSec(N_image) = -E11(N_image)/E22(N_image);
end
% Young tangent
YoungTan = diff(SIGL)./diff(E22);
PoissonTan = -diff(E11)./diff(E22);


% TRACES ------------------------------------------------------------------

figure(10);Range_mon_trace(10,1,1,2,4);
figure(11);Range_mon_trace(11,1,2,2,4);
figure(12);Range_mon_trace(12,2,1,2,4);
figure(13);Range_mon_trace(13,2,2,2,4);
figure(14);Range_mon_trace(14,1,3,2,4);
figure(15);Range_mon_trace(15,2,3,2,4);
figure(16);Range_mon_trace(16,1,4,2,4);


% Contrainte-deformation
figure(10);clf;hold on;grid on
set(gcf,'MenuBar', 'None','ColorMap',gray);
plot(E22,SIGL/1E6,'ko-','LineWidth',trait_gro,'MarkerSize',6,'MarkerFaceColor',3*[1 1 1]/4);
plot(0,0,'k.','markersize',0.01);
% for p=1:size(E22,2)
%     text(E22(p),SIGL(p)/1E6,num2str(p))
% end
xlabel \epsilon_{L}
ylabel ('\sigma (MPa)')
set(gca,'FontSize',myfontsize,'CLim',[0 1]);
set(10,'Name',['Contrainte - deformation'],'Color',9*[1 1 1]/10,'NumberTitle','off');
print Fig_Contrainte_Deformation -deps2 -f10


% Module d'Young Secant
figure(11);clf;hold on;grid on
set(gcf,'MenuBar', 'None','ColorMap',gray);
plot(E22(2:end),YoungSec(2:end)/1E6,'ko-','LineWidth',trait_gro,'MarkerSize',6,'MarkerFaceColor',3*[1 1 1]/4);
plot(0,0,'k.','markersize',0.01);
xlabel \epsilon_{L}
ylabel ('E (MPa)')
set(gca,'FontSize',myfontsize,'CLim',[0 1]);
set(11,'Name',['Module d''Young secant'],'Color',9*[1 1 1]/10,'NumberTitle','off');
print Fig_Young_secant -deps2 -f11


% Coefficient de Poisson Secant
figure(12);clf;hold on;grid on
set(gcf,'MenuBar', 'None','ColorMap',gray);
plot(E22(2:end),PoissonSec(2:end),'ko-','LineWidth',trait_gro,'MarkerSize',6,'MarkerFaceColor',3*[1 1 1]/4);
plot(0,0,'k.','markersize',0.01);
xlabel \epsilon_{L}
ylabel \nu
set(gca,'FontSize',myfontsize,'CLim',[0 1]);
set(12,'Name',['Coefficient de Poisson secant'],'Color',9*[1 1 1]/10,'NumberTitle','off');
print Fig_Poisson_secant -deps2 -f12


% Trajet de deformation
figure(13);clf;hold on;grid on
set(gcf,'MenuBar', 'None','ColorMap',gray);
plot(E22,E11,'ko-','LineWidth',trait_gro,'MarkerSize',6,'MarkerFaceColor',3*[1 1 1]/4);
plot(0,0,'k.','markersize',0.01);
xlabel \epsilon_{L}
ylabel \epsilon_{T}
set(gca,'FontSize',myfontsize,'CLim',[0 1]);
set(13,'Name',['Deformation L - Deformation T'],'Color',9*[1 1 1]/10,'NumberTitle','off');
print Fig_Deformations_L_T -deps2 -f13


% Erreur de correlation
figure(14);clf;hold on;grid on
set(gcf,'MenuBar', 'None');
plot(E22(2:end),psis(2:end),'ko-','LineWidth',trait_gro,'MarkerSize',6,'MarkerFaceColor',3*[1 1 1]/4);
plot(0,0,'k.','markersize',0.01);
xlabel \epsilon_{L}
ylabel('erreur \Psi')
set(gca,'FontSize',myfontsize);%,'CLim',[0 1]);
set(14,'Name',['Erreur de correlation'],'Color',9*[1 1 1]/10,'NumberTitle','off');
print Fig_Erreur_correlation -deps2 -f14


% Module d'Young Tangent
figure(15);clf;hold on;grid on
set(gcf,'MenuBar', 'None','ColorMap',gray);
plot(E22(2:end),YoungTan/1E6,'ko-','LineWidth',trait_gro,'MarkerSize',6,'MarkerFaceColor',3*[1 1 1]/4);
plot(0,0,'k.','markersize',0.01);
xlabel \epsilon_{L}
ylabel ('E (MPa)')
set(gca,'FontSize',myfontsize,'CLim',[0 1]);
set(15,'Name',['Module d''Young tangent'],'Color',9*[1 1 1]/10,'NumberTitle','off');
print Fig_Young_tangeant -deps2 -f15


% Coefficient de Poisson Tangent
figure(16);clf;hold on;grid on
set(gcf,'MenuBar', 'None','ColorMap',gray);
plot(E22(2:end),PoissonTan,'ko-','LineWidth',trait_gro,'MarkerSize',6,'MarkerFaceColor',3*[1 1 1]/4);
plot(0,0,'k.','markersize',0.01);
xlabel \epsilon_{L}
ylabel \nu
set(gca,'FontSize',myfontsize,'CLim',[0 1]);
set(16,'Name',['Coefficient de Poisson tangent'],'Color',9*[1 1 1]/10,'NumberTitle','off');
print Fig_Poisson_tangeant -deps2 -f16


% Enregistrement des resultats obtenus
% save 'Mesure_Young_Poisson_resultats' E11 E22 E12 SIGL YoungSec PoissonSec YoungTan PoissonTan psis Surface Force







































