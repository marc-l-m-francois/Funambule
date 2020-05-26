% function GenImages_DIC;
clear all
path(path,'../Champs');
path(path,'../Routines_graphiques');


% =========================================================================
%
%                             INITS
%
% =========================================================================


% REGLAGES DU CHAMP -------------------------------------------------------

% Nom du champ a considerer
NOMCHAMP = 'EBHOMHPP';
% Les parametres du champ
li = [1 2 3 0.1 0.3 4 5];
% Taille de l'image a generer
n1 = 200; % lignes
n2 = 200; % colonnes
% Pour le champ
global sizeF
sizeF = [n1,n2];


% REGLAGES DE CE PROGRAMME ------------------------------------------------

% Taille moyenne des taches (pixels) 2*R0 = 3 .. 5 pixels optimal
R0 = 2;
% Raffinement de la sous-image
raf = 10;
% Amplitude du bruit sur la position des sous pixels
%nbr = 5;
%br = 1/(nbr*2*raf);


% =========================================================================
%
%                             CORPS
%
% =========================================================================


% INITS DES IMAGES --------------------------------------------------------

% Images initialement blacnhes
F = ones(n1,n2);
G = ones(n1,n2);

% INIT DES IMAGES FINES ---------------------------------------------------

% Coordonnees des points
xf1 = [0.5+1/(2*raf):1/raf:n2+0.5-1/(2*raf)];
xf2 = [0.5+1/(2*raf):1/raf:n1+0.5-1/(2*raf)];
[XF1,XF2] = meshgrid(xf1,xf2);
% Rajout d'un petit bruit pour eviter les soucis de shift de 1/2 image
%XF1 = XF1 + (rand(size(XF1))-0.5)*br;
%XF2 = XF2 + (rand(size(XF2))-0.5)*br;
% Image fine binaire
FF = logical(ones(size(XF1)));
% Mapping : coo image fine deformee
[xf1,xf2] = feval(NOMCHAMP,'champ',li,XF1,XF2);

% GENERATION IMAGES FINES -------------------------------------------------

disp('> Generation image fine');
% Niveau de gris moyen de l'image
F_moy = 1;
while F_moy>0.5
    % Rayon tache noire
    R = R0*rand;
    % Position du centre entre 0.5-R et n+0.5+R (tangente au bords)
    C1 = 0.5-R + (n2+2*R)*rand;
    C2 = 0.5-R + (n1+2*R)*rand;
    % Sous-pixel noir si dans tache
    FF( (XF1-C1).^2 + (XF2-C2).^2 <= R^2 ) = false;
    % Niveau de gris actuel image fine
    F_moy = sum(sum(FF))/prod(size(FF)); % (n1*n2*raf^2)
    disp(F_moy);
end

% GENERATION IMAGES INITIALE F ET FINALE G --------------------------------

disp('> Calcul images F et G');
for p=1:n1
    for q=1:n2
        % Image EI
        F(p,q) = sum(sum(FF(XF1>=q-0.5 & XF1< q+0.5 & XF2>=p-0.5 & XF2< p+0.5)))/raf^2;
        % Image EF
        % nb de sous pixels pas forcement = raf^2
        inds = xf1>=q-0.5 & xf1< q+0.5 & xf2>=p-0.5 & xf2< p+0.5;
        nraf = sum(sum(inds));
        if nraf>0
            G(p,q) = sum(sum(FF(inds)))/nraf;
        end
    end
    Display_avancement(p,n1);
end
% Verification niveau de gris possible disp([mean(mean(F)),F_moy]);


% =========================================================================
%
%                         REPRESENTATION DU RESULTAT
%
% =========================================================================


PlotImage(F,1);set(1,'Name','Image F');
PlotImage(G,2);set(2,'Name','Image G');


% =========================================================================
%
%                       ENREGISTREMENT TIFF 16 bits
%
% =========================================================================

extension = [NOMCHAMP,'_'];
for p=1:size(li,2)
    extension = [extension,num2str(li(p)),'_'];
end
imwrite(uint16(F*65280),[extension,'F.tiff'])
imwrite(uint16(G*65280),[extension,'G.tiff'])


