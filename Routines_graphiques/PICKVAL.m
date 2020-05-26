function PICKVAL(F,DG,X1,X2,x1,x2,P1,P2);
% Sort les valeurs en un point
% X1 X2 coo des pixels dans repere de F
% x1 x2 coo des pixels dans repere de g
% DG domaine de def correlation
% P1 P2 coo des points caracteristiques (bords ou maillage) rep F
% p1 p2 coo des points caracteristiques (bords ou maillage) rep g

% Elimination des points courants hors domaine
x1(~DG) = NaN;
x2(~DG) = NaN;

% Calcul des deplacements
U1 = x1-X1;
U2 = x2-X2;

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

% Tenseur des petites deformations H=F-I
% Almansi-Euler A = (eye(2)-inv(FF*FF'))/2;
% Rotations
% FF = [F11 F12;F21 F22];
% R = FF*(FF'*FF)^(-0.5);
% R11 = R(1,1);R12=R(1,2):R21=R(2,1)
% Angle de rotation
% theta = angle(R(1,1)+i*R(1,2));

disp('======= Deplacement et deformation en un point =======');

PlotImage(F,1,P1,P2);
button = 1;
while button==1
    disp('Cliquer un point. Entree pour sortir.');
    [M1,M2,button] = ginput(1);
    % Sortie aussi en cas d'entree vide
    if isempty(M1)
        button=[];
    else
        PlotImage(F,1,P1,P2);
        plot(M1,M2,'r+','MarkerSize',12);
        IND = X1==round(M1) & X2==round(M2);
        disp(['    Au point (',num2str(X1(IND)),',',num2str(X2(IND)),')']);
        disp(' ');
        disp('    Vecteur deplacement [U1 U2] (pixel) ');
        disp([U1(IND),U2(IND)]);
        disp('    Tenseur des deformations [E11 E12 ; E12 E22] de Green-Lagrange :');
        disp([E11(IND) E12(IND);E12(IND) E22(IND)]);        
        disp('    -----------------------------------');
    end
end
