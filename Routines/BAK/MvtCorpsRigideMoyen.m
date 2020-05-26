function [u1cs,u2cs,theta_cs] = MvtCorpsRigideMoyen(X1,X2,U1,U2,DG);
% Parametres du mouvement de corps rigide le plus proche
% translations u0, v0
% rotation theta
% d'un champ de deplacement quelconque U,V
% attention les points de mesure doivent être régulièrement répartis sur
% une grille
% Calculs sur doc jointe
% Pas utilise car necessite un champplein, sans trous. Voir un jour...

% ===================== TESTS ===================== 
% 
% clear
% x = [1:12];
% y = [1:10];
% [X1,X2] = meshgrid(x,y);
% 
% % Champ mesure --> entree du programme ---> a retrouver
% u1cs = -1;
% u2cs = 12;
% theta_cs = -0.25;
% U1 = u1cs - theta_cs*X2;
% U2 = u2cs + theta_cs*X1;
% DG = logical(ones(size(X1)));
% 
% ===================== CALCULS ===================== 
% Calcul du mouvement de coprs rigide le plus proche 
% voir bibliographie "Calculs specifiques"
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
disp('    Mouvement de corps rigide moyen :'); 
disp(['    translation (',num2str(u1cs),',',...
num2str(u2cs),')']);
disp(['    angle ',num2str(theta_cs),' radians soit ',...
num2str(theta_cs*180/pi),' degres']);
% Pour calculer le champ de deplacement
% U1CS = u1cs-theta_cs*X2;
% U2CS = u2cs+theta_cs*X1;
