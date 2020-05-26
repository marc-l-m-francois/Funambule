function [IM,prof] = PrepImage(IMC);

% Convertit une image IMC en RVB (ou pas) en image noir et blanc IM

% Cas de 4eme couleur transparente (tif)
if size(IMC,3)>3;
    IMC=IMC(:,:,1:3);
    % disp('    Elimination du 4e canal transparent');
end;
% Traitement premilinaire couleur ou nb
if size(IMC,3)==3;
    coul=true;strcoul = 'color';
else;
    coul=false;strcoul = 'black and white';
end
% Traitement preliminaire profondeur de couleurs
if isa(IMC,'logical');
    prof = 1;
elseif isa(IMC,'uint8' );
    prof = 8;
elseif isa(IMC,'uint16');
    prof = 16;
end

if coul
% Ponderation par defaut Norme NTSC 0.2989*rouge + 0.5870*vert + 0.1140*bleu
    rvb = [0.2989 0.5870 0.1140];
    % Premiere conversion
	IM = 	rvb(1)*double(IMC(:,:,1))/(sum(rvb)*(2^prof-1)) + ...
         	rvb(2)*double(IMC(:,:,2))/(sum(rvb)*(2^prof-1)) + ...
         	rvb(3)*double(IMC(:,:,3))/(sum(rvb)*(2^prof-1));
else
	% Images en niveaux de gris
	IM = double(IMC)/(2^prof-1);
end

% Infos a l'utilisateur
disp(['    ',num2str(size(IM,1)),'x',num2str(size(IM,2)),' pixels, ',...
    num2str(prof),' bits, initially in ',strcoul]);
