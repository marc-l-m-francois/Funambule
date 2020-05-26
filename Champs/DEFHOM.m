function [x1,x2] = DEFHOM(choix,li,X1,X2,p);
% -------------------------------------------------------------------------
% Deformation homogene
% li(1) Translation suivant x
% li(2) Translation suivant y
% [li(3) li(4)
%  li(5) li(6)] composantes du tenseur F-I
% -------------------------------------------------------------------------
% x1 x2	position image finale
% X1 X2	position image initiale (entiers) global
% li 	parametres du champ : On doit avoir li(n)=0 ==> x1=X1 et x2=X2
% p     indice  du lambda pour derivation
% -------------------------------------------------------------------------


switch(choix)
    
    
    case 'nombre_param'
        
        
        x1 = 6;
        x2 = [];
        
        
    case 'infos'
        

        % Tenseur gradient FF (pour ne pas confondre avec image F)
        FF = [1+li(3) li(4);li(5) 1+li(6)];
        % Tenseur des petites deformations
        H = FF-eye(2);
        EPS = (H+H')/2;
        % Dilatations ou Cauchy-Green droit
        C = FF'*FF;
        % Cauchy-Green gauche
        B = FF*FF';
        % Tenseur des deformations de Green-Lagrange
    %     E11 = (F11^2+F21^2-1)/2
    %     E12 = (F11*F12+F21*F22)/2
    %     E22 = (F12^2+F22^2-1)/2
        E = (FF'*FF - eye(2))/2;
        [EV,ED] = eig(E);
        % Almansi-Euler
        A = (eye(2)-inv(FF*FF'))/2;
        % Rotations
        R = FF*(FF'*FF)^(-0.5);
        theta = angle(R(1,1)+i*R(1,2));
        % Displays
        infos = {...
%         En repere initial --> illisible
%         ['Translation tx = ',num2str(li(1)),' pixels']; 
%         ['Translation ty = ',num2str(li(2)),' pixels']; 
%         ['Tenseur H = F-I x1 1000 (contient li(3..6))'];
%         num2str(1000*H(1,1:2));
%         num2str(1000*H(2,1:2));
%         ['Tenseur des petites deformations x1 1000 :'];
%         num2str(1000*EPS(1,1:2));
%         num2str(1000*EPS(2,1:2));
        'Tenseur de Grenn-Lagrange x1 1000 : ';
        num2str(1000*E(1,1:2));
        num2str(1000*E(2,1:2));
%         'Deformations principales de Green Lagrange x1 1000 ';
%         num2str(1000*[ED(1,1),ED(2,2)]);
%         'Rotation (radians, degres) : ';
%         num2str([theta,theta*180/pi]);
        };
        disp(infos);
        x1 = [];
        x2 = [];

        
    case 'champ'
        
        
        x1 = (1+li(3))*(X1+li(1)) + (li(4))*(X2+li(2));
        x2 = (li(5))*(X1+li(1)) + (1+li(6))*(X2+li(2));
        
        
    case 'derivees'


        if p==1
            x1 = (1+li(3))*ones(size(X1));
            x2 = (  li(5))*ones(size(X1));
        elseif p==2
            x1 = (  li(4))*ones(size(X1));
            x2 = (1+li(6))*ones(size(X1));
        elseif p==3
            x1 = (X1+li(1));
            x2 = zeros(size(X2));
        elseif p==4
            x1 = (X2+li(2));
            x2 = zeros(size(X2));
        elseif p==5
            x1 = zeros(size(X1));
            x2 = (X1+li(1));
        elseif p==6
            x1 = zeros(size(X1));
            x2 = (X2+li(2));
        end


end



	
% Hessien	
% 	x1 = zeros(size(X1));
% 	x2 = zeros(size(X2));
% 		
% 	if     (p==1 & q==3) | (p==3 & q==1) | (p==2 & q==4) | (p==4 & q==2)
% 		x1 =  ones(size(X1));
% 	elseif (p==1 & q==5) | (p==5 & q==1) | (p==2 & q==6) | (p==6 & q==2)
% 		x2 =  ones(size(X1));
% 	end
