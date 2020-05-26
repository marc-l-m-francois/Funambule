function [x1,x2] = CORSOL(choix,li,X1,X2,p);
% -------------------------------------------------------------------------
% Mouvement de corps solide
% li(1)	Translation suivant x
% li(2)	Translation suivant y
% li(3)	Angle theta de la rotation
% li(4) = e Echelle 1+e
% -------------------------------------------------------------------------
% x1 x2	position image finale
% X1 X2	position image initiale (entiers) global
% li 	parametres du champ : On doit avoir li(n)=0 ==> x1=X1 et x2=X2
% p     indice  du lambda pour derivation
% -------------------------------------------------------------------------


switch(choix)
    
    
    case 'nombre_param'
        
        
        x1 = 4;
        x2 = [];
        
        
    case 'infos'
        

        infos = {...
             ['Translation tx = ',num2str(li(1)),' pixels'];
             ['Translation ty = ',num2str(li(2)),' pixels'];
             ['Rotation    rz = ',num2str(li(3)),' radians = ', num2str(180*li(3)/pi),' degres '];
             ['Echelle    1+e = ',num2str(1+li(4))]
             };
        disp(infos);
        x1 = [];
        x2 = [];
        
        
    case 'champ'
        
        
        x1 = (   X1*cos(li(3)) - X2*sin(li(3)))*(1+li(4)) + li(1);
        x2 = (   X1*sin(li(3)) + X2*cos(li(3)))*(1+li(4)) + li(2);
        
        
    case 'derivees'


        if p==1
            x1 =  ones (size(X1));
            x2 =  zeros(size(X1));
        elseif p==2
            x1 = zeros(size(X1));
            x2 = ones (size(X1)); 
        elseif p==3
            x1 = ( - X1*sin(li(3)) - X2*cos(li(3)))*(1+li(4));
            x2 = (   X1*cos(li(3)) - X2*sin(li(3)))*(1+li(4));
        elseif p==4
            x1 =   X1*cos(li(3)) - X2*sin(li(3));
            x2 =   X1*sin(li(3)) + X2*cos(li(3));
        end

end