function [Diag,theta1] = ValeursPropres(T);
% Donne valeurs propres et base propre d'un tenseur du second ordre

[base,Diag] = eig(T);

% Verif
% u = [V(1,1),V(2,1)];
% v = [V(1,2),V(2,2)];

% T - (D(1,1)*u'*u + D(2,2)*v'*v)
theta1 = angle(base(1,1)+i*base(2,1));
theta2 = angle(base(1,2)+i*base(2,2));
% (theta2-theta1)*180/pi