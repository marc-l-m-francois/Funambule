function D = DetectoDedans(X1,X2,P1,P2);
% fait la meme chose que la fonction inpolygon de Matlab !!
% Pour savoir si les points M(X1,X2) sont au dedans ou au dehors 
% du polygone defini par les points P(P1,P2) (sens indifferend)

% Reglage (peu sensible)
tol = 1E-3;

Q1 = circshift(P1,-1,2);
Q2 = circshift(P2,-1,2);

% Angles thetai entre les vecteurs PMi et PMi+1
alfa = zeros(size(X1));
for p=1:max(size(P1))
    % Points coincidents comptes dedans -1E-12 pour eviter le /0
    alfa = alfa + angle( (Q1(p)-X1+i*(Q2(p)-X2))./(P1(p)-X1+i*(P2(p)-X2)-1E-12) );
end
D = abs(alfa)>tol;