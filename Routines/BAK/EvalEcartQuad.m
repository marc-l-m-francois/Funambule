function Delta = EvalEcartQuad(L,L0,courbe);
% Pour le papier
% Ecart quadratique entre courbe et courbe de reference
% Si tu as un L0 re reference bien sur

raf = 1000;
x1 = [0:1/raf:1];

[X01 ,X02] = feval(courbe,L0,x1);
[Xc1 ,Xc2] = feval(courbe,L ,x1);

d2 = zeros(size(x1));
for p=1:size(x1,2)
   d2(p) = min( (Xc1-X01(p)).^2 + (Xc2-X02(p)).^2 );
end

% Ecart RMS - because ça marche super et que la moyenne n'est pas nulle
Delta = rms(sqrt(d2));