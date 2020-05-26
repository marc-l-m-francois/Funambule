function F = FlouGaussien(F,rg);

if rg>0
	%siz = round(50/10);
	%H = ones(siz,siz)/siz^2;
	% H = fspecial('gaussian',6*ceil(rg)+1,rg);
	
	H = fspecial('gaussian',6*ceil(rg)+1,rg/2);
	
	
	%warning off;figure(5);
	%set(5,'Name','Filtre gaussien utilise');
	%warning on
	%surf(H);shading faceted;
	F = filter2(H,F);
end
