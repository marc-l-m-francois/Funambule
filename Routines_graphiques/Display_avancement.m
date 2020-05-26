function Display_avancement(n,N);
% Trace une barre d'avancement d'un process a N etapes pour l''etape n

% Nombre de poins dans la barre d''avancement
nb = 72;
nbo = floor(nb*n/N);

if n==1
	if N==1
		% Barre entiere d''un coup
		fprintf('%c',46*ones(1,nb));
	end
elseif nbo == floor(nb*(n-1)/N)+1
	% Un point
	fprintf('%c',46);
end

if n==N
	% Retour chariot
	fprintf('\n');
end
