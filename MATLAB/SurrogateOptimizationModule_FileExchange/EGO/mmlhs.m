function X_best = mmlhs(X_start, population, iterations, q)
% Evolutionary operation search for the most space filling Latin hypercube
% of a certain size and dimensionality. There is no need to call this
% directly - use bestlh.m
%
% Copyright 2007 A Sobester
%

n = size(X_start,1);

X_best = X_start; Phi_best = mmphi(X_best);

leveloff = floor(0.85*iterations);

for it = 1:iterations
	if it < leveloff
		mutations = round(1+(0.5*n-1)*(leveloff-it)/(leveloff-1));
	else
		mutations = 1;
	end
	X_improved  = X_best; Phi_improved = Phi_best;
	
	for offspring = 1:population
		X_try = perturb(X_best, mutations);
		Phi_try = mmphi(X_try, q);
		
		if Phi_try < Phi_improved
			X_improved = X_try;
			Phi_improved = Phi_try;
		end
	end
	
	if Phi_improved < Phi_best
		X_best = X_improved;
		Phi_best = Phi_improved;
	end
	%plot(X_best(:,1),X_best(:,2),'o');drawnow;
end