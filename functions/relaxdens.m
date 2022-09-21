function [r, drdx] = relaxdens(x,p,q,rhomin, type)

% [r, drdx] = relaxdens(x,p,q)
%     Compute relaxed density for stress definition
%
%     The argument x must be the raw (design variable) density.
%
%     'type' indicates which type of relaxation function to use, and it can
%     be
%       'stdpq'   : standard pq-relaxation based on raw density
% 	  	'modpq'   : modified pq-relaxation based on modified SIMP densities
%

switch type
    case 'stdpq'
        r = x.^(p-q);
        drdx = (p-q)*x.^(p-q-1);
    case 'modpq'
        rp = rhomin + (1-rhomin)*x.^p;
        rq = rhomin + (1-rhomin)*x.^q;
        r = rp./rq;
        drdx = (1-rhomin)*(p*x.^(p-1) - q*(x.^(q-1)).*r)./rq;
    otherwise
        disp 'Error in relaxdens: unknown relaxation type.';
end

end