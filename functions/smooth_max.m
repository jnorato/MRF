function [S, dSdx] = smooth_max(x,p,form_def)
%
% This function computes a smooth approximation of the maximum of x.  The
% type of smooth approximation (listed below) is given by the argument
% form_def, and the corresponding approximation parameter is given by p.
%
%
%     The optional third argument is a string that indicates the way the 
%	  approximation is defined, possible values are:
% 		'p-norm'   : overestimate using p-norm 
% 		'p-mean'   : underestimate using p-mean 
%		'KS'           : Kreisselmeier-Steinhauser, overestimate
%		'LKS'          : Kreisselmeier-Steinhauser, underestimate
%       'average'      : This is provided only for testing, it is not an
%                        actual max.

    % Do a shift to prevent overflow
    epxm = @(x) exp(p*(x-max(x)));
    sum_epxm = @(x) sum(epxm(x));
    switch form_def
        case 'p-norm'
            S = sum(x.^p).^(1/p);
            dSdx = (x./S).^(p-1);
        case 'p-mean'
            N = size(x,1);
            S = (sum(x.^p)/N).^(1/p);
            dSdx = (1/N)*(x./S).^(p-1);            
        case 'KS'
            S = max(x)+log(sum_epxm(x))/p;
            dSdx = epxm(x)./sum_epxm(x);
        case 'LKS'
            % note: convergence might be fixed with Euler-Gamma
            N = length(x);
            S = max(x) + log(sum_epxm(x)/N)/p; 
            dSdx =epxm(x)./sum_epxm(x);
        case 'avgKS'
            N = length(x);  
            S = max(x)+0.5*(log(sum_epxm(x)) + log(sum_epxm(x)/N))/p;
            dSdx = epxm(x)./sum_epxm(x);     
        case 'softmax'
            S = sum(x.*epxm(x))/sum(epxm(x));
            dSdx = epxm(x).*(1+p*(x - S))/sum_epxm(x);
        case 'average' 
            N = length(x);
            S = sum(x)/N;
            dSdx = 1/N*ones(length(x),1);
        otherwise
            error('smooth_max received invalid form_def.')
    end


end