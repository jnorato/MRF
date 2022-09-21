function [mx0, dmx0] = smooth_max_x0(x, p, type, eps)
%
% This function computes a vector with a smooth maximum approximation of 
% the maximum between 0 and each of the components of vector x, with 
% approximation parameter p.  
%
% Type can be 'softmax', 'KS', 'LKS', 'intH'
%
% The function returns the vector mx0 with the smooth maximum approximations,
% and the derivative of dmx0(i) with respect to its argument x(i)
%
% The parameter eps is optional, and it is the shift for the shiftedKS
% function

if ~exist('eps', 'var')
    eps = 0.01;
end

% m is a shift to prevent overflow due to finite arithmetic
m=max(x);

switch type
    case 'softmax' 
        epxm = exp(p*(x-m));
        empm = exp(-p*m);
        mx0 = x.*epxm./(epxm + empm);
        dmx0 = epxm.*(epxm + empm*(1 + p*x))./((epxm + empm).^2);
        
    case 'SiLU' 
        mx0 = x./(1+exp(-x));
        dmx0 = (1 + exp(-x).*(1+x))./(1+exp(-x)).^2;
        
    case 'softplus'
        % Note that for max(x,0), the KS function is called the softplus
        % function in machine learning
        epxm = exp(p*(x-m));
        empm = exp(-p*m);
        mx0 = m + (1/p)*log(epxm+empm);
        % dmx0 = epxm./(epxm + empm);
        dmx0 = 0.5*(1+tanh(p*x/2));

    case 'shiftedKS'
        xs = (1/p)*log(exp(p*eps)-1);      
        epxm = exp(p*(x+xs-m));
        empm = exp(-p*m);
          mx0 = m + (1/p)*log(epxm+empm); 
        dmx0 = epxm./(epxm + empm);
        
    case 'LKS'
        epxm = exp(p*(x-m));
        empm = exp(-p*m);
        mx0 = m + (1/p)*log(0.5*(epxm+empm))    ;
        dmx0 = epxm./(epxm + empm);
        
    case 'intH' % Integrated Heavise
        px = p*x;
        branch1 = px <= -1;
        branch2 = abs(px) < 1;
        branch3 = px >= 1;
        mx0 = branch1 .* 0 + ...
              branch2 .* (1/32/p).*((px + 1).^4).*(px.^2 - 4*px + 5) + ...
              branch3 .* x;
        dmx0 = branch1 .* 0 + ...
               branch2 .* (1/16).*((px + 1).^3).*(3*px.^2 - 9*px + 8) + ...
               branch3 * 1.0;

    case 'ELU' % Exponential linear unit
        branch1 = x > 0;
        branch2 = x<=0;
        mx0 = branch1.*x + p*branch2.*(exp(x)-1);
        dmx0 = branch1.*1.0 + p*branch2.*exp(x);
        
    case 'smoothELU'
        branch1 = x<=0;
        branch2 = (x>0).*(x <= eps);
        branch3 = x>eps;
        mx0 = branch1.*0 + branch2.*((2/eps)*x.^2 - (1/eps^2)*x.^3) + branch3.*x;
        dmx0 = branch1.*0 + branch2.*((4/eps)*x - (3/eps^2)*x.^2) + branch3.*1;
    case 'PiecePoly' % Piecewise polynomial
        eps = 1e-5;
        branch1 = x<-1; % For safety; element stress constraint will never be less than -1
        branch2 = (x>=-1).*(x <= 0); % Linear portion for -1 <= x <= 1
        branch3 = (x>0).*(x<1); % C^2-continuous, polynomial portion
        branch4 = x >= 1; % Linear portion for x >= 1
        mx0 = branch1.*0 + ...
              eps*branch2.*(x+1) + ...
              branch3.*(eps + eps*x + (6-16*eps)*x.^3 + (23*eps-8)*x.^4 + (3-9*eps)*x.^5) + ...
              branch4.*x;
        dmx0 = branch1.*0 + ...
              eps*branch2 + ...
              branch3.*(eps + 3*(6-16*eps)*x.^2 + 4*(23*eps-8)*x.^3 + 5*(3-9*eps)*x.^4) + ...
              branch4.*1;       
           
    otherwise
        disp('Eror: type of smooth max function not recognized.');
end