function f = ant_fnsf(x1,y1,x2,y2,m3,f_init)
% This is the code to estimate the fundamental matrix using
%  the unconstrained vesrion of the fundamental numerical scheme.
%  It has been written so as to be included in Phil Torr's
%  structure and motion package (see http://research.microsoft.com/downloads/)
%  for the purposes of comparing estimators.  If you really want
%  good estimates of F using FNS, however, you should use the 
%  constrained estimator (which I will write next and should be
%  labelled something like ant_cfnsf.m)
% For a description of the method see 
%  http://www.cs.adelaide.edu.au/users/hengel/Vision/ParameterEstimation/index.html
% Anton van den Hengel - anton@cs.adelaide.edu.au

if nargin < 5 | nargin > 6
    error('Wrong number of args to fns');
elseif nargin == 5
    % We don't have an initial estimate so we need to generate one
    t = torr_estf(x1,y1,x2,y2,length(x1),m3);
else
    t = f_init;
end

iter = 0; bestResSoFar = 1e10;
t_old = ones(size(t));
while ((iter < 5) & (norm(t - t_old) > 1e-8))
    iter = iter + 1;
    t_old = t;
    
    [M,X] = fnsJamlMatrices(t,x1,y1,x2,y2,m3);
    
    % I'm keeping the best estimate (in terms of J_AML) as it goes by
    %  this just adds a little bit more stability.  If you're in a 
    %  hurry you can remove this bit and reduce the number of iterations
    %iter
    thisRes = t'*M*t;
    if thisRes < bestResSoFar
        bestT = t;
        bestResSoFar = thisRes;
    end
    
    [v, l] = eig(X);
    
    t = v(:, 1);
%    t = froNormalise(t);
    t = t./norm(t,'fro');
    
end

    M = fnsJamlMatrices(t,x1,y1,x2,y2,m3);    
    thisRes = t'*M*t;
    if thisRes < bestResSoFar
        bestT = t;
    end

    f = bestT;

return



function [M,X] = fnsJamlMatrices(theta,x1,y1,x2,y2,m3)
%
% Computes M_theta and X_theta (see the paper for details)
%  It is possible to make this faster, but at the cost of 
%  making it harder to understand
%  Also I wouldn't usually use m3, but if we don't use it
%  here then the error measures in the rest of this package
%  need to be modified (or really just have m3 set to 1).

t = reshape(theta', 9, 1);
Xbits = zeros(9,9);
M = zeros(9,9);

u = zeros(9,1);
du = zeros(9, 4);    du(7, 1) = m3;    du(8, 2) = m3;    du(3, 3) = m3;    du(6, 4) = m3;

for i = 1:length(x1),
    u = [x1(i) * x2(i), y1(i) * x2(i), x2(i)*m3, x1(i) * y2(i), ...
            y1(i) * y2(i), y2(i)*m3, x1(i)*m3, y1(i)*m3, m3^2]';
    
    A = u * u';
    % Maybe it would be faster to do these as a block?
    du(1, 1) = x2(i);
    du(4, 1) = y2(i);
    
    du(2, 2) = x2(i);
    du(5, 2) = y2(i);
    
    du(1, 3) = x1(i);
    du(2, 3) = y1(i);
    
    du(4, 4) = x1(i);
    du(5, 4) = y1(i);
    
    % use identity covs ...
    B = du * du';
    theta_b_theta = t' * B * t;
    M = M + (A / theta_b_theta);
    if nargout > 1
        Xbits = Xbits - B * ((t' * A * t) / (theta_b_theta^2));
    end
end

if nargout > 1
    X = M + Xbits;
end

