function t = ant_cfnsf(x1,y1,x2,y2,m3,f_init)
% This is the code to estimate the fundamental matrix using
%  the constrained vesrion of the fundamental numerical scheme.
%  It has been written so as to be included in Phil Torr's
%  structure and motion package (see http://research.microsoft.com/downloads/)
%  for the purposes of comparing estimators.  
% For a description of the method see 
%  http://www.cs.adelaide.edu.au/users/hengel/Vision/ParameterEstimation/index.html
% Anton van den Hengel - anton@cs.adelaide.edu.au

% The m3 thing makes things a little bit more complicated.  The error measures all use
%  m3 set to some number which usually isn't 1.  So we need to do the estimation
%  taking that into account.  The problem is that we need to (Hartley) normalise the
%  data before we do cfns or it won't work, and I'm not sure that the m3 thing is equivalent
%  as it's only a scaling and not a recentering, and I haven't tested it.
[x1,y1,x2,y2] = HartleyNormaliseFData(x1,y1,x2,y2);

if nargin < 5 | nargin > 6
    error('Wrong number of args to fns');
elseif nargin == 5
    % We don't have an initial estimate so we need to generate one
        t = torr_estf(x1,y1,x2,y2,length(x1),1);
else
    % Need to be careful with initial estimates because they have to 
    %  be calculated with m3 = 1
    t = f_init;
end

iter = 0;
t_old = ones(size(t));
while ((iter < 5) & (norm(t - t_old) > 1e-8))
    iter = iter + 1;
    t_old = t;
    
    ZZ = constFnsFactorMatrix(t,x1,y1,x2,y2);
    %[M,ZZ] = fnsJamlMatrices(t,x1,y1,x2,y2);
    [v, l] = eig(ZZ);
    l = abs(diag(l));
    [l, indices] = sort(l);
    t = v(:, indices(1));
    
    % Normalise just for numerical stability
    t = t/norm(t,'fro');
    
end

% Undo the Hartley Normalisation
t = unNormalise(t);

% And put the result into the right form re m3
T = eye(3);
T(3,3)=m3;
PhilT=inv(T);
t = reshape((PhilT' * reshape(t,3,3)' * PhilT)',9,1);
t = t/norm(t,'fro');
return



function [M,X,H] = fnsJamlMatrices(t,x1,y1,x2,y2)
%
% Computes M_theta and X_theta (see the paper for details)
%  It is possible to make this faster, but at the cost of 
%  making it harder to understand

Xbits = zeros(9,9);
M = zeros(9,9);
T = zeros(9,9);
tt = t*t'; % Just to speed things up a little
u = zeros(9,1);
du = zeros(9, 4);    du(7, 1) = 1;    du(8, 2) = 1;    du(3, 3) = 1;    du(6, 4) = 1;

for i = 1:length(x1),
    u = [x1(i) * x2(i), y1(i) * x2(i), x2(i), x1(i) * y2(i), ...
            y1(i) * y2(i), y2(i), x1(i), y1(i), 1]';
    
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
        tAt = t' * A * t;
        Xbits = Xbits - B * ((tAt) / (theta_b_theta^2));
        if nargout == 3
            T = T + (2/theta_b_theta^2)*(A*tt*B + B*tt*A - (2/theta_b_theta)*(tAt * B*tt*B));
        end
    end
end

if nargout > 1
    X = M + Xbits;
    if nargout == 3
        H = 2*(X-T);
    end
end

return



function [ZZ] = constFnsFactorMatrix(t,x1,y1,x2,y2)
%
% Computes Z'Z_theta 
% 

tt = t*t';
e = eye(9);

a = dphif(t);
aa = a*a';
norma = norm(a);
Phi = H_phi_f(t);

[M,X,H] = fnsJamlMatrices(t,x1,y1,x2,y2);

P = eye(9)-aa*norma^-2;

A_t = P*H*(2*tt-norm(t)^2*eye(9));

% Now we construct B_t from its B_tits
B_t = zeros(9, 9);
for i =1:9
    B_t = B_t + (Phi*e(:,i)*a' + a*e(:,i)'*Phi)*X*t*e(:,i)';
end
B_t = B_t - 2*norma^-2*aa*X*t*a'*Phi;
B_t = norm(t)^2*(norma^-2)*B_t;

% Construct C
phi_t = phiF(t);
C_t = 3*norma^-2*((phi_t/4)*Phi + aa - (phi_t/2)*norma^-2*aa*Phi);

% And finaly Z
Z = A_t + B_t + C_t;

% And Z'Z
ZZ = Z'*Z;

return



function p = phiF(f)
%
% Makes constraint psi from funmatrix in F or theta form
%
%  t = froNormalise(vectorise(f));
  t = f./norm(f);
  p = (+ t(1) * (t(5)*t(9) - t(6)*t(8)) ...
       - t(2) * (t(4)*t(9) - t(6)*t(7)) ...
       + t(3) * (t(4)*t(8) - t(7)*t(5)) );

function a = dphif(t)

% a_phi is the Jacobian of phi for F estimation (d Phi / d theta)
a = [t(5)*t(9)-t(6)*t(8);
    t(6)*t(7)-t(4)*t(9);
    t(4)*t(8)-t(5)*t(7);
    t(3)*t(8)-t(2)*t(9);
    t(1)*t(9)-t(3)*t(7);
    t(2)*t(7)-t(1)*t(8);
    t(2)*t(6)-t(3)*t(5);
    t(3)*t(4)-t(1)*t(6);
    t(1)*t(5)-t(2)*t(4)]./2;
return


function Phi = H_phi_f(t)

%Calculate the Hessian of \phi at \Theta
PhiBit = [0 0 0     0  t(9) -t(8)    0  -t(6)  t(5);
    0 0 0 -t(9)    0   t(7)  t(6)    0  -t(4);
    0 0 0  t(8) -t(7)    0  -t(5)  t(4)    0 ;
    0 0 0    0     0     0     0   t(3) -t(2);
    0 0 0    0     0     0  -t(3)    0   t(1);
    0 0 0    0     0     0   t(2) -t(1)    0 ;
    0 0 0    0     0     0     0     0     0 ;
    0 0 0    0     0     0     0     0     0 ;
    0 0 0    0     0     0     0     0     0 ];
Phi = PhiBit + PhiBit';



function [nx1,ny1,nx2,ny2] = HartleyNormaliseFData(x1,y1,x2,y2)
global NormalisationT

% Hartleyesq normalisation of the data ...
% If a f-matrix F_n is computed from the normalised data, then the
% f-matrix, F, corresponding to the un-normalised data is given by:
% F =  T2' * F_n * T1;

%THis was written for tensor data so I'll just munge the vectors into a tensor
nPoints = length(x1);
data = zeros(2,nPoints,2);
data(1,:,1)=x1;data(1,:,2)=x2;data(2,:,1)=y1;data(2,:,2)=y2;
norm_data = zeros(size(data));

% find centroid of the image points
centers = mean(data,2);
centered = data-repmat(centers,[1,nPoints,1]);
points1 = centered(:,:,1);
points1 = points1 .* points1;
total_dist = sum(sqrt(sum(points1)));
factor = nPoints * sqrt(2) / total_dist;

normedData(:,:,1) = factor*centered(:,:,1);
T1 = diag([factor, factor, 1]) * ...
     [1 0 -centers(1,1,1); 0 1 -centers(2,1,1); 0 0 1];

points2 = centered(:,:,2);
points2 = points2 .* points2;
total_dist = sum(sqrt(sum(points2)));
factor = nPoints * sqrt(2) / total_dist;

normedData(:,:,2) = factor*centered(:,:,2);
T2 = diag([factor, factor, 1]) * ...
     [1 0 -centers(1,1,2); 0 1 -centers(2,1,2); 0 0 1];

%total_dist1 = sum(sum(sqrt(nPoints*std(points,2))))
%scale = sqrt(2)/mean(std(points,2))

NormalisationT(:,:,1) = T1;
NormalisationT(:,:,2) = T2;
nx1=normedData(1,:,1);nx2=normedData(1,:,2);ny1=normedData(2,:,1);ny2=normedData(2,:,2);

return


function normalTheta = unNormalise(theta)
global NormalisationT

normalTheta = reshape((NormalisationT(:,:,2)' * reshape(theta,3,3)' * NormalisationT(:,:,1))',9,1);