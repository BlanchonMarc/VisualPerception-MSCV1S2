%% PROJECTIVE RECONSTRUCTION
clear all; close all; clc;


%% Part 1

%% Scene creation and scale adding on fourth row

randomize = 2;


if randomize == 1
    Scene = rand(4, 100) * 50;                  %Random 3D points and scale (100)
    Scene(2,:) = Scene(2,:) + 25 ;              %Pushing the scene on Y axis
    ScaleMatrixGeneration = ones(1,100);        %Ones for scale vector

    Scene(4, : )= ScaleMatrixGeneration ;       %Complete scene generation
elseif randomize == 2
        PTCloud = pcread('teapot.ply');
        Scene = [ PTCloud.Location' ; zeros(1 , size(PTCloud.Location,1))];
        Scene(1,:) = Scene(1,:) + 25 ;  %Pushing the scene on X axis
        Scene(2,:) = Scene(2,:) + 25 ;  %Pushing the scene on Y axis
        Scene(3,:) = Scene(3,:) + 25 ;  %Pushing the scene on Z axis
        ScaleMatrixGeneration = ones(1,size(PTCloud.Location,1));        %Ones for scale vector

        Scene(4, : ) = ScaleMatrixGeneration ;       %Complete scene generation
    else
        
        [X , Y , Z] = sphere(16);
        X = X + 25;
        Y = Y + 25;
        Z = Z + 25;
        
        
        X = [X(:)];
        Y = [Y(:)];
        Z = [Z(:)];
        
        Scene = [X' ; Y' ; Z' ; ones(1 , size(X,1))];
        
    
end

%% Stereo System Creation

fxfy = 1;                                  %fx and fy declaration
u_zero = 0;                                %u0
v_zero = 0;                                %v0 

fx = fxfy;
fy=fxfy;

fandpos_Matrix = [fx 0 u_zero;              %Matrix Computation
                  0 fy v_zero;              %with fx fy
                  0 0 1];                   %and u0 v0

Extended_Identity = [1 0 0 0 ; 0 1 0 0 ; 0 0 1 0 ]; %extended identity

%% Camera 1

if randomize == 1
    ViewPosition = [25 25 1.5]';
else
    ViewPosition = [25 25 25+1.5]';
end
CameraPos = [0 0 0]';
Orientation = [0 1 0]';

Rotation = eye(3,3);                           %rotation camera 1
Rotation(3, :) = ((ViewPosition - CameraPos) / norm(ViewPosition - CameraPos))'; % Z axis
Rotation(1, :) = cross(Orientation,Rotation(3,:));
Rotation(1, :) = Rotation(1, :) / norm(Rotation(1,:));
Rotation(2, :) = cross(Rotation(3, :) ,Rotation(1,:));

Translation = [0; 0; 0];                       %translation camera 1
RT_Matrix = [Rotation(1,1) Rotation(1,2) Rotation(1,3) Translation(1);
             Rotation(2,1) Rotation(2,2) Rotation(2,3) Translation(2);
             Rotation(3,1) Rotation(3,2) Rotation(3,3) Translation(3);
             0 0 0 1];                         %Applying to matrix

Camera_one = fandpos_Matrix * Extended_Identity * RT_Matrix;                 %3x4 Pinhole Matrix camera 1

fprintf('Camera one Matrix Computed \n');


%% Camera 2

if randomize == 1
    ViewPosition = [-25 -25 1.5]';
else
    ViewPosition = [-25 -25 25+1.5]';
end
CameraPos = [25 0 0]';
Orientation = [0 1 0]';

Rotation = eye(3,3);                          %rotation camera 2
Rotation(3, :) = ((ViewPosition - CameraPos) / norm(ViewPosition - CameraPos))'; % Z axis
Rotation(1, :) = cross(Orientation,Rotation(3,:));
Rotation(1, :) = Rotation(1, :) / norm(Rotation(1,:));
Rotation(2, :) = cross(Rotation(3, :) ,Rotation(1,:));

Translation = [50; 0; 0];                     %translation camera 2
RT_Matrix = [Rotation(1,1) Rotation(1,2) Rotation(1,3) Translation(1);
             Rotation(2,1) Rotation(2,2) Rotation(2,3) Translation(2);
             Rotation(3,1) Rotation(3,2) Rotation(3,3) Translation(3);
             0 0 0 1];                         %Applying to matrix
  
Camera_two = fandpos_Matrix * Extended_Identity * RT_Matrix;                 %3x4 Pinhole Matrix camera 2


fprintf('Camera Two Matrix Computed \n');


%% Plot Scene + Cameras
figure
hold on
sc1 = scatter3(Camera_one(1,4),Camera_one(2,4),Camera_one(3,4), 'filled');
sc2 = scatter3(Camera_two(1,4),Camera_two(2,4),Camera_two(3,4), 'filled');

if randomize == 1 %if teapot need to print only small points because filled gives ba render
    sc3 = scatter3(Scene(1,:),Scene(2,:),Scene(3,:), 'filled');
else
    sc3 = scatter3(Scene(1,:),Scene(2,:),Scene(3,:), '.');
end


view(-30,10)
grid on
xlabel('x');ylabel('y');zlabel('z');
title('3D points and cameras');
legend([sc1 sc2 sc3] , 'Camera 1' , 'Camera 2' , '3D Scene');
hold off

%% Part 2

%Computing ground truth essential matrix

T = Translation;                           %Copy translation between two camera ( for us the second camera translation till the first camera is foxed at 0 0)
T_cross = [0 -T(3) T(2);                 %Transformation
           T(3) 0 -T(1);
           -T(2) T(1) 0];
R = -2 .*Rotation;                       %till second rotation is inverse of first one, the general rotation between 1 and 2 is two time the inverse of the second rotation
GT_Essential_Matrix = T_cross * R        %Ground Truth Essential Matrix


[u,eps,v] = svd(GT_Essential_Matrix)  %such that E  = U * Eps * V' // Reference: Wiki


% Computing R
%Reference :  http://www.maths.lth.se/matematiklth/personal/calle/datorseende13/notes/forelas6.pdf

W = [0 -1 0 ;   %W a rotation
     1 0 0 ;
     0 0 1];

Z = [0 1 0 ;    %Z a skew symmetric ( square matrix that it's transpose is its negation) meaning -Z = transpose(Z)
     -1 0 0 ;
     0 0 0];

% ZW = diag([1 0 0])
% transpose(ZW) = -diag([1 0 0])

E_test = u * diag([1    1    0]) *v'

Scale_Factor_Ess = GT_Essential_Matrix ./ E_test;

Scale_Factor_Ess(isnan(Scale_Factor_Ess)) = 0;

  
fprintf('Essential Matrix recovered from svd with scale : \n');
Scale_Factor_Ess
Scale_Factor_Value = Scale_Factor_Ess(2,1)

% Till E = S1R1 = S2R2

S1 = -u * Z * u'
S2 = u * Z * u'
R1 = u * W' * v'
R2 = u * W * v'

fprintf('\n Checking if results are good coherent and R are trully rotations ( identity results of Rt * R ) \n')


R1Check = R1'*R1
R2Check = R2'*R2


fprintf('Results are coherent and R correspond to a rotations \n')


% Computation of translation

fprintf('Computation of the estimated translation \n')

T_est = [u(1,3) u(2,3) u(3,3)]'

%% Part 3

fprintf('\n Checking if results are coherent, meaning translation are in the same direction. \n')
fprintf('This is the case the real translation belong on the x axis and we have a scale factor of -50 \n')

%Projections matrix estimation

P1 = [eye(3) [0 0 0]']

P2_one = [u * W * v'  u( : , 3)]
P2_two = [u * W' * v' u( : , 3)]
P2_three = [u * W * v'  (-1 .* u( : , 3))]
P2_four = [u * W' * v'  (-1 .* u( : , 3))]

%Good projection matrix estimation

% Select the correct second camera matrix

P1est=P1;

%To choose the true second matrix, four possible configurations

P2est1=P2_one;
P2est2=P2_two;
P2est3=P2_three;
P2est4=P2_four;


x1 = P1 * Scene;
x1(1,:) = x1(1,:) ./ x1(3,:);           %Divide all lines by the third one
x1(2,:) = x1(2,:) ./ x1(3,:);           %Values should be up to scale
x1(3,:) = x1(3,:) ./ x1(3,:);           %Third line should be scale factor (1)
x1 = x1(1:2 , :);

x2 = P2_one * Scene;
x2(1,:) = x2(1,:) ./ x2(3,:);           %Divide all lines by the third one
x2(2,:) = x2(2,:) ./ x2(3,:);           %Values should be up to scale
x2(3,:) = x2(3,:) ./ x2(3,:);           %Third line should be scale factor (1)
x2 = x2(1:2 , :);

%Test with first pair:
ptnb = 0;

Q1 = triangulate(x1(:,1)',x2(:,1)',P1est',P2est1');
Q1 = [Q1 1];

c1=Q1(3)*Q1(4)
c2=(P2est1*Q1');
c2=c2(3)*Q1(4)

%Check which matrix is the good one for each c1 c2 values
if c1 > 0 && c2 > 0
    ptnb=1;
end

%Test with second pair:

Q1 = triangulate(x1(:,1)',x2(:,1)',P1est',P2est2');
Q1 = [Q1 1];

c1=Q1(3)*Q1(4)
c2=(P2est2*Q1');
c2=c2(3)*Q1(4)

if c1 > 0 && c2 > 0
    ptnb=2;
end

%Test with third pair:

Q1 = triangulate(x1(:,1)',x2(:,1)',P1est',P2est3');
Q1 = [Q1 1];

c1=Q1(3)*Q1(4)
c2=(P2est3*Q1');
c2=c2(3)*Q1(4)

if c1 > 0 && c2 > 0
    ptnb=3;
end

%Test with fourth pair:

Q1 = triangulate(x1(:,1)',x2(:,1)',P1est',P2est4');
Q1 = [Q1 1];

c1=Q1(3)*Q1(4)
c2=(P2est4*Q1');
c2=c2(3)*Q1(4)

if c1 > 0 && c2 > 0
    ptnb=4;
end

fprintf('The Matrix is the good one if both c1 and c2 are positive , the result will be displayed in green automatically according to values \n')


%% Part 4

fprintf('Estimation of projection first camera : Processing \n')
tic
Proj1 = P1 * Scene;
Proj1(1,:) = Proj1(1,:) ./ Proj1(3,:);           %Divide all lines by the third one
Proj1(2,:) = Proj1(2,:) ./ Proj1(3,:);           %Values should be up to scale
Proj1(3,:) = Proj1(3,:) ./ Proj1(3,:);           %Third line should be scale factor (1)
Proj1 = Proj1(1:2 , :)';
toc
fprintf('Estimation of projection first camera : Done \n')
%TEST1
fprintf('Estimation of projection Second camera with first matrix : Processing \n')
tic
Proj2 = P2_one * Scene;
Proj2(1,:) = Proj2(1,:) ./ Proj2(3,:);           %Divide all lines by the third one
Proj2(2,:) = Proj2(2,:) ./ Proj2(3,:);           %Values should be up to scale
Proj2(3,:) = Proj2(3,:) ./ Proj2(3,:);           %Third line should be scale factor (1)
Proj2 = Proj2(1:2 , :)';

PTCloud = triangulate(Proj1 , Proj2 , P1' , P2_one');
toc
fprintf('Estimation of projection Second camera : Done \n')


hold on
figure
subplot(2,2,1);
if ptnb == 1
    scatter3(PTCloud(:,1),PTCloud(:,2),PTCloud(:,3), 'g.')
else
    scatter3(PTCloud(:,1),PTCloud(:,2),PTCloud(:,3), 'r.')
end


%TEST2
fprintf('Estimation of projection Second camera with second matrix : Processing \n')
tic
Proj3 = P2_two * Scene;
Proj3(1,:) = Proj3(1,:) ./ Proj3(3,:);           %Divide all lines by the third one
Proj3(2,:) = Proj3(2,:) ./ Proj3(3,:);           %Values should be up to scale
Proj3(3,:) = Proj3(3,:) ./ Proj3(3,:);           %Third line should be scale factor (1)
Proj3 = Proj3(1:2 , :)';

PTCloud2 = triangulate(Proj1 , Proj3 , P1' , P2_two');
toc
fprintf('Estimation of projection Second camera : Done \n')
tic
subplot(2,2,2);
if ptnb == 2
    scatter3(PTCloud2(:,1),PTCloud2(:,2),PTCloud2(:,3), 'g.')
else
    scatter3(PTCloud2(:,1),PTCloud2(:,2),PTCloud2(:,3), 'r.')
end

%TEST3
fprintf('Estimation of projection Second camera with third matrix : Processing \n')
Proj4 = P2_three * Scene;
Proj4(1,:) = Proj4(1,:) ./ Proj4(3,:);           %Divide all lines by the third one
Proj4(2,:) = Proj4(2,:) ./ Proj4(3,:);           %Values should be up to scale
Proj4(3,:) = Proj4(3,:) ./ Proj4(3,:);           %Third line should be scale factor (1)
Proj4 = Proj4(1:2 , :)';

PTCloud3 = triangulate(Proj1 , Proj4 , P1' , P2_three');
toc
fprintf('Estimation of projection Second camera : Done \n')
subplot(2,2,3);

if ptnb == 3
    scatter3(PTCloud3(:,1),PTCloud3(:,2),PTCloud3(:,3), 'g.')
else
    scatter3(PTCloud3(:,1),PTCloud3(:,2),PTCloud3(:,3), 'r.')
end


%TEST3
fprintf('Estimation of projection Second camera with fourth matrix : Processing \n')
tic
Proj5 = P2_four * Scene;
Proj5(1,:) = Proj5(1,:) ./ Proj5(3,:);           %Divide all lines by the third one
Proj5(2,:) = Proj5(2,:) ./ Proj5(3,:);           %Values should be up to scale
Proj5(3,:) = Proj5(3,:) ./ Proj5(3,:);           %Third line should be scale factor (1)
Proj5 = Proj5(1:2 , :)';

PTCloud3 = triangulate(Proj1 , Proj5 , P1' , P2_four');
toc
fprintf('Estimation of projection Second camera : Done \n')
subplot(2,2,4);

if ptnb == 4
    scatter3(PTCloud3(:,1),PTCloud3(:,2),PTCloud3(:,3), 'g.')
else
    scatter3(PTCloud3(:,1),PTCloud3(:,2),PTCloud3(:,3), 'r.')
end


title('Display of every estimated point cloud, green point cloud is the good one')
hold off


P1GT = [eye(3) [0 0 0]'];
P2GT = [Rotation Translation];

fprintf('GROUND TRUTH : Estimation of projection first camera : Processing \n')
tic
Proj1 = P1GT * Scene;
Proj1(1,:) = Proj1(1,:) ./ Proj1(3,:);           %Divide all lines by the third one
Proj1(2,:) = Proj1(2,:) ./ Proj1(3,:);           %Values should be up to scale
Proj1(3,:) = Proj1(3,:) ./ Proj1(3,:);           %Third line should be scale factor (1)
Proj1 = Proj1(1:2 , :)';
toc
fprintf('GROUND TRUTH : Estimation of projection first camera : Done \n')

fprintf('GROUND TRUTH : projection Second camera: Processing \n')
tic
Proj2 = P2GT * Scene;
Proj2(1,:) = Proj2(1,:) ./ Proj2(3,:);           %Divide all lines by the third one
Proj2(2,:) = Proj2(2,:) ./ Proj2(3,:);           %Values should be up to scale
Proj2(3,:) = Proj2(3,:) ./ Proj2(3,:);           %Third line should be scale factor (1)
Proj2 = Proj2(1:2 , :)';

PTCloud4 = triangulate(Proj1 , Proj2 , P1GT' , P2GT');
fprintf('GROUND TRUTH : projection Second camera : Done \n')
toc
figure
hold on
scatter3(PTCloud4(:,1),PTCloud4(:,2),PTCloud4(:,3), 'b.')
scatter3(PTCloud2(:,1),PTCloud2(:,2),PTCloud2(:,3), 'g.')
title('Display of the estimated Point Cloud and The Real Point Cloud')
hold off


fprintf('The result between Ground truth and estimation is almost the same \n')


choice = input('Do you want to display the therm by therm error? (1/0) \n');

if choice == 1
    Difference = PTCloud4 - PTCloud2
else
    Difference = PTCloud4 - PTCloud2;
end

%% Part 5
% Yamid / Marc / Stack Overflow

%  p are homogeneous coordinates of the first image of size 3 by n
%  q are homogeneous coordinates of the second image of size 3 by n
p = [x1; ones(1, size(x1,2))];
q = [x2; ones(1, size(x2,2))];

% pq=[(p(1,1:5).*q(1,1:5))' (p(2,1:5).*q(1,1:5))' (p(3,1:5).*q(1,1:5))' (p(1,1:5).*q(2,1:5))' (p(2,1:5).*q(2,1:5))' (p(3,1:5).*q(2,1:5))' (p(1,1:5).*q(3,1:5))' (p(2,1:5).*q(3,1:5))' (p(3,1:5).*q(3,1:5))'];
% 
% nullspace=null(pq);

n = size(p);
NPOINTS = 1000;

% set up matrix A such that A*[v1,v2,v3,s1,s2,s3,s4,s5,s6]' = 0
A = zeros(NPOINTS, 9);

for i = 1:NPOINTS
  A(i,:) = kron(p(:,i),q(:,i))';
end

[U,S,V] = svd(A);

% pick the eigenvector corresponding to the smallest eigenvalue
e = V(:,9);
e = (round(1.0e+10*e))*(1.0e-10);
% essential matrix 
Eest = reshape(e, 3, 3);


%% 
[U,S,V] = svd(Eest);

tu=[U(1,3) U(2,3) U(3,3)]';

D=[0 1 0;-1 0 0;0 0 1];

Ra=U*D*V';

Rb=U*D'*V'; %Correct one

%% Select the correct second camera matrix:

P1est=[eye(3) zeros(3,1)];

%To choose the true second matrix, four possible configurations:

P2est1=[Ra tu];
P2est2=[Ra -tu];
P2est3=[Rb tu];
P2est4=[Rb -tu];

%Test with first pair:

Q1 = triangulate(x1(:,1)',x2(:,1)',P1est',P2est1');
Q1 = [Q1 1];

c1=Q1(3)*Q1(4)
c2=(P2est1*Q1');
c2=c2(3)*Q1(4)

%Test with second pair:

Q1 = triangulate(x1(:,1)',x2(:,1)',P1est',P2est2');
Q1 = [Q1 1];

c1=Q1(3)*Q1(4)
c2=(P2est2*Q1');
c2=c2(3)*Q1(4)

%Test with third pair:

Q1 = triangulate(x1(:,1)',x2(:,1)',P1est',P2est3');
Q1 = [Q1 1];

c1=Q1(3)*Q1(4)
c2=(P2est3*Q1');
c2=c2(3)*Q1(4)

%Test with fourth pair:

Q1 = triangulate(x1(:,1)',x2(:,1)',P1est',P2est4');
Q1 = [Q1 1];

c1=Q1(3)*Q1(4)
c2=(P2est4*Q1');
c2=c2(3)*Q1(4)


%% 3D reconstruction
Q = triangulate(x1',x2',P1est',P2est2');
figure;
scatter3(Q(:,1),Q(:,2),Q(:,3),'.');

%To compare with the ground truth:


figure;
hold on
scatter3(Q(:,1),Q(:,2),Q(:,3),'g.');


scatter3(PTCloud2(:,1),PTCloud2(:,2),PTCloud2(:,3), 'b.') % the real ground truth
%scatter3(Qtruth(:,1),Qtruth(:,2),Qtruth(:,3),'r.');
hold off

Difference = PTCloud2 - Q;

Variance = var(Difference)




























