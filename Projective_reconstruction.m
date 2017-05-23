%% PROJECTIVE RECONSTRUCTION
clear all; close all; clc;

first=1;

if first == 1
fprintf('If this is the first time running the code make sure to change the folders path line 309 and 319 and then change the value of variable in confition line 3');
pause;
end

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

%% Projection

%Projection 1
P1 = Camera_one * Scene;                %Multiply the scene by the camera matrix
P1(1,:) = P1(1,:) ./ P1(3,:);           %Divide all lines by the third one
P1(2,:) = P1(2,:) ./ P1(3,:);           %Values should be up to scale
P1(3,:) = P1(3,:) ./ P1(3,:);           %Third line should be scale factor (1)

fprintf('Projection from Camera One Computed \n');

%Pojection 2
P2 = Camera_two * Scene;
P2(1,:) = P2(1,:) ./ P2(3,:);
P2(2,:) = P2(2,:) ./ P2(3,:);
P2(3,:) = P2(2,:) ./ P2(3,:);

fprintf('Projection from Camera Two Computed \n');


%% Plot the 2D Projection for each cameras

scattering = 1;


if scattering ~= 1
    if randomize == 1 %if teapot need to print only small points because filled gives ba render
        figure
        
        sc1 = plot(P1(1,:), P1(2,:), 'o');
        
        grid on
        xlabel('x');ylabel('y');
        title('2D Projection for First Camera');
        legend('Camera 1' );
        
        
        figure
        
        sc2 = plot(P2(1,:), P2(2,:), 'or');
        
        grid on
        xlabel('x');ylabel('y');
        title('2D Projection for Second Camera');
        legend('Camera 2' );
    else
        figure
        
        sc1 = plot(P1(1,:), P1(2,:), '.');
        
        grid on
        xlabel('x');ylabel('y');
        title('2D Projection for First Camera');
        legend('Camera 1' );
        
        
        figure
        
        sc2 = plot(P2(1,:), P2(2,:), '.r');
        
        grid on
        xlabel('x');ylabel('y');
        title('2D Projection for Second Camera');
        legend('Camera 2' );
    end
else

    
    if randomize == 1 %if teapot need to print only small points because filled gives ba render
        figure
        
        sc1 = plot(P1(1,:), P1(2,:), 'filled');
        
        grid on
        xlabel('x');ylabel('y');
        title('2D Projection for First Camera');
        legend('Camera 1' );
        
        
        figure
        
        sc2 = plot(P2(1,:), P2(2,:), 'filled');
        
        grid on
        xlabel('x');ylabel('y');
        title('2D Projection for Second Camera');
        legend('Camera 2' );
    else
        figure
        
        sc1 = scatter(P1(1,:), P1(2,:), '.');
        
        grid on
        xlabel('x');ylabel('y');
        title('2D Projection for First Camera');
        legend('Camera 1' );
        
        
        figure
        
        sc2 = scatter(P2(1,:), P2(2,:), '.');
        
        grid on
        xlabel('x');ylabel('y');
        title('2D Projection for Second Camera');
        legend('Camera 2' );
    end
    
    
end

%% Part 2

%% Conversion of translation to apply cross product
T=Translation;                           %Copy translation between two camera ( for us the second camera translation till the first camera is foxed at 0 0)
T_cross = [0 -T(3) T(2);                 %Transformation
           T(3) 0 -T(1);
           -T(2) T(1) 0];
R = -2 .*Rotation;                       %till second rotation is inverse of first one, the general rotation between 1 and 2 is two time the inverse of the second rotation
E = T_cross * R;                         %Essential Matrix


%Estimation of fundamental matrix using intrinsec parameters
%E = (K')T * F * K
%F = (K'T)-1 * E * K-1
fprintf('Fundamental Matrix Computation \n');
F = inv(transpose(fandpos_Matrix)) * E * inv(fandpos_Matrix)


%% Part 3

%% Fundamental Matrix from Two Images ( Matlab functions )

% Adaptation of the projections in M x 2 dimension
Pm1 = [P1(1 , :); P1(2 , :)]'; 
Pm2 = [P2(1 , :); P2(2 , :)]';

%Example taken from doc estimateFundamentalMatrix

%Fundamental matrix estimation using ransac ( matlab function )
fprintf('Fundamental Matrix Computation with RANSAC (Matlab) \n');
fRANSAC = estimateFundamentalMatrix(Pm1,Pm2,'Method', 'RANSAC', 'NumTrials', 2000, 'DistanceThreshold', 1e-4)

% Fundamental matrix estimation using no specific parameters but ability to
% extract the inliers
[FMatlab,inliersIndex] = estimateFundamentalMatrix(Pm1,Pm2);

fprintf('Fundamental Matrix Computation : Succes \n');
FMatlab

fprintf('Index Inliers Computation : Succes \n');
inliersIndex


%Trying to plot and siplay the correspondancies, not working, images are
%too big even with warning desactivation and auto adjust option for images
isWorking = 0;
if isWorking ~= 0
    warning('off', 'Images:initSize:adjustingMag');
    figure;
    showMatchedFeatures(Pm1(:), Pm2(:), Pm1(inliersIndex,:),Pm2(inliersIndex,:),'montage','PlotOptions',{'ro','go','y--'});
    title('Point matches after outliers were removed');
end


%% Fundamental Matrix from Two Images ( Slavi's Toolbox )

% In Progress ----------------------------------------------------------

%add path to the toolbox to avoid copy paste everything and keep folder
%clean
addpath(genpath('/Users/marc/Documents/MsCV/Semester 2/Visual Perception/Lab/Toolbox_Fundamental_matrix'));




% END In Progress ------------------------------------------------------


%% 3D Scene from cannonical representation

addpath(genpath('/Users/marc/Documents/MsCV/Semester 2/Visual Perception/Lab/matlab5-4-17'));

fail = 0;%Try the fail for the next computation

% Computation of the estimation using invited teacher code
try
    [U,S,V]=svd(FMatlab);
    epipole=V(:,end);
    translation = [0 -epipole(3) epipole(2); epipole(3) 0 -epipole(1); -epipole(2) epipole(1) 0];
    M=FMatlab*inv(translation)
    
    

    P1est=[eye(3) zeros(3,1)];
    P2est=[M epipole];

    X_est = f_intersection(P1est, P2est, Pm1', Pm2');

    figure
    sc3 = scatter3(X_est(1,:),X_est(2,:),X_est(3,:), '.');
    title('Estimated 3D Scene');
    fail=0;
    fprintf('Matrix Computation : Succes \n');
catch
    fail=1;
    fprintf('Matrix Computation : FAIL !! Rerun the code \n');
end

%% 2D Residual Error

if fail == 0
fprintf('Residual Error \n');
ResError = sqrt( mse(X_est - Scene(1:3,:)) )


%% Comparison

figure
PTCloud = pcread('teapot.ply');
Scene2 = [ PTCloud.Location' ; zeros(1 , size(PTCloud.Location,1))];
hold on
scatter3(X_est(1,:),X_est(2,:),X_est(3,:), 'r.'); %VERY VERY VERY SMALL
scatter3(Scene2(1,:),Scene2(2,:),Scene2(3,:), 'g.');
title('Try comparison between estimation and base model')
hold off

else
    fprintf('Computation fail from Fundamental Matrix Computed from Matlab, Relaunch the code !! \n');
end



