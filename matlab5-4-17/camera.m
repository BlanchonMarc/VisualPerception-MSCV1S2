function P = camera( K, look, cop, up )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

R(3,:) = ((look - cop)/norm(look - cop))'; %asse Z che passa tra il centro ottico e il look point
R(1,:) = cross(up',R(3,:));
R(1,:) = R(1,:)/norm(R(1,:)); %asse X ortogonale a Z e all'asse verticale
R(2,:) = cross(R(3,:),R(1,:));%asse Y ortogonale al piano XZ


P = K*[R, -R*cop ];


end

