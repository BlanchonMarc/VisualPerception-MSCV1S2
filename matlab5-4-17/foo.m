clear all
close all

X = rand(3, 12);
P = rand(3,4);
x = projection(P,X);

[K,R,t] = art(P); %extract K,R, t from P
x = inv(K) * [x; ones(1, size(x,2))];  % normalized image points
X = rigid([R,t],X); % transform 3D point into camera frame

for i = 1: size(x,2)
    plot3([x(1,i), X(1,i)],  [x(2,i), X(2,i)], [x(3,i), X(3,i)])
    hold on    
end
