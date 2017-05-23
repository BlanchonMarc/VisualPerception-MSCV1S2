function x = projection( P, X )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

X = [X; ones(1, size(X,2))];

x = P * X;

x = x ./ repmat(x(3,:),3,1);

x = x(1:2, :);

end
