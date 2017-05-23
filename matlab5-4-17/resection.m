function P = resection( x, X )

X = [X; ones(1, size(X,2))];
x = [x; ones(1, size(x,2))];

L = [];
for i = 1: size(X,2)
    L = [L; kron( X(:,i)', star(x(:,i)) )];
end

[~,~,V] = svd(L);
P = reshape(V(:,end),3,4);

