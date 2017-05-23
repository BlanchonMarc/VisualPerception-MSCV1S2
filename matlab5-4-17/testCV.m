clear all

% one view stuff

X = rand(3, 12);
P = rand(3,4);
x = projection(P,X);

P_est = resection(x, X);
rank([P(:), P_est(:)])

x_est = projection(P_est,X);

norm(x-x_est)

% two views stuff

% left
P1 = rand(3,4);
x1 = projection(P1,X);

% right
P2 = rand(3,4);
x2 = projection(P2,X);

X_est = f_intersection(P1,P2, x1, x2);

norm(X-X_est)










