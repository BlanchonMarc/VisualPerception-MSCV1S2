function X = star( x )
% STAR  Returns the skew-symmetric matrix X s.t. cross(x,y) = X*y

 X=[   0    -x(3)  x(2)
      x(3)    0   -x(1)
     -x(2)  x(1)   0   ];

end

