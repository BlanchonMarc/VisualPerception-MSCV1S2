%funmatrank2   Matriu Fonemantal imposant que rank(F)=2
%    [F] = funmatrank2(Fini) 
%
%    Fini Matriu fonamental inicial
%
%    F matriu fonamental de 3-per-3 amb rank(F)=2 (det(F)=0)
%
% by X. Armangue
% (c) Mr3D - University of Girona, September 2002
%
function [F]=funmatrank2(Fini)

if (size(Fini,1)~=3) | (size(Fini,2)~=3),
   disp('Error: parametres incorrectes')
else
   [U,S,V] = svd(Fini);
   S(3,3)=0;
   F = U*S*V';
end

