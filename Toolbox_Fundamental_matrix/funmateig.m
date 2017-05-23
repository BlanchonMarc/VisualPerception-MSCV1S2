%funmateig    Matriu Fonemantal amb el metode dels 8 o més punts utilitzant valors propis
%    [F] = funmateig(M) 
%
%    M matriu de 4-per-n amb les coordendes de n punts en dues imatges
%      on n és mes gran o igual que 8
%        primera fila: coordenada x de la primera imatge
%        segona fila:  coordenada y de la primera imatge
%        tercera fila: coordenada x de la segona imatge
%        quarta fila:  coordenada y de la segona imatge
%
%    F matriu fonamental de 3-per-3 amb l'origen de coordenades del món a la segona camera
%
% by X. Armangue
% (c) Mr3D - University of Girona, September 2002
%
function [F]=funmateig(M)

if (size(M,1)~=4) | (size(M,2)<8),
   disp('Error: parametres incorrectes')
else
	U=[];
	for i=1:size(M,2),
  	   U=[U; M(1,i)*M(3,i) M(1,i)*M(4,i) M(1,i) M(2,i)*M(3,i) M(2,i)*M(4,i) M(2,i) M(3,i) M(4,i) 1];
	end
	
	[V,D]=eig(U'*U);
      
   [minim,ind]=min(sum(D));
   f=V(:,ind);
	
   F=[f(1) f(2) f(3); f(4) f(5) f(6); f(7) f(8) f(9)];
end



