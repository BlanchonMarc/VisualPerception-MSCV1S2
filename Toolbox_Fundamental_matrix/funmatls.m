%funmatls    Matriu Fonemantal amb el metode dels 8 o més punts utilitzant minims quadrats
%    [F] = funmatls(M) 
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
function [F]=funmatls(M)

if (size(M,1)~=4) | (size(M,2)<8),
   disp('Error: parametres incorrectes')
else
	Up=[];
	for i=1:size(M,2),
	   Up=[Up; M(1,i)*M(3,i) M(1,i)*M(4,i) M(1,i) M(2,i)*M(3,i) M(2,i)*M(4,i) M(2,i) M(3,i) M(4,i)];      
   end

	fp=-(Up'*Up)^(-1)*Up'*ones(size(M,2),1);
	
   F=[fp(1) fp(2) fp(3); fp(4) fp(5) fp(6); fp(7) fp(8) 1];
end





