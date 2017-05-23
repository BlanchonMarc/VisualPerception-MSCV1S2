%funmatgradeig Matriu Fonamental metode del gradien amb els valors propis
%    [F] = funmatgradeig(M,n,error) 
%
%    M matriu de 4-per-n amb les coordendes de n punts en dues imatges
%      on n és mes gran o igual que 8
%        primera fila: coordenada x de la primera imatge
%        segona fila:  coordenada y de la primera imatge
%        tercera fila: coordenada x de la segona imatge
%        quarta fila:  coordenada y de la segona imatge
%    n 		Numero de iteracións maxima
%    error 	Error maxim de la matriu estimada
%
%    F matriu fonamental de 3-per-3 amb l'origen de coordenades del món a la segona camera
%
% by X. Armangue
% (c) Mr3D - University of Girona, September 2002
%
function [F]=funmatgradeig(M,n,error)

if (size(M,1)~=4) | (size(M,2)<8),
   disp('Error: parametres incorrectes')
else
   Fini=funmateig(M);
   
 	U=[];
	for i=1:size(M,2),
	   U=[U; M(1,i)*M(3,i) M(1,i)*M(4,i) M(1,i) M(2,i)*M(3,i) M(2,i)*M(4,i) M(2,i) M(3,i) M(4,i) 1];      
   end
   
   w=[];
   F=Fini;
   fant=[0; 0; 0; 0; 0; 0; 0; 0;0];
   for j=1:n,
      for i=1:size(M,2),
      	m1 = [M(1,i) ; M(2,i); 1];
	      m2 = [M(3,i) ; M(4,i); 1];
	      l1=F*m2;
         l2=F'*m1;
         g(i)=sqrt(l1(1)^2+l1(2)^2+l2(1)^2+l2(2)^2);
   	   UU(i,:)=U(i,:)./g(i);
      end
      
      [V,D]=eig(UU'*UU);
   
   	[minim,ind]=min(sum(D));
		f=V(:,ind);
      F=[f(1) f(2) f(3); f(4) f(5) f(6); f(7) f(8) f(9)];
      
      F=funmatrank2(F);

      if max(abs(f-fant))<error,
         break;
      end
      fant=f;
   end
   
end





