%funmatgradls Matriu Fonamental metode del gradien amb minims quadrats
%    [F] = funmatgradls(M,m,error) 
%
%    M matriu de 4-per-n amb les coordendes de n punts en dues imatges
%      on n és mes gran o igual que 8
%        primera fila: coordenada x de la primera imatge
%        segona fila:  coordenada y de la primera imatge
%        tercera fila: coordenada x de la segona imatge
%        quarta fila:  coordenada y de la segona imatge
%    n 		Numero de iteracions maxim
%    error 	Error maxim de la matriu estimada
%
%    F matriu fonamental de 3-per-3 amb l'origen de coordenades del món a la segona camera
%
% by X. Armangue
% (c) Mr3D - University of Girona, September 2002
%
function [F]=funmatgradls(M,n,error)

if (size(M,1)~=4) | (size(M,2)<8),
   disp('Error: parametres incorrectes')
else
   Fini=funmatls(M);
   
 	Up=[];
	for i=1:size(M,2),
	   Up=[Up; M(1,i)*M(3,i) M(1,i)*M(4,i) M(1,i) M(2,i)*M(3,i) M(2,i)*M(4,i) M(2,i) M(3,i) M(4,i)];      
   end
   
   w=[];
   F=Fini;
   UU=Up;
   fpant=[0; 0; 0; 0; 0; 0; 0; 0];
   for j=1:n,
      for i=1:size(M,2),
      	m1 = [M(1,i) ; M(2,i); 1];
	      m2 = [M(3,i) ; M(4,i); 1];
	      l1=F*m2;
         l2=F'*m1;
         g(i)=sqrt(l1(1)^2+l1(2)^2+l2(1)^2+l2(2)^2);
   	   UU(i,:)=Up(i,:)./g(i);
      end
      fp=-(UU'*UU)^(-1)*UU'*(ones(size(M,2),1)./g');
      F=[fp(1) fp(2) fp(3); fp(4) fp(5) fp(6); fp(7) fp(8) 1];
      if max(abs(fp-fpant))<error,
         break;
      end
      fpant=fp;
   end
   F=[fp(1) fp(2) fp(3); fp(4) fp(5) fp(6); fp(7) fp(8) 1];
end





