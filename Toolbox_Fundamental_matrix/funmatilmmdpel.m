%funmatilmmdpel  Iterative Linear Method Minimizing Distances of Points to Epipolar Lines 
%    [F] = funmatilmmdpel(M,Fini,n,error) 
%
%    M matriu de 4-per-n amb les coordendes de n punts en dues imatges
%      on n és mes gran o igual que 8
%        primera fila: coordenada x de la primera imatge
%        segona fila:  coordenada y de la primera imatge
%        tercera fila: coordenada x de la segona imatge
%        quarta fila:  coordenada y de la segona imatge
%    Fini 	Matriu fonamental inicial
%    n 		numero de iteracións maxima
%    error 	error maxim de la matriu estimada
%
%    F 	matriu fonamental de 3-per-3 amb l'origen de coordenades del món a la segona camera
%
% by X. Armangue
% (c) Mr3D - University of Girona, September 2002
%
function [F]=funmatilmmdpel(M,Fini,n,error)

if (size(M,1)~=4) | (size(M,2)<8),
   disp('Error: parametres incorrectes')
else
   F=Fini;
   
  	Up=[];
	for i=1:size(M,2),
	   Up=[Up; M(1,i)*M(3,i) M(1,i)*M(4,i) M(1,i) M(2,i)*M(3,i) M(2,i)*M(4,i) M(2,i) M(3,i) M(4,i)];      
   end
   
   w=[];
   UU=Up;
   fpant=[0; 0; 0; 0; 0; 0; 0; 0];
   for j=1:n,
      for i=1:size(M,2),
      	m1 = [M(1,i) ; M(2,i); 1];
	      m2 = [M(3,i) ; M(4,i); 1];
	      l1=F*m2;
	      l2=F'*m1;
	      w(i)=sqrt((l1(1)^2+l1(2)^2+l2(1)^2+l2(2)^2)/((l1(1)^2+l1(2)^2)*(l2(1)^2+l2(2)^2)));
   	   UU(i,:)=Up(i,:).*w(i);
      end
      fp=-(UU'*UU)^(-1)*UU'*w';
      F=[fp(1) fp(2) fp(3); fp(4) fp(5) fp(6); fp(7) fp(8) 1];
      if max(abs(fp-fpant))<error,
         break;
      end
      fpant=fp;
   end
   F=[fp(1) fp(2) fp(3); fp(4) fp(5) fp(6); fp(7) fp(8) 1];
end





