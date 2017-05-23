%funmatMestTorr Matriu Fonamental metode M-Estimators proposat per Torr
%    [F] = funmatMestTorr(M) 
%
%    M matriu de 4-per-n amb les coordendes de n punts en dues imatges
%      on n és mes gran o igual que 8
%        primera fila: coordenada x de la primera imatge
%        segona fila:  coordenada y de la primera imatge
%        tercera fila: coordenada x de la segona imatge
%        quarta fila:  coordenada y de la segona imatge
%    n 		numero de iteracions maxim
%    error 	error maxim de la matriu estimada
%
%    F matriu fonamental de 3-per-3 amb l'origen de coordenades del món a la segona camera
%
% by X. Armangue
% (c) Mr3D - University of Girona, September 2002
%
function [F,ga]=funmatMestTorr(M,n,error)

if (size(M,1)~=4) | (size(M,2)<8),
   disp('Error: parametres incorrectes')
else
   w=ones(1,size(M,2));
   ga=ones(1,size(M,2));
   
   Fant=[0 0 0; 0 0 0; 0 0 0];
	U=[];
	for i=1:size(M,2),
  	   U=[U; M(1,i)*M(3,i) M(1,i)*M(4,i) M(1,i) M(2,i)*M(3,i) M(2,i)*M(4,i) M(2,i) M(3,i) M(4,i) 1];
	end
     
   for iter=1:n,      
      ww=w.*ga;
      for i=1:size(M,2),
   	   UU(i,:)=U(i,:).*ww(i);         
      end      
      
		[V,D]=eig(UU'*UU);
      
      [minim,i]=min(sum(D));
		f=V(:,i);
	
      F=[f(1) f(2) f(3); f(4) f(5) f(6); f(7) f(8) f(9)];
      
      F=funmatrank2(F);

      if max(max(abs(F-Fant)))<error,
         break;
      end
      Fant=F;
      
      for i=1:size(M,2),
         r(i)=[M(1:2,i) ; 1]'*F*[M(3:4,i) ; 1];
         x1=M(1,i);
         y1=M(2,i);
         x2=M(3,i);
         y2=M(4,i);
         rx2=F(1,1)*x1+F(2,1)*y1+F(3,1);
         ry2=F(1,2)*x1+F(2,2)*y1+F(3,2);
         rx1=F(1,1)*x2+F(1,2)*y2+F(1,3);
         ry1=F(2,1)*x2+F(2,2)*y2+F(2,3);
         w(i)=sqrt(1/(rx2^2+ry2^2+rx1^2+ry1^2));
         d(i)=w(i)*r(i);
      end;
      
      sig=median(abs(d))/0.6745;
      
      for i=1:size(M,2),
         if abs(d(i))<sig,
            ga(i)=1;
         elseif abs(d(i))<3*sig,
            ga(i)=sig/abs(d(i));
         else
            ga(i)=0;
         end
      end
   end
end

   
