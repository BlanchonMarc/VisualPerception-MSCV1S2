%funmatMesteig Matriu Fonamental metode M-Estimators amb valors propis
%    [F] = funmatMesteig(M,n,error) 
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
function [F,w]=funmatMesteig(M,n,error)

if (size(M,1)~=4) | (size(M,2)<8),
   disp('Error: parametres incorrectes')
else
   Fini=funmateig(M);
   c=4.6851;
   
	U=[];
	for i=1:size(M,2),
  	   U=[U; M(1,i)*M(3,i) M(1,i)*M(4,i) M(1,i) M(2,i)*M(3,i) M(2,i)*M(4,i) M(2,i) M(3,i) M(4,i) 1];
	end
      
   F=Fini;
   r=[ones(1,size(M,2))];
   Fant=[0  0  0; 0  0  0; 0 0 0];   
   
   for iter=1:n,
      rant=r;
      for i=1:size(M,2),
	      r(i)=[M(1:2,i) ; 1]'*F*[M(3:4,i) ; 1];      
	   end;
	   
      switch 2, % sigma relao estimada (2)
      case 1,
         sig=std(r);
      case 2,
         sig=1.4826*(1+5/(size(M,2)-8))*median(abs(r));
      end      
	   
      for i=1:size(M,2),
         switch 2, % pesos (-)
         case 1,
	         if abs(r(i))<c*sig,
		         w(i)=(1-(r(i)/c)^2)^2;
		      else
		         w(i)=0;
            end   
         case 2,
            if abs(r(i))<=sig,
	         	w(i)=1;
		      elseif abs(r(i))<=3*sig,
		         w(i)=sig/abs(r(i));
		      else
		         w(i)=0;
            end
         end
         
   	   UU(i,:)=U(i,:).*w(i);         
		end;

		[V,D]=eig(UU'*UU);
      
      [minim,i]=min(sum(D));
		f=V(:,i);
	
	   F=[f(1) f(2) f(3); f(4) f(5) f(6); f(7) f(8) f(9)];
      
      if max(max(abs(F-Fant)))<error,
         break;
      end
      Fant=F;
   end
end

   
