%funmatError    Calcul de l'error amb la distància del punts a les seves linies epipolars
%    [mean,stddev,minimum,maximum]=distpuntslin(M,F) 
%
%    M matriu de 4-per-7 amb les coordendes de 7 punts en dues imatges
%        primera fila: coordenada x de la primera imatge
%        segona fila:  coordenada y de la primera imatge
%        tercera fila: coordenada x de la segona imatge
%        quarta fila:  coordenada y de la segona imatge
%    F matriu fonamental de 3-per-3 amb l'origen de coordenades del món a la segona camera
%
%    mitjana distància mitjana
%    stddev  desviació tipus de les diatàncies
%    min     distància minima
%    max     distància màxima
%
% by X. Armangue
% (c) Mr3D - University of Girona, September 2002
%
function [mean,stddev,minimum,maximum]=funmatError(M,F)

if (size(M,1)~=4) | (size(F,1)~=3)| (size(F,2)~=3),
   disp('Error: parametres incorrectes')
else
   aux=[];
   for i=1:size(M,2),
      m1 = [M(1,i) ; M(2,i); 1];
      m2 = [M(3,i) ; M(4,i); 1];
      
      l=F*m2;
      d=m1'*l/sqrt(l(1)^2+l(2)^2);
      dmod=abs(d);
      aux=[aux dmod];
      
      l=F'*m1;
      d=m2'*l/sqrt(l(1)^2+l(2)^2);
      dmod=abs(d);    
      aux=[aux dmod];
   end
   aux_2=aux.^2;
	mean=sum(aux)/(size(M,2)*2);
	stddev=sqrt((sum(aux_2)-sum(aux)^2/(size(M,2)*2))/(size(M,2)*2-1));
	maximum=max(aux);
	minimum=min(aux);
end

      
