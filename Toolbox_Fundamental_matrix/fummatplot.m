%fummatplot    Dibuixar les linies epipolars 
%    Representacio de la linea epipolar d'un punt donat sobre la primera imatge en el pla imatge de la segona.
%    fummatplot(M,F) 
%
%    M matriu de 4-per-n amb les coordendes de n punts en dues imatges
%        primera fila: coordenada x de la primera imatge
%        segona fila:  coordenada y de la primera imatge
%        tercera fila: coordenada x de la segona imatge
%        quarta fila:  coordenada y de la segona imatge
%    F Matriu fonamental amb l'origen de coordenades del món a la segona camera
%
% by X. Armangue
% (c) Mr3D - University of Girona, September 2002
%
function fummatplot(M,F)

if (size(M,1)~=4) | (size(F,1)~=3)| (size(F,2)~=3),
   disp('Error: parametres incorrectes')
else
   
   P1=[ M(1:2,:) ; ones(1,size(M,2))];
   P2=[ M(3:4,:) ; ones(1,size(M,2))];

	% Representacio Geometria epipolar Camera 2.
	figure;
	hold on; 
%	axis([-512,512,-100,512]);
	axis([-512,1024,-512,1024]);
	axis ij;
	title('Camera 2');
	for i=1:size(M,2),
		lm2(1:3,i)=F'*P1(1:3,i);
%		plot(P2(1,i),P2(2,i),'g*');
		plot(P2(1,i),P2(2,i),'ko');
%		plot([-1024;1024],[(-1024*lm2(1,i)+lm2(3,i))/(-lm2(2,i)),(1024*lm2(1,i)+lm2(3,i))/(-lm2(2,i))]);  
		plot([-1024;1024],[(-1024*lm2(1,i)+lm2(3,i))/(-lm2(2,i)),(1024*lm2(1,i)+lm2(3,i))/(-lm2(2,i))],'k');  
	end;
	zoom on;

	% Representacio Geometria epipolar Camera 1.
	figure;
	hold on; 
%	axis([-512,512,-100,512]);
	axis([-512,1024,-512,1024]);
	axis ij;
	title('Camera 1');
	for i=1:size(M,2),
 		lm1(1:3,i)=F*P2(1:3,i);
%		plot(P2(1,i),P2(2,i),'g*');
		plot(P1(1,i),P1(2,i),'ko');
%		plot([-1024;1024],[(-1024*lm1(1,i)+lm1(3,i))/(-lm1(2,i)),(1024*lm1(1,i)+lm1(3,i))/(-lm1(2,i))]);  
		plot([-1024;1024],[(-1024*lm1(1,i)+lm1(3,i))/(-lm1(2,i)),(1024*lm1(1,i)+lm1(3,i))/(-lm1(2,i))],'k');  
	end;
	zoom on;

end
