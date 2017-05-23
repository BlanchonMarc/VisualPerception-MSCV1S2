%funmatnr    Matriu Fonemantal amb el metode dels 8 o més punts utilitzant el metode iteratiu de Newton-Raphson
%    [F] = funmatnr(M,Fini,n,error) 
%
%    M matriu de 4-per-n amb les coordendes de n punts en dues imatges
%      on n és mes gran o igual que 8
%        primera fila: coordenada x de la primera imatge
%        segona fila:  coordenada y de la primera imatge
%        tercera fila: coordenada x de la segona imatge
%        quarta fila:  coordenada y de la segona imatge
%    Fini 	Matriu fonamental inicial
%    n 		Numero de iteracións maxim
%    error 	Error maxim de la matriu estimada
%
%    F matriu fonamental de 3-per-3 amb l'origen de coordenades del món a la segona camera
%
% by X. Armangue
% (c) Mr3D - University of Girona, September 2002
%
function [F]=funmatnr(M,Fini,n,error)

if (size(M,1)~=4) | (size(M,2)<8) | (size(Fini,1)~=3) | (size(Fini,1)~=3) | n<=0,
   disp('Error: parametres incorrectes')
else
   
   vf11=Fini(1,1);
	vf12=Fini(1,2);
	vf13=Fini(1,3);
	vf21=Fini(2,1);
	vf22=Fini(2,2);
	vf23=Fini(2,3);
	vf31=Fini(3,1);
	vf32=Fini(3,2);
	vf33=Fini(3,3);
	
	VG=[];
	J=[];
	dk=[];
	VG2=[];
	JT=[];
	JTJ=[];
	IJTJ=[];
	IJTJJT=[];
	CV=[];
	
		
	for k=1:n,
	   for i=1:size(M,2),
	      vu1=M(1,i); vv1=M(2,i); vu2=M(3,i); vv2=M(4,i);
       
         VG(i)=((vu1*vf11+vv1*vf21+vf31)*vu2+(vu1*vf12+vv1*vf22+vf32)*vv2+vu1*vf13+vv1*vf23+vf33)';
	      J(i,1)=vu1*vu2;
	      J(i,2)=vu1*vv2;
	      J(i,3)=vu1;
	      J(i,4)=vv1*vu2;
	      J(i,5)=vv1*vv2;
	      J(i,6)=vv1;
	      J(i,7)=vu2;
	      J(i,8)=vv2;
	      J(i,9)=1;
	   end;
	   nVGk= VG*VG';
		JT = J';
		JTJ = JT*J;
		IJTJ=inv(JTJ);
		IJTJJT=IJTJ*JT;
	   dk =IJTJJT*VG';
   
	   xk=[vf11,vf12,vf13,vf21,vf22,vf23,vf31,vf32,vf33]';
	   j=0;
	   acabar=0;  
	   while acabar < 1,
	      xk2=xk - 2^(-j).*dk;
	      for i=1:size(M,2),
	         vu1=M(1,i); vv1=M(2,i); vu2=M(3,i); vv2=M(4,i);
	         vf11=xk2(1); vf12=xk2(2); vf13=xk2(3); vf21=xk2(4); vf22=xk2(5); vf23=xk2(6); vf31=xk2(7); vf32=xk2(8); vf33=xk2(9);  
		      VG2(i)=((vu1*vf11+vv1*vf21+vf31)*vu2+(vu1*vf12+vv1*vf22+vf32)*vv2+vu1*vf13+vv1*vf23+vf33)';
	         i=i+1;
	      end;
	      nVG2k=VG2'*VG2;
	      if nVG2k < nVGk,  
	         acabar=1;
	      end;
	      j=j+1;
	   end;
	   
	   j=j-1;
		xk2=xk - 2^(-j)*dk;
		vf11=xk2(1); vf12=xk2(2); vf13=xk2(3); vf21=xk2(4); vf22=xk2(5); vf23=xk2(6); vf31=xk2(7); vf32=xk2(8); vf33=xk2(9);
	   CV(k,1)=vf11; CV(k,2)=vf12; CV(k,3)=vf13; CV(k,4)=vf21; CV(k,5)=vf22; CV(k,6)=vf23; CV(k,7)=vf31; CV(k,8)=vf32; CV(k,9)=vf33; CV(k,10)=nVGk;
	   if nVGk < error, 
	      break;
	   end;
	end
   F=[vf11 vf12 vf13; vf21 vf22 vf23; vf31 vf32 vf33];
end
