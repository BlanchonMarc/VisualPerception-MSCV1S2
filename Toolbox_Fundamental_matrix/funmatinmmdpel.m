%funmatinmmdpel Iterative Nonlinear Method Minimizing Distances of Points to Epipolar Lines 
%    [F] = funmatinmmdpel(M,Fini,n) 
%
%    M matriu de 4-per-n amb les coordendes de n punts en dues imatges
%      on n és mes gran o igual que 8
%        primera fila: coordenada x de la primera imatge
%        segona fila:  coordenada y de la primera imatge
%        tercera fila: coordenada x de la segona imatge
%        quarta fila:  coordenada y de la segona imatge
%    Fini 	Matriu fonamental inicial
%    n 		numero de iteracións maxima
%
%    F 	matriu fonamental de 3-per-3 amb l'origen de coordenades del món a la segona camera
%
% by X. Armangue
% (c) Mr3D - University of Girona, September 2002
%
function [F]=funmatinmmdpel(M,Fini,n)

if (size(M,1)~=4) | (size(M,2)<8),
   disp('Error: parametres incorrectes')
else
   indexs=[1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];
   vmini=0;
   
%   if exist('minpuntslin2.mat','file'),
%      load minpuntslin2.mat FF JJ;
%   else
%		syms x1 y1 x2 y2 a b c d;
%	   FF=[];
%	   JJ=[];
%	   comb=1;
%	   for in1=1:6,
%			j0= indexs(in1,1);
%		   j1= indexs(in1,2);
%		   j2= indexs(in1,3);
%			for in2=1:6,
%				i0= indexs(in2,1);
%		   	i1= indexs(in2,2);
%		      i2= indexs(in2,3);
%        
%	         FF(i0,j0,comb)=(a*x1+b*y1)*x2+(c*x1+d*y1)*y2;
%		      FF(i0,j1,comb)=-a*x2-c*y2;
%		      FF(i0,j2,comb)=-b*x2-d*y2;
%		      FF(i1,j0,comb)=-a*x1-b*y1;
%		      FF(i1,j1,comb)=a;
%		      FF(i1,j2,comb)=b;
%		      FF(i2,j0,comb)=-c*x1-d*y1;
%		      FF(i2,j1,comb)=c;
%		      FF(i2,j2,comb)=d;
%         
%	         for ff=1:3,
%		      	for fc=1:3,
%		            JJ(1,3*(ff-1)+fc,comb)=diff(FF(ff,fc,comb),x1);
%		         end
%		      end
%		      for ff=1:3,
%		         for fc=1:3,
%		            JJ(2,3*(ff-1)+fc,comb)=diff(FF(ff,fc,comb),y1);
%		         end
%		      end
%		      for ff=1:3,
%		         for fc=1:3,
%		            JJ(3,3*(ff-1)+fc,comb)=diff(FF(ff,fc,comb),x2);
%		         end
%		      end
%		      for ff=1:3,
%		         for fc=1:3,
%		            JJ(4,3*(ff-1)+fc,comb)=diff(FF(ff,fc,comb),y2);
%		         end
%		      end
%		      for ff=1:3,
%		         for fc=1:3,
%		            JJ(5,3*(ff-1)+fc,comb)=diff(FF(ff,fc,comb),a);
%		         end
%	   	   end
%		      for ff=1:3,
%		         for fc=1:3,
%		            JJ(6,3*(ff-1)+fc,comb)=diff(FF(ff,fc,comb),b);
%		         end
%		      end
%		      for ff=1:3,
%		         for fc=1:3,
%		            JJ(7,3*(ff-1)+fc,comb)=diff(FF(ff,fc,comb),c);
%		         end
%		      end
%		      for ff=1:3,
%		         for fc=1:3,
%		            JJ(8,3*(ff-1)+fc,comb)=diff(FF(ff,fc,comb),d);
%		         end
%	         end
%		      comb=comb+1;         
%	      end
%   	end
%		save minpuntslin2.mat FF JJ;
%   end
   for ii=1:n;   
	   nor =[];
      Jact=[];   
      comb=1;
	   [vect1,val1]=eig(Fini); 
	   [vect2,val2]=eig(Fini');
	   for in1=1:6,
	      j0= indexs(in1,1);
	      j1= indexs(in1,2);
	      j2= indexs(in1,3);
      
	      if imag(sum(val1(:,1)))==0,
	         vec1=vect1(:,1)./vect1(j0,1);
	      else
	         vec1=vect1(:,3)./vect1(j0,3);
	      end
			x1=vec1(j1);
			y1=vec1(j2);
    
	      for in2=1:6,
		      i0= indexs(in2,1);
	   	   i1= indexs(in2,2);
	      	i2= indexs(in2,3);
	         
		      if imag(sum(val2(:,1)))==0,
	   	      vec2=vect2(:,1)./vect2(i0,1);
	      	else
	         	vec2=vect2(:,3)./vect2(i0,3);
            end           
				x2=vec2(i1);
				y2=vec2(i2);
	         
				a=Fini(i1,j1);
				b=Fini(i1,j2);
				c=Fini(i2,j1);
	         d=Fini(i2,j2);
	         
            switch comb
            case 1,
               Jact=	[ a*x2+c*y2,         0,         0,        -a,         0,         0,        -c,         0,         0;
							  b*x2+d*y2,         0,         0,        -b,         0,         0,        -d,         0,         0;
							  a*x1+b*y1,        -a,        -b,         0,         0,         0,         0,         0,         0;
							  c*x1+d*y1,        -c,        -d,         0,         0,         0,         0,         0,         0;
							      x1*x2,       -x2,         0,       -x1,         1,         0,         0,         0,         0;
							      y1*x2,         0,       -x2,       -y1,         0,         1,         0,         0,         0;
							      x1*y2,       -y2,         0,         0,         0,         0,       -x1,         1,         0;
							      y1*y2,         0,       -y2,         0,         0,         0,       -y1,         0,         1];
            case 2,
               Jact= [ a*x2+c*y2,         0,         0,        -c,         0,         0,        -a,         0,         0;
							  b*x2+d*y2,         0,         0,        -d,         0,         0,        -b,         0,         0;
							  a*x1+b*y1,        -a,        -b,         0,         0,         0,         0,         0,         0;
							  c*x1+d*y1,        -c,        -d,         0,         0,         0,         0,         0,         0;
							      x1*x2,       -x2,         0,         0,         0,         0,       -x1,         1,         0;
							      y1*x2,         0,       -x2,         0,         0,         0,       -y1,         0,         1;
							      x1*y2,       -y2,         0,       -x1,         1,         0,         0,         0,         0;
							      y1*y2,         0,       -y2,       -y1,         0,         1,         0,         0,         0];
            case 3,
               Jact= [        -a,         0,         0, a*x2+c*y2,         0,         0,        -c,         0,         0;
							         -b,         0,         0, b*x2+d*y2,         0,         0,        -d,         0,         0;
							          0,         0,         0, a*x1+b*y1,        -a,        -b,         0,         0,         0;
							          0,         0,         0, c*x1+d*y1,        -c,        -d,         0,         0,         0;
							        -x1,         1,         0,     x1*x2,       -x2,         0,         0,         0,         0;
							        -y1,         0,         1,     y1*x2,         0,       -x2,         0,         0,         0;
							          0,         0,         0,     x1*y2,       -y2,         0,       -x1,         1,         0;
							          0,         0,         0,     y1*y2,         0,       -y2,       -y1,         0,         1];
            case 4,
               Jact= [        -c,         0,         0, a*x2+c*y2,         0,         0,        -a,         0,         0;
							         -d,         0,         0, b*x2+d*y2,         0,         0,        -b,         0,         0;
							          0,         0,         0, a*x1+b*y1,        -a,        -b,         0,         0,         0;
							          0,         0,         0, c*x1+d*y1,        -c,        -d,         0,         0,         0;
							          0,         0,         0,     x1*x2,       -x2,         0,       -x1,         1,         0;
							          0,         0,         0,     y1*x2,         0,       -x2,       -y1,         0,         1;
							        -x1,         1,         0,     x1*y2,       -y2,         0,         0,         0,         0;
							        -y1,         0,         1,     y1*y2,         0,       -y2,         0,         0,         0];
            case 5,
               Jact= [        -a,         0,         0,        -c,         0,         0, a*x2+c*y2,         0,         0;
							         -b,         0,         0,        -d,         0,         0, b*x2+d*y2,         0,         0;
							          0,         0,         0,         0,         0,         0, a*x1+b*y1,        -a,        -b;
							          0,         0,         0,         0,         0,         0, c*x1+d*y1,        -c,        -d;
							        -x1,         1,         0,         0,         0,         0,     x1*x2,       -x2,         0;
							        -y1,         0,         1,         0,         0,         0,     y1*x2,         0,       -x2;
							          0,         0,         0,       -x1,         1,         0,     x1*y2,       -y2,         0;
							          0,         0,         0,       -y1,         0,         1,     y1*y2,         0,       -y2];
				case 6,
               Jact= [        -c,         0,         0,        -a,         0,         0, a*x2+c*y2,         0,         0;
							         -d,         0,         0,        -b,         0,         0, b*x2+d*y2,         0,         0;
							          0,         0,         0,         0,         0,         0, a*x1+b*y1,        -a,        -b;
							          0,         0,         0,         0,         0,         0, c*x1+d*y1,        -c,        -d;
							          0,         0,         0,       -x1,         1,         0,     x1*x2,       -x2,         0;
							          0,         0,         0,       -y1,         0,         1,     y1*x2,         0,       -x2;
							        -x1,         1,         0,         0,         0,         0,     x1*y2,       -y2,         0;
							        -y1,         0,         1,         0,         0,         0,     y1*y2,         0,       -y2];
            case 7,
               Jact= [ a*x2+c*y2,         0,         0,        -a,         0,         0,        -c,         0,         0;
							  b*x2+d*y2,         0,         0,        -b,         0,         0,        -d,         0,         0;
							  a*x1+b*y1,        -b,        -a,         0,         0,         0,         0,         0,         0;
							  c*x1+d*y1,        -d,        -c,         0,         0,         0,         0,         0,         0;
							      x1*x2,         0,       -x2,       -x1,         0,         1,         0,         0,         0;
							      y1*x2,       -x2,         0,       -y1,         1,         0,         0,         0,         0;
							      x1*y2,         0,       -y2,         0,         0,         0,       -x1,         0,         1;
							      y1*y2,       -y2,         0,         0,         0,         0,       -y1,         1,         0];
            case 8,
               Jact= [ a*x2+c*y2,         0,         0,        -c,         0,         0,        -a,         0,         0;
							  b*x2+d*y2,         0,         0,        -d,         0,         0,        -b,         0,         0;
							  a*x1+b*y1,        -b,        -a,         0,         0,         0,         0,         0,         0;
							  c*x1+d*y1,        -d,        -c,         0,         0,         0,         0,         0,         0;
							      x1*x2,         0,       -x2,         0,         0,         0,       -x1,         0,         1;
							      y1*x2,       -x2,         0,         0,         0,         0,       -y1,         1,         0;
							      x1*y2,         0,       -y2,       -x1,         0,         1,         0,         0,         0;
							      y1*y2,       -y2,         0,       -y1,         1,         0,         0,         0,         0];
            case 9,
               Jact= [        -a,         0,         0, a*x2+c*y2,         0,         0,        -c,         0,         0;
							         -b,         0,         0, b*x2+d*y2,         0,         0,        -d,         0,         0;
							          0,         0,         0, a*x1+b*y1,        -b,        -a,         0,         0,         0;
							          0,         0,         0, c*x1+d*y1,        -d,        -c,         0,         0,         0;
							        -x1,         0,         1,     x1*x2,         0,       -x2,         0,         0,         0;
							        -y1,         1,         0,     y1*x2,       -x2,         0,         0,         0,         0;
							          0,         0,         0,     x1*y2,         0,       -y2,       -x1,         0,         1;
							          0,         0,         0,     y1*y2,       -y2,         0,       -y1,         1,         0];
            case 10,
               Jact= [        -c,         0,         0, a*x2+c*y2,         0,         0,        -a,         0,         0;
							         -d,         0,         0, b*x2+d*y2,         0,         0,        -b,         0,         0;
							          0,         0,         0, a*x1+b*y1,        -b,        -a,         0,         0,         0;
							          0,         0,         0, c*x1+d*y1,        -d,        -c,         0,         0,         0;
							          0,         0,         0,     x1*x2,         0,       -x2,       -x1,         0,         1;
							          0,         0,         0,     y1*x2,       -x2,         0,       -y1,         1,         0;
							        -x1,         0,         1,     x1*y2,         0,       -y2,         0,         0,         0;
							        -y1,         1,         0,     y1*y2,       -y2,         0,         0,         0,         0];
            case 11,
               Jact= [        -a,         0,         0,        -c,         0,         0, a*x2+c*y2,         0,         0;
							         -b,         0,         0,        -d,         0,         0, b*x2+d*y2,         0,         0;
							          0,         0,         0,         0,         0,         0, a*x1+b*y1,        -b,        -a;
							          0,         0,         0,         0,         0,         0, c*x1+d*y1,        -d,        -c;
							        -x1,         0,         1,         0,         0,         0,     x1*x2,         0,       -x2;
							        -y1,         1,         0,         0,         0,         0,     y1*x2,       -x2,         0;
							          0,         0,         0,       -x1,         0,         1,     x1*y2,         0,       -y2;
							          0,         0,         0,       -y1,         1,         0,     y1*y2,       -y2,         0];
            case 12,
               Jact= [        -c,         0,         0,        -a,         0,         0, a*x2+c*y2,         0,         0;
							         -d,         0,         0,        -b,         0,         0, b*x2+d*y2,         0,         0;
							          0,         0,         0,         0,         0,         0, a*x1+b*y1,        -b,        -a;
							          0,         0,         0,         0,         0,         0, c*x1+d*y1,        -d,        -c;
							          0,         0,         0,       -x1,         0,         1,     x1*x2,         0,       -x2;
							          0,         0,         0,       -y1,         1,         0,     y1*x2,       -x2,         0;
							        -x1,         0,         1,         0,         0,         0,     x1*y2,         0,       -y2;
							        -y1,         1,         0,         0,         0,         0,     y1*y2,       -y2,         0];
            case 13,
               Jact= [         0, a*x2+c*y2,         0,         0,        -a,         0,         0,        -c,         0;
							          0, b*x2+d*y2,         0,         0,        -b,         0,         0,        -d,         0;
							         -a, a*x1+b*y1,        -b,         0,         0,         0,         0,         0,         0;
							         -c, c*x1+d*y1,        -d,         0,         0,         0,         0,         0,         0;
							        -x2,     x1*x2,         0,         1,       -x1,         0,         0,         0,         0;
							          0,     y1*x2,       -x2,         0,       -y1,         1,         0,         0,         0;
							        -y2,     x1*y2,         0,         0,         0,         0,         1,       -x1,         0;
							          0,     y1*y2,       -y2,         0,         0,         0,         0,       -y1,         1];
            case 14,
               Jact= [         0, a*x2+c*y2,         0,         0,        -c,         0,         0,        -a,         0;
							          0, b*x2+d*y2,         0,         0,        -d,         0,         0,        -b,         0;
							         -a, a*x1+b*y1,        -b,         0,         0,         0,         0,         0,         0;
							         -c, c*x1+d*y1,        -d,         0,         0,         0,         0,         0,         0;
							        -x2,     x1*x2,         0,         0,         0,         0,         1,       -x1,         0;
							          0,     y1*x2,       -x2,         0,         0,         0,         0,       -y1,         1;
							        -y2,     x1*y2,         0,         1,       -x1,         0,         0,         0,         0;
							          0,     y1*y2,       -y2,         0,       -y1,         1,         0,         0,         0];
            case 15,
               Jact= [         0,        -a,         0,         0, a*x2+c*y2,         0,         0,        -c,         0;
							          0,        -b,         0,         0, b*x2+d*y2,         0,         0,        -d,         0;
							          0,         0,         0,        -a, a*x1+b*y1,        -b,         0,         0,         0;
							          0,         0,         0,        -c, c*x1+d*y1,        -d,         0,         0,         0;
							          1,       -x1,         0,       -x2,     x1*x2,         0,         0,         0,         0;
							          0,       -y1,         1,         0,     y1*x2,       -x2,         0,         0,         0;
							          0,         0,         0,       -y2,     x1*y2,         0,         1,       -x1,         0;
							          0,         0,         0,         0,     y1*y2,       -y2,         0,       -y1,         1];
            case 16,
               Jact= [         0,        -c,         0,         0, a*x2+c*y2,         0,         0,        -a,         0;
							          0,        -d,         0,         0, b*x2+d*y2,         0,         0,        -b,         0;
							          0,         0,         0,        -a, a*x1+b*y1,        -b,         0,         0,         0;
							          0,         0,         0,        -c, c*x1+d*y1,        -d,         0,         0,         0;
							          0,         0,         0,       -x2,     x1*x2,         0,         1,       -x1,         0;
							          0,         0,         0,         0,     y1*x2,       -x2,         0,       -y1,         1;
							          1,       -x1,         0,       -y2,     x1*y2,         0,         0,         0,         0;
							          0,       -y1,         1,         0,     y1*y2,       -y2,         0,         0,         0];
            case 17,
               Jact= [         0,        -a,         0,         0,        -c,         0,         0, a*x2+c*y2,         0;
							          0,        -b,         0,         0,        -d,         0,         0, b*x2+d*y2,         0;
							          0,         0,         0,         0,         0,         0,        -a, a*x1+b*y1,        -b;
							          0,         0,         0,         0,         0,         0,        -c, c*x1+d*y1,        -d;
							          1,       -x1,         0,         0,         0,         0,       -x2,     x1*x2,         0;
							          0,       -y1,         1,         0,         0,         0,         0,     y1*x2,       -x2;
							          0,         0,         0,         1,       -x1,         0,       -y2,     x1*y2,         0;
							          0,         0,         0,         0,       -y1,         1,         0,     y1*y2,       -y2];
            case 18,
               Jact= [         0,        -c,         0,         0,        -a,         0,         0, a*x2+c*y2,         0;
							          0,        -d,         0,         0,        -b,         0,         0, b*x2+d*y2,         0;
							          0,         0,         0,         0,         0,         0,        -a, a*x1+b*y1,        -b;
							          0,         0,         0,         0,         0,         0,        -c, c*x1+d*y1,        -d;
							          0,         0,         0,         1,       -x1,         0,       -x2,     x1*x2,         0;
							          0,         0,         0,         0,       -y1,         1,         0,     y1*x2,       -x2;
							          1,       -x1,         0,         0,         0,         0,       -y2,     x1*y2,         0;
							          0,       -y1,         1,         0,         0,         0,         0,     y1*y2,       -y2];
            case 19,
               Jact= [         0, a*x2+c*y2,         0,         0,        -a,         0,         0,        -c,         0;
							          0, b*x2+d*y2,         0,         0,        -b,         0,         0,        -d,         0;
							         -b, a*x1+b*y1,        -a,         0,         0,         0,         0,         0,         0;
							         -d, c*x1+d*y1,        -c,         0,         0,         0,         0,         0,         0;
							          0,     x1*x2,       -x2,         0,       -x1,         1,         0,         0,         0;
							        -x2,     y1*x2,         0,         1,       -y1,         0,         0,         0,         0;
							          0,     x1*y2,       -y2,         0,         0,         0,         0,       -x1,         1;
							        -y2,     y1*y2,         0,         0,         0,         0,         1,       -y1,         0];
            case 20,
               Jact= [         0, a*x2+c*y2,         0,         0,        -c,         0,         0,        -a,         0;
							          0, b*x2+d*y2,         0,         0,        -d,         0,         0,        -b,         0;
							         -b, a*x1+b*y1,        -a,         0,         0,         0,         0,         0,         0;
							         -d, c*x1+d*y1,        -c,         0,         0,         0,         0,         0,         0;
							          0,     x1*x2,       -x2,         0,         0,         0,         0,       -x1,         1;
							        -x2,     y1*x2,         0,         0,         0,         0,         1,       -y1,         0;
							          0,     x1*y2,       -y2,         0,       -x1,         1,         0,         0,         0;
							        -y2,     y1*y2,         0,         1,       -y1,         0,         0,         0,         0];
            case 21,
               Jact= [         0,        -a,         0,         0, a*x2+c*y2,         0,         0,        -c,         0;
							          0,        -b,         0,         0, b*x2+d*y2,         0,         0,        -d,         0;
							          0,         0,         0,        -b, a*x1+b*y1,        -a,         0,         0,         0;
							          0,         0,         0,        -d, c*x1+d*y1,        -c,         0,         0,         0;
							          0,       -x1,         1,         0,     x1*x2,       -x2,         0,         0,         0;
							          1,       -y1,         0,       -x2,     y1*x2,         0,         0,         0,         0;
							          0,         0,         0,         0,     x1*y2,       -y2,         0,       -x1,         1;
							          0,         0,         0,       -y2,     y1*y2,         0,         1,       -y1,         0];
            case 22,
               Jact= [         0,        -c,         0,         0, a*x2+c*y2,         0,         0,        -a,         0;
							          0,        -d,         0,         0, b*x2+d*y2,         0,         0,        -b,         0;
							          0,         0,         0,        -b, a*x1+b*y1,        -a,         0,         0,         0;
							          0,         0,         0,        -d, c*x1+d*y1,        -c,         0,         0,         0;
							          0,         0,         0,         0,     x1*x2,       -x2,         0,       -x1,         1;
							          0,         0,         0,       -x2,     y1*x2,         0,         1,       -y1,         0;
							          0,       -x1,         1,         0,     x1*y2,       -y2,         0,         0,         0;
							          1,       -y1,         0,       -y2,     y1*y2,         0,         0,         0,         0];
            case 23,
               Jact= [         0,        -a,         0,         0,        -c,         0,         0, a*x2+c*y2,         0;
							          0,        -b,         0,         0,        -d,         0,         0, b*x2+d*y2,         0;
							          0,         0,         0,         0,         0,         0,        -b, a*x1+b*y1,        -a;
							          0,         0,         0,         0,         0,         0,        -d, c*x1+d*y1,        -c;
							          0,       -x1,         1,         0,         0,         0,         0,     x1*x2,       -x2;
							          1,       -y1,         0,         0,         0,         0,       -x2,     y1*x2,         0;
							          0,         0,         0,         0,       -x1,         1,         0,     x1*y2,       -y2;
							          0,         0,         0,         1,       -y1,         0,       -y2,     y1*y2,         0];
            case 24,
               Jact= [         0,        -c,         0,         0,        -a,         0,         0, a*x2+c*y2,         0;
							          0,        -d,         0,         0,        -b,         0,         0, b*x2+d*y2,         0;
							          0,         0,         0,         0,         0,         0,        -b, a*x1+b*y1,        -a;
							          0,         0,         0,         0,         0,         0,        -d, c*x1+d*y1,        -c;
							          0,         0,         0,         0,       -x1,         1,         0,     x1*x2,       -x2;
							          0,         0,         0,         1,       -y1,         0,       -x2,     y1*x2,         0;
							          0,       -x1,         1,         0,         0,         0,         0,     x1*y2,       -y2;
							          1,       -y1,         0,         0,         0,         0,       -y2,     y1*y2,         0];
            case 25,
               Jact= [         0,         0, a*x2+c*y2,         0,         0,        -a,         0,         0,        -c;
							          0,         0, b*x2+d*y2,         0,         0,        -b,         0,         0,        -d;
							         -a,        -b, a*x1+b*y1,         0,         0,         0,         0,         0,         0;
							         -c,        -d, c*x1+d*y1,         0,         0,         0,         0,         0,         0;
							        -x2,         0,     x1*x2,         1,         0,       -x1,         0,         0,         0;
							          0,       -x2,     y1*x2,         0,         1,       -y1,         0,         0,         0;
							        -y2,         0,     x1*y2,         0,         0,         0,         1,         0,       -x1;
							          0,       -y2,     y1*y2,         0,         0,         0,         0,         1,       -y1];
            case 26,
               Jact= [         0,         0, a*x2+c*y2,         0,         0,        -c,         0,         0,        -a;
							          0,         0, b*x2+d*y2,         0,         0,        -d,         0,         0,        -b;
							         -a,        -b, a*x1+b*y1,         0,         0,         0,         0,         0,         0;
							         -c,        -d, c*x1+d*y1,         0,         0,         0,         0,         0,         0;
							        -x2,         0,     x1*x2,         0,         0,         0,         1,         0,       -x1;
							          0,       -x2,     y1*x2,         0,         0,         0,         0,         1,       -y1;
							        -y2,         0,     x1*y2,         1,         0,       -x1,         0,         0,         0;
							          0,       -y2,     y1*y2,         0,         1,       -y1,         0,         0,         0];
            case 27,
               Jact= [         0,         0,        -a,         0,         0, a*x2+c*y2,         0,         0,        -c;
							          0,         0,        -b,         0,         0, b*x2+d*y2,         0,         0,        -d;
							          0,         0,         0,        -a,        -b, a*x1+b*y1,         0,         0,         0;
							          0,         0,         0,        -c,        -d, c*x1+d*y1,         0,         0,         0;
							          1,         0,       -x1,       -x2,         0,     x1*x2,         0,         0,         0;
							          0,         1,       -y1,         0,       -x2,     y1*x2,         0,         0,         0;
							          0,         0,         0,       -y2,         0,     x1*y2,         1,         0,       -x1;
							          0,         0,         0,         0,       -y2,     y1*y2,         0,         1,       -y1];
            case 28,
               Jact= [         0,         0,        -c,         0,         0, a*x2+c*y2,         0,         0,        -a;
							          0,         0,        -d,         0,         0, b*x2+d*y2,         0,         0,        -b;
							          0,         0,         0,        -a,        -b, a*x1+b*y1,         0,         0,         0;
							          0,         0,         0,        -c,        -d, c*x1+d*y1,         0,         0,         0;
							          0,         0,         0,       -x2,         0,     x1*x2,         1,         0,       -x1;
							          0,         0,         0,         0,       -x2,     y1*x2,         0,         1,       -y1;
							          1,         0,       -x1,       -y2,         0,     x1*y2,         0,         0,         0;
							          0,         1,       -y1,         0,       -y2,     y1*y2,         0,         0,         0];
            case 29,
               Jact= [         0,         0,        -a,         0,         0,        -c,         0,         0, a*x2+c*y2;
							          0,         0,        -b,         0,         0,        -d,         0,         0, b*x2+d*y2;
							          0,         0,         0,         0,         0,         0,        -a,        -b, a*x1+b*y1;
							          0,         0,         0,         0,         0,         0,        -c,        -d, c*x1+d*y1;
							          1,         0,       -x1,         0,         0,         0,       -x2,         0,     x1*x2;
							          0,         1,       -y1,         0,         0,         0,         0,       -x2,     y1*x2;
							          0,         0,         0,         1,         0,       -x1,       -y2,         0,     x1*y2;
							          0,         0,         0,         0,         1,       -y1,         0,       -y2,     y1*y2];
            case 30,
               Jact= [         0,         0,        -c,         0,         0,        -a,         0,         0, a*x2+c*y2;
							          0,         0,        -d,         0,         0,        -b,         0,         0, b*x2+d*y2;
							          0,         0,         0,         0,         0,         0,        -a,        -b, a*x1+b*y1;
							          0,         0,         0,         0,         0,         0,        -c,        -d, c*x1+d*y1;
							          0,         0,         0,         1,         0,       -x1,       -x2,         0,     x1*x2;
							          0,         0,         0,         0,         1,       -y1,         0,       -x2,     y1*x2;
							          1,         0,       -x1,         0,         0,         0,       -y2,         0,     x1*y2;
							          0,         1,       -y1,         0,         0,         0,         0,       -y2,     y1*y2];
            case 31,
               Jact= [         0,         0, a*x2+c*y2,         0,         0,        -a,         0,         0,        -c;
							          0,         0, b*x2+d*y2,         0,         0,        -b,         0,         0,        -d;
							         -b,        -a, a*x1+b*y1,         0,         0,         0,         0,         0,         0;
							         -d,        -c, c*x1+d*y1,         0,         0,         0,         0,         0,         0;
							          0,       -x2,     x1*x2,         0,         1,       -x1,         0,         0,         0;
							        -x2,         0,     y1*x2,         1,         0,       -y1,         0,         0,         0;
							          0,       -y2,     x1*y2,         0,         0,         0,         0,         1,       -x1;
							        -y2,         0,     y1*y2,         0,         0,         0,         1,         0,       -y1];
            case 32,
               Jact= [         0,         0, a*x2+c*y2,         0,         0,        -c,         0,         0,        -a;
							          0,         0, b*x2+d*y2,         0,         0,        -d,         0,         0,        -b;
							         -b,        -a, a*x1+b*y1,         0,         0,         0,         0,         0,         0;
							         -d,        -c, c*x1+d*y1,         0,         0,         0,         0,         0,         0;
							          0,       -x2,     x1*x2,         0,         0,         0,         0,         1,       -x1;
							        -x2,         0,     y1*x2,         0,         0,         0,         1,         0,       -y1;
							          0,       -y2,     x1*y2,         0,         1,       -x1,         0,         0,         0;
							        -y2,         0,     y1*y2,         1,         0,       -y1,         0,         0,         0];
            case 33,
               Jact= [         0,         0,        -a,         0,         0, a*x2+c*y2,         0,         0,        -c;
							          0,         0,        -b,         0,         0, b*x2+d*y2,         0,         0,        -d;
							          0,         0,         0,        -b,        -a, a*x1+b*y1,         0,         0,         0;
							          0,         0,         0,        -d,        -c, c*x1+d*y1,         0,         0,         0;
							          0,         1,       -x1,         0,       -x2,     x1*x2,         0,         0,         0;
							          1,         0,       -y1,       -x2,         0,     y1*x2,         0,         0,         0;
							          0,         0,         0,         0,       -y2,     x1*y2,         0,         1,       -x1;
							          0,         0,         0,       -y2,         0,     y1*y2,         1,         0,       -y1];
            case 34,
               Jact= [         0,         0,        -c,         0,         0, a*x2+c*y2,         0,         0,        -a;
							          0,         0,        -d,         0,         0, b*x2+d*y2,         0,         0,        -b;
							          0,         0,         0,        -b,        -a, a*x1+b*y1,         0,         0,         0;
							          0,         0,         0,        -d,        -c, c*x1+d*y1,         0,         0,         0;
							          0,         0,         0,         0,       -x2,     x1*x2,         0,         1,       -x1;
							          0,         0,         0,       -x2,         0,     y1*x2,         1,         0,       -y1;
							          0,         1,       -x1,         0,       -y2,     x1*y2,         0,         0,         0;
							          1,         0,       -y1,       -y2,         0,     y1*y2,         0,         0,         0];
            case 35,
               Jact= [         0,         0,        -a,         0,         0,        -c,         0,         0, a*x2+c*y2;
							          0,         0,        -b,         0,         0,        -d,         0,         0, b*x2+d*y2;
							          0,         0,         0,         0,         0,         0,        -b,        -a, a*x1+b*y1;
							          0,         0,         0,         0,         0,         0,        -d,        -c, c*x1+d*y1;
							          0,         1,       -x1,         0,         0,         0,         0,       -x2,     x1*x2;
							          1,         0,       -y1,         0,         0,         0,       -x2,         0,     y1*x2;
							          0,         0,         0,         0,         1,       -x1,         0,       -y2,     x1*y2;
							          0,         0,         0,         1,         0,       -y1,       -y2,         0,     y1*y2];
            case 36,
               Jact= [         0,         0,        -c,         0,         0,        -a,         0,         0, a*x2+c*y2;
							          0,         0,        -d,         0,         0,        -b,         0,         0, b*x2+d*y2;
							          0,         0,         0,         0,         0,         0,        -b,        -a, a*x1+b*y1;
							          0,         0,         0,         0,         0,         0,        -d,        -c, c*x1+d*y1;
							          0,         0,         0,         0,         1,       -x1,         0,       -x2,     x1*x2;
							          0,         0,         0,         1,         0,       -y1,       -x2,         0,     y1*x2;
							          0,         1,       -x1,         0,         0,         0,         0,       -y2,     x1*y2;
							          1,         0,       -y1,         0,         0,         0,       -y2,         0,     y1*y2];
            case 37,
					disp('Error')
            end            
            vectJ=0;
	         for jc=1:9,
	            vectJ(jc,comb)=det([Jact(:,1:(jc-1)) Jact(:,(jc+1):9)]);
	         end
				nor=[nor norm(vectJ(:,comb))];
            
%            nor2=[nor2 (a*d-b*c)^2*sqrt(x1^2+y1^2+1)*sqrt(x2^2+y2^2+1)];

	         comb=comb+1;
	      end;
	   end;
	   [vm,im]=max(nor);
      if vmini<vm,
			in1=fix((im-1)/6)+1;
		   in2=mod(im,6);
		   if in2==0,
		      in2=6;
		   end
         
		   j0= indexs(in1,1);
		   j1= indexs(in1,2);
		   j2= indexs(in1,3);
		   i0= indexs(in2,1);
		   i1= indexs(in2,2);
		  	i2= indexs(in2,3);
	      
		   if imag(sum(val1(:,1)))==0,
		      vec1=vect1(:,1)./vect1(j0,1);
		   else
		      vec1=vect1(:,3)./vect1(j0,3);
		   end
	      
	   	x1=vec1(j1);
		  	y1=vec1(j2);
		         
			if imag(sum(val2(:,1)))==0,
		   	vec2=vect2(:,1)./vect2(i0,1);
		   else
		      vec2=vect2(:,3)./vect2(i0,3);
		   end   
		     
		  	x2=vec2(i1);
		   y2=vec2(i2);
		         
		   a=Fini(i1,j1);
		   b=Fini(i1,j2);
		   c=Fini(i2,j1);
		   d=Fini(i2,j2);
		
			switch im,
         case 1,
            F= [ (a*x1+b*y1)*x2+(c*x1+d*y1)*y2,                    -a*x2-c*y2,                    -b*x2-d*y2;
					                     -a*x1-b*y1,                             a,                             b;
					                     -c*x1-d*y1,                             c,                             d];
         case 2,
            F= [ (a*x1+b*y1)*x2+(c*x1+d*y1)*y2,                    -a*x2-c*y2,                    -b*x2-d*y2;
					                     -c*x1-d*y1,                             c,                             d;
					                     -a*x1-b*y1,                             a,                             b];
         case 3,
            F= [                    -a*x1-b*y1,                             a,                             b;
					  (a*x1+b*y1)*x2+(c*x1+d*y1)*y2,                    -a*x2-c*y2,                    -b*x2-d*y2;
					                     -c*x1-d*y1,                             c,                             d];
         case 4,
            F= [                    -c*x1-d*y1,                             c,                             d;
					  (a*x1+b*y1)*x2+(c*x1+d*y1)*y2,                    -a*x2-c*y2,                    -b*x2-d*y2;
					                     -a*x1-b*y1,                             a,                             b];
         case 5,
            F= [                    -a*x1-b*y1,                             a,                             b;
					                     -c*x1-d*y1,                             c,                             d;
					  (a*x1+b*y1)*x2+(c*x1+d*y1)*y2,                    -a*x2-c*y2,                    -b*x2-d*y2];
         case 6,
            F= [                    -c*x1-d*y1,                             c,                             d;
					                     -a*x1-b*y1,                             a,                             b;
					  (a*x1+b*y1)*x2+(c*x1+d*y1)*y2,                    -a*x2-c*y2,                    -b*x2-d*y2];
         case 7,
            F= [ (a*x1+b*y1)*x2+(c*x1+d*y1)*y2,                    -b*x2-d*y2,                    -a*x2-c*y2;
					                     -a*x1-b*y1,                             b,                             a;
					                     -c*x1-d*y1,                             d,                             c];
         case 8,
            F= [ (a*x1+b*y1)*x2+(c*x1+d*y1)*y2,                    -b*x2-d*y2,                    -a*x2-c*y2;
					                     -c*x1-d*y1,                             d,                             c;
					                     -a*x1-b*y1,                             b,                             a];
         case 9,
            F= [                    -a*x1-b*y1,                             b,                             a;
					  (a*x1+b*y1)*x2+(c*x1+d*y1)*y2,                    -b*x2-d*y2,                    -a*x2-c*y2;
					                     -c*x1-d*y1,                             d,                             c];
         case 10,
            F= [                    -c*x1-d*y1,                             d,                             c;
					  (a*x1+b*y1)*x2+(c*x1+d*y1)*y2,                    -b*x2-d*y2,                    -a*x2-c*y2;
					                     -a*x1-b*y1,                             b,                             a];
         case 11,
            F= [                    -a*x1-b*y1,                             b,                             a;
					                     -c*x1-d*y1,                             d,                             c;
					  (a*x1+b*y1)*x2+(c*x1+d*y1)*y2,                    -b*x2-d*y2,                    -a*x2-c*y2];
         case 12,
            F= [                    -c*x1-d*y1,                             d,                             c;
					                     -a*x1-b*y1,                             b,                             a;
					  (a*x1+b*y1)*x2+(c*x1+d*y1)*y2,                    -b*x2-d*y2,                    -a*x2-c*y2];
         case 13,
            F= [                    -a*x2-c*y2, (a*x1+b*y1)*x2+(c*x1+d*y1)*y2,                    -b*x2-d*y2;
					                              a,                    -a*x1-b*y1,                             b;
					                              c,                    -c*x1-d*y1,                             d];
         case 14,
            F= [                    -a*x2-c*y2, (a*x1+b*y1)*x2+(c*x1+d*y1)*y2,                    -b*x2-d*y2;
					                              c,                    -c*x1-d*y1,                             d;
					                              a,                    -a*x1-b*y1,                             b];
         case 15,
            F= [                             a,                    -a*x1-b*y1,                             b;
					                     -a*x2-c*y2, (a*x1+b*y1)*x2+(c*x1+d*y1)*y2,                    -b*x2-d*y2;
					                              c,                    -c*x1-d*y1,                             d];
         case 16,
            F= [                             c,                    -c*x1-d*y1,                             d;
					                     -a*x2-c*y2, (a*x1+b*y1)*x2+(c*x1+d*y1)*y2,                    -b*x2-d*y2;
					                              a,                    -a*x1-b*y1,                             b];
         case 17,
            F= [                             a,                    -a*x1-b*y1,                             b;
					                              c,                    -c*x1-d*y1,                             d;
					                     -a*x2-c*y2, (a*x1+b*y1)*x2+(c*x1+d*y1)*y2,                    -b*x2-d*y2];
         case 18,
            F= [                             c,                    -c*x1-d*y1,                             d;
					                              a,                    -a*x1-b*y1,                             b;
					                     -a*x2-c*y2, (a*x1+b*y1)*x2+(c*x1+d*y1)*y2,                    -b*x2-d*y2];
         case 19,
            F= [                    -b*x2-d*y2, (a*x1+b*y1)*x2+(c*x1+d*y1)*y2,                    -a*x2-c*y2;
					                              b,                    -a*x1-b*y1,                             a;
					                              d,                    -c*x1-d*y1,                             c];
         case 20,
            F= [                    -b*x2-d*y2, (a*x1+b*y1)*x2+(c*x1+d*y1)*y2,                    -a*x2-c*y2;
					                              d,                    -c*x1-d*y1,                             c;
					                              b,                    -a*x1-b*y1,                             a];
         case 21,
            F= [                             b,                    -a*x1-b*y1,                             a;
					                     -b*x2-d*y2, (a*x1+b*y1)*x2+(c*x1+d*y1)*y2,                    -a*x2-c*y2;
					                              d,                    -c*x1-d*y1,                             c];
         case 22,
            F= [                             d,                    -c*x1-d*y1,                             c;
					                     -b*x2-d*y2, (a*x1+b*y1)*x2+(c*x1+d*y1)*y2,                    -a*x2-c*y2;
					                              b,                    -a*x1-b*y1,                             a];
         case 23,
            F= [                             b,                    -a*x1-b*y1,                             a;
					                              d,                    -c*x1-d*y1,                             c;
					                     -b*x2-d*y2, (a*x1+b*y1)*x2+(c*x1+d*y1)*y2,                    -a*x2-c*y2];
         case 24,
            F= [                             d,                    -c*x1-d*y1,                             c;
					                              b,                    -a*x1-b*y1,                             a;
					                     -b*x2-d*y2, (a*x1+b*y1)*x2+(c*x1+d*y1)*y2,                    -a*x2-c*y2];
         case 25,
            F= [                    -a*x2-c*y2,                    -b*x2-d*y2, (a*x1+b*y1)*x2+(c*x1+d*y1)*y2;
					                              a,                             b,                    -a*x1-b*y1;
					                              c,                             d,                    -c*x1-d*y1];
         case 26,
            F= [                    -a*x2-c*y2,                    -b*x2-d*y2, (a*x1+b*y1)*x2+(c*x1+d*y1)*y2;
					                              c,                             d,                    -c*x1-d*y1;
					                              a,                             b,                    -a*x1-b*y1];
         case 27,
            F= [                             a,                             b,                    -a*x1-b*y1;
					                     -a*x2-c*y2,                    -b*x2-d*y2, (a*x1+b*y1)*x2+(c*x1+d*y1)*y2;
					                              c,                             d,                    -c*x1-d*y1];
         case 28,
            F= [                             c,                             d,                    -c*x1-d*y1;
					                     -a*x2-c*y2,                    -b*x2-d*y2, (a*x1+b*y1)*x2+(c*x1+d*y1)*y2;
					                              a,                             b,                    -a*x1-b*y1];
         case 29,
            F= [                             a,                             b,                    -a*x1-b*y1;
					                              c,                             d,                    -c*x1-d*y1;
					                     -a*x2-c*y2,                    -b*x2-d*y2, (a*x1+b*y1)*x2+(c*x1+d*y1)*y2];
         case 30,
            F= [                             c,                             d,                    -c*x1-d*y1;
					                              a,                             b,                    -a*x1-b*y1;
					                     -a*x2-c*y2,                    -b*x2-d*y2, (a*x1+b*y1)*x2+(c*x1+d*y1)*y2];
         case 31,
            F= [                    -b*x2-d*y2,                    -a*x2-c*y2, (a*x1+b*y1)*x2+(c*x1+d*y1)*y2;
					                              b,                             a,                    -a*x1-b*y1;
					                              d,                             c,                    -c*x1-d*y1];
         case 32,
            F= [                    -b*x2-d*y2,                    -a*x2-c*y2, (a*x1+b*y1)*x2+(c*x1+d*y1)*y2;
					                              d,                             c,                    -c*x1-d*y1;
					                              b,                             a,                    -a*x1-b*y1];
         case 33,
            F= [                             b,                             a,                    -a*x1-b*y1;
					                     -b*x2-d*y2,                    -a*x2-c*y2, (a*x1+b*y1)*x2+(c*x1+d*y1)*y2;
					                              d,                             c,                    -c*x1-d*y1];
         case 34,
            F= [                             d,                             c,                    -c*x1-d*y1;
					                     -b*x2-d*y2,                    -a*x2-c*y2, (a*x1+b*y1)*x2+(c*x1+d*y1)*y2;
					                              b,                             a,                    -a*x1-b*y1];
         case 35,
            F= [                             b,                             a,                    -a*x1-b*y1;
					                              d,                             c,                    -c*x1-d*y1;
					                     -b*x2-d*y2,                    -a*x2-c*y2, (a*x1+b*y1)*x2+(c*x1+d*y1)*y2];
         case 36,
            F= [                             d,                             c,                    -c*x1-d*y1;
					                              b,                             a,                    -a*x1-b*y1;
					                     -b*x2-d*y2,                    -a*x2-c*y2, (a*x1+b*y1)*x2+(c*x1+d*y1)*y2];
         case 37,
				disp('Error')
         end         
         Fini=F;    
         vmini=vm;
	   else
         F=Fini;   
         break;
	   end	      
	end;
end;





