%funmatLMedSMest    Matriu Fonemantal amb el metode Least Median of Squares amb valors propis
%    [F] = funmatLMedSeig(M,b,P,ep) 
%
%    M matriu de 4-per-7 amb les coordendes de 7 punts en dues imatges
%        primera fila: coordenada x de la primera imatge
%        segona fila:  coordenada y de la primera imatge
%        tercera fila: coordenada x de la segona imatge
%        quarta fila:  coordenada y de la segona imatge
%    b  buckets
%    P  probabilitat que hi hagui un F sense outliers
%    ep ratio d'ouliers
%
%    F matriu fonamental de 3-per-3 amb l'origen de coordenades del m�n a la segona camera
%
% by X. Armangue
% (c) Mr3D - University of Girona, September 2002
%
function [F,wout]=funmatLMedSMest(M,b,P,ep,n,error)

if (size(M,1)~=4) | (size(M,2)<8),
   disp('Error: parametres incorrectes')
else
   maxx1=max(M(1,:));
   minx1=min(M(1,:));
   bx1=maxx1-minx1;
   maxy1=max(M(2,:));
   miny1=min(M(2,:));
   by1=maxy1-miny1;
   
   bucket=zeros(2,b*b);
   pointb=[];
   
   for i=1:size(M,2),
      bucketx=fix((M(1,i)-minx1)*b/bx1);
      if bucketx==b,
         bucketx=b-1;
      end
      buckety=fix((M(2,i)-miny1)*b/by1);
      if buckety==b,
         buckety=b-1;
      end      
      pointb(i)=b*buckety+bucketx+1;
      bucket(1,pointb(i))=bucket(1,pointb(i))+1;
      bucket(2,pointb(i))=0;
   end
   
   
   mat=round(log(1-P)/log(1-(1-ep)^7));
   for inf=1:mat,
		for i=1:b*b,
      	bucket(2,i)=0;
      end;
	   Mtemp=[];
	   for i=1:7,
	      bo=0;
	      while bo==0,
	         inpu=fix(rand(1,1)*size(M,2)+1);
	         if inpu==size(M,2)+1,
	            inpu=size(M,2);
	         end         
	            
            if (bucket(2,pointb(inpu))==0) & (bucket(1,pointb(inpu))>0),
               bucket(2,pointb(inpu))=1;
               Mtemp=[Mtemp M(:,inpu)];
               bo=1;
            end            
   	   end
	   end            
      FF(:,:,inf)=funmat7p(Mtemp);
      
      aux=[];
      r=[];
	   for i=1:size(M,2),
	      m1 = [M(1,i) ; M(2,i); 1];
	      m2 = [M(3,i) ; M(4,i); 1];
         
			l=FF(:,:,inf)*m2;
	      d=m1'*l/sqrt(l(1)^2+l(2)^2);
   	   dmod1=abs(d);
         
			l=FF(:,:,inf)'*m1;
	      d=m2'*l/sqrt(l(1)^2+l(2)^2);
   	   dmod2=abs(d);    
         
         r=[r dmod1^2+dmod2^2];
	   end
      Mj(inf)=median(r);
   end
   [Jm,J]=min(Mj);
   Fini=FF(:,:,J);
   
   
   switch 1, % diat�ncies (1) o residus (2)
   case 1
      dmod=[];
		for i=1:size(M,2),
	      m1 = [M(1,i) ; M(2,i); 1];
	      m2 = [M(3,i) ; M(4,i); 1];
	      
	      l=Fini*m2;
	      d=m1'*l/sqrt(l(1)^2+l(2)^2);
	      dmod1=abs(d);
	      
	      l=Fini'*m1;
	      d=m2'*l/sqrt(l(1)^2+l(2)^2);
	      dmod2=abs(d);    
	      
	      dmod=[dmod dmod1^2+dmod2^2];      
	   end
      
      switch 2, % sigma real (1) o estimada (2)
      case 1,
         sig=std(sqrt(dmod));
      case 2,
         sig=1.4826*(1+5/(size(M,2)-7))*sqrt(Mj(J));
      end
      
      switch 2, % pesos (2)
      case 1,
         for i=1:size(M,2),
		   	if sqrt(dmod(i))<=3*sig,
		         w(i)=1;
		      else
		         w(i)=0;
		      end
         end      
      case 2,
         for i=1:size(M,2),
				if dmod(i)^2<=((2.5*sig)^2)*1,
			       w(i)=1;
		      else
			       w(i)=0;
		      end
         end
      end 
      
   case 2,
      r=[];
	   for i=1:size(M,2),
			m1 = [M(1,i) ; M(2,i); 1];
	      m2 = [M(3,i) ; M(4,i); 1];
	      r=[r m1'*Fini*m2];
	   end
	   sig=std(r);

   	for i=1:size(M,2),
			if abs(r(i))<=(3*sig),
	      	w(i)=1;
	      else
	         w(i)=0;
	      end      
      end     
   end
   
   
%   U=[];
%   UU=[];
%	for i=1:size(M,2),
%      if w(i)~=0,
%      	U=[M(1,i)*M(3,i) M(1,i)*M(4,i) M(1,i) M(2,i)*M(3,i) M(2,i)*M(4,i) M(2,i) M(3,i) M(4,i) 1];      
%         UU=[UU ; U];         
%      end      
%  	end
	
%	[V,D]=eig(UU'*UU);    
%   [minim,ind]=min(sum(D));
%   f=V(:,ind);	
%   F=[f(1) f(2) f(3); f(4) f(5) f(6); f(7) f(8) f(9)];





Mini=M;
M=goodpoints(M,w,0.5);
www=w;
   w=ones(1,size(M,2));
   ga=ones(1,size(M,2));

%		F=Fini;
%		for i=1:size(M,2),
%         r(i)=[M(1:2,i) ; 1]'*F*[M(3:4,i) ; 1];
%         x1=M(1,i);
%         y1=M(2,i);
%         x2=M(3,i);
%         y2=M(4,i);
%         rx2=F(1,1)*x1+F(2,1)*y1+F(3,1);
%         ry2=F(1,2)*x1+F(2,2)*y1+F(3,2);
%         rx1=F(1,1)*x2+F(1,2)*y2+F(1,3);
%         ry1=F(2,1)*x2+F(2,2)*y2+F(2,3);
%         w(i)=sqrt(1/(rx2^2+ry2^2+rx1^2+ry1^2));
%         d(i)=w(i)*r(i);
%      end;
%      
%      sig=median(abs(d))/0.6745;
%      
%      for i=1:size(M,2),
%         if abs(d(i))<sig,
%            ga(i)=1;
%         elseif abs(d(i))<3*sig,
%            ga(i)=sig/abs(d(i));
%         else
%            ga(i)=0;
%         end
%      end



	[mitjanaini,stddev,minim,maxim]=funmatError(goodpoints(M,w,0.75),Fini);
   Fant=Fini;
%   Fini=[0 0 0; 0 0 0; 0 0 0];
	U=[];
	for i=1:size(M,2),
  	   U=[U; M(1,i)*M(3,i) M(1,i)*M(4,i) M(1,i) M(2,i)*M(3,i) M(2,i)*M(4,i) M(2,i) M(3,i) M(4,i) 1];
	end
     
   for iter=1:n,      
%      fummatplot(goodpoints(M,w,0.75),F)
%     	[mitjanaLMedSeig,stddevLMedSeig,minimLMedSeig,maximLMedSeig]=funmatError(goodpoints(M,w,0.75),F)
%      pause
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


	j=1;
	for i=1:size(www,2),
      if www(i)==1,
         wout(i)=ga(j);
         j=j+1;
      end      
   end
      
	[mitjanafi,stddev,minim,maxim]=funmatError(goodpoints(M,w,0.75),F);
   if mitjanaini<mitjanafi,
      [F,wout]=funmatLMedSMest(Mini,b,P,ep,n,error);
   end        
end
