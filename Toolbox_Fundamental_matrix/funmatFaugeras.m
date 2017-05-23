%funmatFaugeras   Matriu Fonemantal imposant que rank(F)=2, metode proposat per Faugeras (1995)
%    [F] = funmatFaugeras(M) 
%
%    M matriu de 4-per-n amb les coordendes de n punts en dues imatges
%      on n és mes gran o igual que 8
%        primera fila: coordenada x de la primera imatge
%        segona fila:  coordenada y de la primera imatge
%        tercera fila: coordenada x de la segona imatge
%        quarta fila:  coordenada y de la segona imatge
%
%    F matriu fonamental de 3-per-3 amb rank(F)=2 (det(F)=0)
%
% by X. Armangue
% (c) Mr3D - University of Girona, September 2002
%
function [F]=funmatFaugeras(M)

if (size(M,1)~=4) | (size(M,2)<8),
    disp('Error: parametres incorrectes')
else
    U=[];
    for i=1:2:14,
        U=[U; M(1,i)*M(3,i) M(1,i)*M(4,i) M(1,i) M(2,i)*M(3,i) M(2,i)*M(4,i) M(2,i) M(3,i) M(4,i) 1];
    end
    
    [UU,SS,VV]=svd(U);
    
    F1=[VV(1,8) VV(2,8) VV(3,8); VV(4,8) VV(5,8) VV(6,8); VV(7,8) VV(8,8) VV(9,8)];
    F2=[VV(1,9) VV(2,9) VV(3,9); VV(4,9) VV(5,9) VV(6,9); VV(7,9) VV(8,9) VV(9,9)];
    
    F111=F1(1,1);	F112=F1(1,2);	F113=F1(1,3);
    F121=F1(2,1);	F122=F1(2,2);	F123=F1(2,3);
    F131=F1(3,1);	F132=F1(3,2);	F133=F1(3,3);
    
    F211=F2(1,1);	F212=F2(1,2);	F213=F2(1,3);
    F221=F2(2,1);	F222=F2(2,2);	F223=F2(2,3);
    F231=F2(3,1);	F232=F2(3,2);	F233=F2(3,3);
    
    a=1;
    d=F221*F213*F232-F211*F223*F232+F211*F222*F233-F221*F212*F233+F231*F212*F223-F231*F213*F222;
    c=+F231*a*F112*F223+F231*F212*a*F123-3*F231*F212*F223*a-F231*a*F113*F222-F231*F213*a*F122+3*F231*F213*F222*a+a*F131*F212*F223-a*F131*F213*F222+a*F121*F213*F232-F221*a*F112*F233-F221*F212*a*F133+3*F221*F212*F233*a+F221*a*F113*F232+F221*F213*a*F132-3*F221*F213*F232*a+F211*a*F122*F233+F211*F222*a*F133-3*F211*F222*F233*a-F211*a*F123*F232-F211*F223*a*F132+3*F211*F223*F232*a-a*F121*F212*F233+a*F111*F222*F233-a*F111*F223*F232;
    b=+2*F231*F213*a^2*F122-3*F231*F213*a^2*F222+F231*a^2*F112*F123-2*F231*a^2*F112*F223-2*F231*F212*a^2*F123+3*F231*F212*a^2*F223-F231*a^2*F113*F122+2*F231*a^2*F113*F222+a^2*F131*F212*F123-2*a^2*F131*F212*F223-a^2*F131*F113*F222-a^2*F131*F213*F122+2*a^2*F131*F213*F222+a^2*F121*F113*F232+a^2*F121*F213*F132-2*a^2*F121*F213*F232-F221*a^2*F112*F133+2*F221*a^2*F112*F233+2*F221*F212*a^2*F133-3*F221*F212*a^2*F233+F221*a^2*F113*F132-2*F221*a^2*F113*F232-2*F221*F213*a^2*F132+3*F221*F213*a^2*F232+a^2*F131*F112*F223-2*F211*a^2*F122*F233-2*F211*F222*a^2*F133+3*F211*F222*a^2*F233-F211*a^2*F123*F132+2*F211*a^2*F123*F232+2*F211*F223*a^2*F132-3*F211*F223*a^2*F232-a^2*F121*F112*F233-a^2*F121*F212*F133+2*a^2*F121*F212*F233-2*a^2*F111*F222*F233-a^2*F111*F123*F232-a^2*F111*F223*F132+2*a^2*F111*F223*F232+a^2*F111*F122*F233+a^2*F111*F222*F133+F211*a^2*F122*F133;
    a=+F231*a^3*F112*F223+F231*a^3*F212*F123-F231*a^3*F212*F223+F231*a^3*F113*F122-F231*a^3*F113*F222-F231*a^3*F213*F122+F231*a^3*F213*F222-a^3*F131*F112*F223-a^3*F131*F212*F123+a^3*F131*F212*F223-a^3*F131*F113*F122+a^3*F131*F113*F222+a^3*F131*F213*F122-a^3*F131*F213*F222+a^3*F121*F212*F133-a^3*F121*F212*F233+a^3*F121*F113*F132-a^3*F121*F113*F232-a^3*F121*F213*F132+a^3*F121*F213*F232+F221*a^3*F112*F133-F221*a^3*F112*F233-F221*a^3*F212*F133+F221*a^3*F212*F233-F221*a^3*F113*F132+F221*a^3*F113*F232+F221*a^3*F213*F132-F221*a^3*F213*F232+a^3*F121*F112*F233-a^3*F111*F222*F133+a^3*F111*F222*F233-a^3*F111*F123*F132+a^3*F111*F123*F232+a^3*F111*F223*F132-a^3*F111*F223*F232-F211*a^3*F122*F133+F211*a^3*F122*F233+F211*a^3*F222*F133-F211*a^3*F222*F233+F211*a^3*F123*F132-F211*a^3*F123*F232-F211*a^3*F223*F132+F211*a^3*F223*F232+a^3*F111*F122*F133-a^3*F111*F122*F233-a^3*F121*F112*F133+a^3*F131*F112*F123-F231*a^3*F112*F123;
    
    sol=roots([a b c d]);
    
    FFF=[];
    j=1;
    for i=1:3,
        if imag(sol(i))==0,
            FFF(:,:,j)=sol(i)*F1+(1-sol(i))*F2;
            j=j+1;
        end
    end
    
    U=[];
    for i=1:size(M,2),
        U=[U; M(1,i)*M(3,i) M(1,i)*M(4,i) M(1,i) M(2,i)*M(3,i) M(2,i)*M(4,i) M(2,i) M(3,i) M(4,i) 1];
    end
    c8=U(:,8);
    c9=U(:,9);
    B=U(:,1:7);
    
    FF=[];
    nor=[];
    dete=[];
    comb=1;
    for n=1:size(FFF,3),
        Fini=FFF(:,:,n);
        for i=1:9,
            for j=i+1:9,
                c1=mod(i,3);
                if c1==0,
                    c1=3;
                end;
                f1=fix((i-1)/3)+1;
                c2=mod(j,3); 
                if c2==0, 
                    c2=3;
                end;
                f2=fix((j-1)/3)+1;
                f8=Fini(f1,c1)/Fini(f2,c2);
                f9=1;
                
                g=-f8*inv(B'*B)*B'*c8-f9*inv(B'*B)*B'*c9;
                nor=[nor norm([g;f8;f9])];
                FF(:,:,comb)=[g(1) g(2) g(3); g(4) g(5) g(6); g(7) f8 f9]./nor(comb);
                dete=[dete det(FF(:,:,comb))];
                
                comb=comb+1;
            end
        end
    end
    [vm,im]=min(abs(dete));
    F=FF(:,:,im);
end

