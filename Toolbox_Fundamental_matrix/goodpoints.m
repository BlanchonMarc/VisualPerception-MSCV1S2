%goodpoints    
%    [M]=goodpoints(Mini,w,threshold)
%
% by X. Armangue
% (c) Mr3D - University of Girona, September 2002
%
function [M]=goodpoints(Mini,w,threshold)

if ((size(Mini,1)~=4) | (size(w,1)~=1) | (size(Mini,2)~=size(w,2))),
   disp('Error: parametres incorrectes')
else
   M=[];
   for i=1:size(Mini,2),
      if w(i)>=threshold,
         M=[M Mini(:,i)];
      end      
   end
end
