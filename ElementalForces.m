function [F,Disl] = ElementalForces(Element,DEG,NE,U,Tel,Kel)
Dis=zeros(2*DEG,1);        
X=zeros(1,2*DEG);
F=zeros(2*DEG,NE);
Disl=zeros(2*DEG,1);

[t,k]=size(U);
i=1;
for l=1:1:NE
    
 for j= DEG-1:-1:0
    X(DEG-j)=DEG*Element(l,2)-j;
    X(2*DEG-j)=DEG*Element(l,3)-j;
 end
 
 for m=1:1:2*DEG
    for n=1:1:t
        if(U(n,1)==X(1,m))
            Dis(i)=U(n,2);
            i=i+1;
        end
    end
 end
 i=1;
 
  Disl(:,l)=Tel(:,:,l)\Dis;
  F(:,l)= Kel(:,:,l)*Disl(:,l);
  
end

        