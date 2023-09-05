function Kadd=FrameAssemble(Keg,NE,DEG,NN,Element)

% Assemble the Element Stiffness/Mass Matrix 
% member to Global Stiffness/Mass Matrix
% Kadd    Global Stiffness additive/Mass Matrix
% Keg    GlobalElement Stiffness/Mass Matrix
% X    Row Matrix for DOF of each element
                 

X=zeros(1,2*DEG);
Kadd=zeros(DEG*NN,DEG*NN);

for l=1:1:NE
    
 for j= DEG-1:-1:0
    X(DEG-j)=DEG*Element(l,2)-j;
    X(2*DEG-j)=DEG*Element(l,3)-j;
 end
 
 for m=1:1:2*DEG
    for n=1:1:2*DEG
        p=X(1,m);
        q=X(1,n);
        Kadd(p,q)=Kadd(p,q)+Keg(m,n,l);
    end
 end
 
end
    
    
