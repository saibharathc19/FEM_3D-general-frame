function [Kel,Keg,Tel] =FrameElementStiffness(Element,Node,DEG,NE)

% FrameElementStiffness    Creats the Element Stiffness Matrix for a
%                                Frame Element
%                          E    Modulous of Elasticity of the member
%                          I    Moment of Inertia of Frame Element
%                          A    Cross Sectional Area of the member
%                          L    Length of the member
%                     (x1,y1)   Coordinate of Node at start end
%                     (x2,y2)   Coordinate of Node at start end



% T,Tr     Transformation Matrix
%Tel    Final element Transformation matrix for each element
%Tfin    Final element Transformation matrix
% Kel   Element Stiffness in Local Coordinate System
% Keg   Element Stiffness in Global Coordinate System
% 

%preallocate matrixes

T=zeros(3,3);
Tel=zeros(2*DEG,2*DEG,NE);
Kel=zeros(2*DEG,2*DEG,NE);
Keg=zeros(2*DEG,2*DEG,NE);

%alloting zeros for transformation matrix
if(DEG==3)
    dim=2;
   
else
    dim=3;
   
end
%Formation of all global axis
Xglo=[1;0;0];
Yglo=[0;1;0];
Zglo=[0;0;1];




for l=1:1:NE
%S is that to store all the values of ae/l....
s(1)=(Element(l,6)*Element(l,5))/Element(l,4);
s(2)=(Element(l,6)*Element(l,12)*12)/(Element(l,4).^3);
s(3)=(Element(l,6)*Element(l,13)*12)/(Element(l,4).^3);
s(4)=(Element(l,15)*Element(l,14))/Element(l,4);
s(5)=s(3)*(Element(l,4).^2)/3;
s(6)=s(2)*(Element(l,4).^2)/3;
s(7)=-1*s(3)*(Element(l,4))/2;
s(8)=s(2)*(Element(l,4))/2;
%dividing K matrix into 4 matrixes Ka,Kb,Kc,Kd

Ka=zeros(6,6);
for i=1:1:6
    for j=1:1:6
        if(i==j)
            Ka(i,j)=s(i);
        end
        if(i==3)
            if(j==5)
            Ka(i,j)=s(7);
            Ka(j,i)=Ka(i,j);
            end
        end
        if(i==2)
            if(j==6)
            Ka(i,j)=s(8);
            Ka(j,i)=Ka(i,j);
            end
        end
    end
end
%as there is a change in sign for values of K a and Kb introducing - signs
s(1)=-1*s(1);
s(2)=-1*s(2);
s(3)=-1*s(3);
s(4)=-1*s(4);
s(5)=-1*s(3)*(Element(l,4).^2)/6;
s(6)=-1*s(2)*(Element(l,4).^2)/6;
Kb=zeros(6,6);
for i=1:1:6
    for j=1:1:6
        if(i==j)
            Kb(i,j)=s(i);
        end
        if(i==3 && j==5)
            Kb(i,j)=s(7);
            Kb(j,i)=-1*Kb(i,j);
        end
        if(i==2 && j==6)
            Kb(i,j)=s(8);
            Kb(j,i)=-1*Kb(i,j);
        end
    end
end
%changing the values accordingly for (3,5) and(2,6)
Kc=Kb;
Kc(3,5)=-1*Kb(3,5);
Kc(2,6)=-1*Kb(2,6);
Kc(5,3)=-1*Kb(5,3);
Kc(6,2)=-1*Kb(6,2);
Kd=Ka;
Kd(3,5)=-1*Ka(3,5);
Kd(2,6)=-1*Ka(2,6);
Kd(5,3)=-1*Ka(5,3);
Kd(6,2)=-1*Ka(6,2);
Ka;
Kb;
Kc;
Kd;
K=[Ka Kb;Kc Kd];
isequal(K,K.');
if(dim==2)
    K(:,3)=[];
    K(:,3)=[];
    K(:,3)=[];
    K(:,6)=[];
    K(:,6)=[];
    K(:,6)=[];
    K(3,:)=[];
    K(3,:)=[];
    K(3,:)=[];
    K(6,:)=[];
    K(6,:)=[];
    K(6,:)=[];
end
Kel(:,:,l)=K;
%Formation of transformation matrix
    for i=1:1:3
       if(i==1)
           
        T(1,i)= (Node(Element(l,3),2)-Node(Element(l,2),2))/(Element(l,4));
        T(2,i)= (Node(Element(l,3),3)-Node(Element(l,2),3))/(Element(l,4));
        T(3,i)= (Node(Element(l,3),4)-Node(Element(l,2),4))/(Element(l,4));
       end
             
        if(i==2)
           Z=cross(T(:,1),Xglo);
           if(Z==0)
              if(Xglo==T(:,1))
                  Z=Zglo;
              else
                  Z=-Zglo;
              end
           end
           
           T(1,3)=(dot(Z,Xglo))/norm(Z);
           T(2,3)=(dot(Z,Yglo))/norm(Z);
           T(3,3)=(dot(Z,Zglo))/norm(Z);
           
        end
           
           if(i==3) 
             Z=cross(T(:,3),T(:,1)); 
            if(Z==0)
              if(Xglo==T(:,1))
                  Z=Yglo;
              else
                  Z=-Yglo;
              end
           end
             T(1,2)=(dot(Z,Xglo))/norm(Z);
           T(2,2)=(dot(Z,Yglo))/norm(Z);
           T(3,2)=(dot(Z,Zglo))/norm(Z);
             
           end
       
    end
    
    if(dim==2)
        o=zeros(3,3);
        T(3,:)=[0 0 1];
        T(:,3)=[0; 0; 1];
        Tfin=[T o;o T];
    else
       o=zeros(3,3);
       m=zeros(6,6);
       Tr=[T o;o T];
       Tfin=[Tr m;m Tr];
       
    end
     
    Tel(:,:,l)=Tfin;
    Keg(:,:,l)=(Tfin)*K*(Tfin.');
end
    
end
