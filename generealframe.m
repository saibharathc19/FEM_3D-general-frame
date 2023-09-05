%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             MATLAB PROGRAM FOR THE ANALYSIS OF 3D GENERAL FRAME         %
%                ===========================================              %
%                   WRITTEN BY:   CHAPPARAM SAI BHARATH                   %
%                    DEPT. OF CIVIL ENGINEERING                           %
%              INDIAN INSTITUTE OF TECHNOLOGY, KHARAGPUR                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
format long
DEG=6;
%DEG= Number of DOF at every node
%All dimensions are in m or N/m2 unless otherwise specified in input file
%% INPUT DATA REQUIRED FOR THE PROGRAM %%
[NN NE NR Node Element RestrainDoF FreeDoF Load]=infun(DEG);
% NN = Number of Nodes
% NE = Number of Elements
% NR = Number of Restrained Degree of Freedoms
% Node= node matrix
    % Node(i,1) = Node Number
    % Node(i,2) = X-Coordinate of ith Node
    % Node(i,3) = Y-Coordinate of ith node
    % Node(i,4) = Z-Coordinate of ith node
% Element = Element connectivity matrix
    % Element(i,1) = Element Number
    % Element(i,2) = Start Node of ith Element
    % Element(i,3) = End Node of ith Elementnt
    % Element(i,4) = Cross Sectional Area of ith Element
    % Element(i,5) = Stiffness of ith Element
    % Element(i,6) = Moment of Inertia of ith Element
    % Element(i,7) = Weight/per unit length of ith Element
    % Element(i,8) = depth of(y) ith Element
    % Element(i,9) = Width of  (Z) of ith Element  
    % Element(i,10) = piossons ratio
    
%% ELEMENT STIFFNESS MATRIX %%
[Kel,Keg,Tel] =FrameElementStiffness(Element,Node,DEG,NE);
EStiffness= Keg;
% EStiffness =  frame Element Stiffness(2*DEG,2*DEG)
% Kel        =  frame Element Stiffness(2*DEG,2*DEG,l)-local
% Keg        =  frame Element Stiffness(2*DEG,2*DEG,l)-Global
% Tel=          frame Element tranformation matrix(2*DEG,2*DEG,element Num)

%% ASSEMBLY STIFFNESS MATRIX %%
Kadd =FrameAssemble(Keg,NE,DEG,NN,Element);
GStiffness=Kadd;
isequal(Kadd,Kadd.');
% GStiffness =  frame Assembeled Stiffness(NE*DEG,NE*DEG)
% Kadd=         frame assembeled Stiffness(NE*DEG,NE*DEG)

%% ASSIGNINMENT OF BOUNDARY CONDITIONS %%
[GStiffnessBC,Kff,Kfc,Kcf,Kcc,X] =FrameStiffnessAssign(Kadd,DEG,NN,RestrainDoF,FreeDoF);
tf = isequal(X,GStiffnessBC);
%tf            =verifying whether bothGstiffness and X are equal
% GStiffnessBC = Global Stiffness Matrix Applying Boundary Condition;
% X            = Global Stiffness Matrix after Applying Kff,Kcc,Kfc,Kcf;

 %% ASSIGNMENT OF LOADS ON THE RESPECTIVE DOF %%  
LoadBC=zeros(DEG*NN,2);

for j=1:1:DEG*NN
LoadBC(j,1)=j;
end
r=1;
for j=1:1:NN
for i=1:1:DEG
    LoadBC(r,2)=Load(i,j);
    r=r+1;
end
end
LoadBC;
LoadBC(RestrainDoF(1,:),:)=[];
% LoadBC = Loadings After Applying Boundary Conditions

%% ANALYSIS OF THE STRUCTURE %%
u=(Kfc*(RestrainDoF(2,:).'));
DispUnConstrained= Kff\(LoadBC(:,2)-u);
Uf(:,1)=(FreeDoF(1,:).');
Uf(:,2)=(DispUnConstrained);
Uc(:,1)=(RestrainDoF(1,:).');
Uc(:,2)=(RestrainDoF(2,:).');

%Uf= Matrix containg free displacement to respective DOF 
%Uf(:,1)=Free DOF
%Uf(:,2)=Free DOF Displacement
%Uc(:,1)=Restrained DOF
%Uc(:,2)=Restrained DOF Displacement

U=[Uf;Uc];        
Displacement=zeros(DEG*NN,2);
m=1;
for i=1:1:NN*DEG
    for j=1:1:NN*DEG
        if(U(j,1)==i)
            Displacement(m,2)=U(j,2);
            m=m+1;
        end
    end
end
% Displacement= GLobal displacement at each DOF

Rf(:,2)=X*U(:,2);
Rf(:,1)=U(:,1);
ReactionForce=zeros(DEG*NN,2);
m=1;
for i=1:1:NN*DEG
    for j=1:1:NN*DEG
        if(Rf(j,1)==i)
            ReactionForce(m,2)=Rf(j,2);
            m=m+1;
        end
    end
end


[F,Disl] = ElementalForces(Element,DEG,NE,U,Tel,Kel);
% F size is NE X DEG*2
% Disl size is NE X DEG*2
% F containts the forces of elements at their nodes
% Disl contains the local displacements of elements at their nodes


        
        
        
        
        
        
        
        
        