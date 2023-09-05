
function [NN, NE NR Node Element RestrainDoF FreeDoF Load]=infun(DEG)

finput = fopen('PlaneFrameInp.in','r');
Dummy = fgets(finput); Dummy = fgets(finput); Dummy = fgets(finput); 
Dummy = fgets(finput); Dummy = fgets(finput); Dummy = fgets(finput); 
Dummy = fgets(finput); Dummy = fgets(finput); Dummy = fgets(finput); 
% Preallocate all the matrixes

Temp = str2num(fgets(finput));

% NN = Number of Nodes
% NE = Number of Elements
% NR = Number of Restrained Degree of Freedoms
%Element= 

[NN, NE, NR] = deal(Temp(1), Temp(2), Temp(3));

Element=zeros(NE,7);
Node=zeros(NN,4);
RestrainDoF =zeros(2,NR);
FreeDoF=zeros(1,(DEG*NN)-NR);

% INPUT COORDINATES OF EACH NODE
Dummy = fgets(finput); Dummy = fgets(finput); Dummy = fgets(finput);
Dummy = fgets(finput); Dummy = fgets(finput);
for i = 1:NN
    Temp = str2num(fgets(finput));
    % Node= node matrix
    % Node(i,1) = Node Number
    % Node(i,2) = X-Coordinate of ith Node
    % Node(i,3) = Y-Coordinate of ith node
    Node(i,:)=deal(Temp);
    
end

% INPUT CONNECTIVITY, GEOMETRY DATA AND MATERIAL PROPERTIES
Dummy = fgets(finput); Dummy = fgets(finput); Dummy = fgets(finput);
Dummy = fgets(finput); Dummy = fgets(finput); 

for i=1:NE
    Temp = str2num(fgets(finput));
    % Element(i,1) = Element Number
    % Element(i,2) = Start Node of ith Element
    % Element(i,3) = End Node of ith Element
    % Element(i,4) = Length of ith Element
    % Element(i,5) = Cross Sectional Area of ith Element
    % Element(i,6) = Stiffness of ith Element
    % Element(i,7) = Moment of Inertia of ith Element
    % Element(i,8) = Weight/per unit length of ith Element
    % Element(i,9) = piossons ratio
    % Element(i,10) = depth of(y) ith Element
    % Element(i,11) = Width of  (Z) of ith Element     
    % Element(i,12) = Moment of Inertia(Iz) of ith Element
    % Element(i,13) = Moment of Inertia(Iy) of ith Element
    % Element(i,14) = Moment of Inertia(Ix) of ith Element
    % Element(i,15) = Shear Modulus of ith Element
    [Element(i,1), Element(i,2), Element(i,3), r,Element(i,6), ...
        Element(i,7),Element(i,8),Element(i,9),Element(i,10),Element(i,11)] = deal(Temp(1), Temp(2), Temp(3), Temp(4), ...
        Temp(5), Temp(6), Temp(7),Temp(8),Temp(9),Temp(10));
    Element(i,6)=(Element(i,6)*25000*1000)/(28500000000);
    Element(i,8)=(25/9.8);
    
    Element(i,5)=Element(i,11)*Element(i,10);
    
    Element(i,4) = PlaneElementLength(Node(Element(i,2),2), Node(Element(i,2),3),Node(Element(i,2),4),...
    Node(Element(i,3),2),Node(Element(i,3),3),Node(Element(i,3),4));
    Element(i,12)= (Element(i,10)*Element(i,9).^3)/12;
    Element(i,13)= (Element(i,9)*Element(i,10).^3)/12;
    %Element(i,12)=0.00006;
    %Element(i,13)=0.00006;
    Element(i,14)= Element(i,12)+Element(i,13);
    Element(i,15)= Element(i,6)/(2*(1+Element(i,9)));
    %THe following are used for specified problem
    
   %Element(i,14)=1000;
   % Element(i,15)=84000000000;
end
% INPUT RESTRAINED DEGREE OF FREEDOM
Dummy = fgets(finput); Dummy = fgets(finput); Dummy = fgets(finput);
%RestrainDoF(1,:) is for free DOF Numbers
%RestrainDoF(2,:) is for free DOF Displacments
RestrainDoF(1,:) = str2num(fgets(finput));
RestrainDoF(2,:) = str2num(fgets(finput));
% EVALUATIONG FREE DEGREE OF FREEDOM
%they are used to get the free DOF for evaluationg the final K matrix
[m,n]=size(RestrainDoF);
loc=1;
k=1;
for i=1:1:DEG*NN
    for j=1:1:n
        if(i==RestrainDoF(1,j))
            loc=0;
        end
    end
    if(loc~=0)
        FreeDoF(1,k)=i;
            k=k+1;
    end
    loc=1;
end

% INPUT EXTERNAL LOADS ACTING ON PLANE TRUSS
        % Load = Loading on all Degree of Freedom
        % Load(rows) = The nmber of rows will be num ber of DOF
        %lOAD(Column)=The number of columns will be number of nodes

        Dummy = fgets(finput); Dummy = fgets(finput); Dummy = fgets(finput);
        Dummy = fgets(finput); Dummy = fgets(finput); Dummy = fgets(finput); 
        Dummy = fgets(finput);
        
        Load = zeros(DEG,NN);
        LoadEl = zeros(2*DEG,NE);

        for i = 1:NN
            Temp = deal(str2num(fgets(finput)));
            Ele = Temp(1);
            Load(:,Ele) = Temp(2:(DEG)+1);
            %for j=DEG-1:-1:0
            %Load(DEG*Element(Ele,2)-j,1) = Load(DEG*Element(Ele,2)-j,1)+LoadEl(DEG-j,i);
            %Load(DEG*Element(Ele,3)-j,1) = Load(DEG*Element(Ele,3)-j,1)+LoadEl(2*DEG-j,i);
            %end
        end
        Element(:,5)
        Element(:,6)
end






