function length = PlaneElementLength(x1,y1,z1,x2,y2,z2)

%PlaneElementLength         Returns the Length of a 2D Frame Element
%           (x1,y1)         Coordinate of the First Node
%           (x2,y2)         Coordinate of the Second Node   
length=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);