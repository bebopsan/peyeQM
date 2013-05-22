// Gmsh project created on Mon Apr  8 10:05:31 2013
Point(1) = {2, -1, 0, 1.0};
Point(2) = {2, 1, 0, 1.0};
Point(3) = {-2, 1, 0, 1.0};
Point(4) = {-2, -1, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};
Physical Line(7) = {1, 2, 3, 4};
Physical Surface(8) = {6};
Transfinite Line {2, 4} = 20 Using Progression 1;
Transfinite Line {1, 3} = 10 Using Progression 1;
Transfinite Surface {6};
Recombine Surface {6};
