// Gmsh project created on Mon Apr  8 10:05:31 2013
Point(1) = {1, -1, 0, 1.0};
Point(2) = {1, 1, 0, 1.0};
Point(3) = {-1, 1, 0, 1.0};
Point(4) = {-1, -1, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 4};
Line(5) = {4, 1};
Line Loop(6) = {3, 5, 1, 2};
Plane Surface(7) = {6};
Physical Line(8) = {1, 2, 3, 5};
Physical Surface(9) = {7};
Transfinite Line {5, 1, 2, 3} = 10 Using Progression 1;
Transfinite Surface {7};
Recombine Surface {7};
