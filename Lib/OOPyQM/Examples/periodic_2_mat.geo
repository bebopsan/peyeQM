// Gmsh project created on Sun Apr 14 10:49:30 2013
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 0.5, 0, 1.0};
Point(4) = {1, 1, 0, 1.0};
Point(5) = {0, 1, 0, 1.0};
Point(6) = {0, 0.5, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line(7) = {3, 6};
Line Loop(8) = {1, 2, 7, 6};
Plane Surface(9) = {8};
Line Loop(10) = {7, -5, -4, -3};
Plane Surface(11) = {10};
Physical Line(12) = {1};
Physical Line(13) = {2, 3};
Physical Line(14) = {4};
Physical Line(15) = {5, 6};
Physical Surface(16) = {9};
Physical Surface(17) = {11};
Transfinite Line {1} = 2 Using Progression 1;
Transfinite Line {4} = 2 Using Progression 1;
Transfinite Line {7} = 2 Using Progression 1;
Transfinite Line {2} = 1 Using Progression 1;
Transfinite Line {3} = 1 Using Progression 1;
Transfinite Line {5} = 1 Using Progression 1;
Transfinite Line {6} = 1 Using Progression 1;
Transfinite Surface {9};
Transfinite Surface {11};
Recombine Surface {9, 11};