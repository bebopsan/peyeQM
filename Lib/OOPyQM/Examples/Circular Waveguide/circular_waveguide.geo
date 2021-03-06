// Gmsh project created on Mon Apr  8 14:27:31 2013
Point(1) = {0, 0, 0, 1.0};
Point(2) = {0.8, -0.8, 0, 1.0};
Point(3) = {0.8, 0.8, 0, 1.0};
Point(4) = {-0.8, 0.8, 0, 1.0};
Point(5) = {-0.8, -0.8, 0, 1.0};
Point(6) = {1.41421356237, -1.41421356237, 0, 1.0};
Point(7) = {1.41421356237, 1.41421356237, 0, 1.0};
Point(8) = {-1.41421356237, 1.41421356237, 0, 1.0};
Point(9) = {-1.41421356237, -1.41421356237, 0, 1.0};
Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 5};
Line(4) = {5, 2};
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};
Line(9) = {2, 6};
Line(10) = {3, 7};
Line(11) = {4, 8};
Line(12) = {5, 9};
Physical Line(13) = {4, 1, 2, 3};
Physical Line(14) = {5, 6, 7, 8};
Physical Line(15) = {9, 10, 11, 12};
Line Loop(16) = {9, 5, -10, -1};
Plane Surface(17) = {16};
Line Loop(18) = {10, 6, -11, -2};
Plane Surface(19) = {18};
Line Loop(20) = {11, 7, -12, -3};
Plane Surface(21) = {20};
Line Loop(22) = {12, 8, -9, -4};
Plane Surface(23) = {22};
Line Loop(24) = {4, 1, 2, 3};
Plane Surface(25) = {24};
Physical Surface(26) = {17, 19, 21, 23, 25};
Transfinite Line {4, 1, 2, 3} = 10 Using Progression 1;
Transfinite Line {9, 10, 11, 12} = 5 Using Progression 1;
Transfinite Line {5, 6, 7, 8} = 10 Using Progression 1;
Transfinite Surface {17};
Transfinite Surface {19};
Transfinite Surface {21};
Transfinite Surface {25};
Transfinite Surface {23};
Recombine Surface {17, 19, 21, 23, 25};
