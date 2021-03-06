// Gmsh project created on Mon Apr  8 14:27:31 2013
Point(1) = {0, 0, 0, 1.0};
Point(2) = {0.4, 0, 0, 1.0};
Point(3) = {0, 0.4, 0, 1.0};
Point(4) = {-0.4, 0, 0, 1.0};
Point(5) = {0, -0.4, 0, 1.0};
Point(6) = {0, -1.5, 0, 1.0};
Point(7) = {1, 0, 0, 1.0};
Point(8) = {0, 1.5, 0, 1.0};
Point(9) = {-1, 0, 0, 1.0};
Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 5};
Line(4) = {5, 2};
Ellipse(5) = {6, 1, 8, 7};
Ellipse(6) = {7, 1, 9, 8};
Ellipse(7) = {8, 1, 6, 9};
Ellipse(8) = {9, 1, 7, 6};
Line(9) = {6, 5};
Line(10) = {7, 2};
Line(11) = {8, 3};
Line(12) = {9, 4};
Physical Line(13) = {5, 6, 7, 8};
Physical Line(14) = {9, 10, 11, 12};
Physical Line(15) = {4, 1, 2, 3};
Line Loop(16) = {5, 10, -4, -9};
Plane Surface(17) = {16};
Line Loop(18) = {6, 11, -1, -10};
Plane Surface(19) = {18};
Line Loop(20) = {7, 12, -2, -11};
Plane Surface(21) = {20};
Line Loop(22) = {8, 9, -3, -12};
Plane Surface(23) = {22};
Line Loop(24) = {4, 1, 2, 3};
Plane Surface(25) = {24};
Physical Surface(26) = {17, 19, 21, 23, 25};
Transfinite Line {5, 6, 7, 8} = 4 Using Progression 1;
Transfinite Line {9, 10, 11, 12} = 4 Using Progression 1;
Transfinite Line {4, 1, 2, 3} = 4 Using Progression 1;
Transfinite Surface {17};
Transfinite Surface {19};
Transfinite Surface {21};
Transfinite Surface {23};
Transfinite Surface {25};
Recombine Surface {17, 19, 21, 23, 25};
