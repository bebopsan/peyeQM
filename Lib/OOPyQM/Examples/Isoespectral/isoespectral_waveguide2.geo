// Gmsh project created on Tue Apr  9 11:45:19 2013
Point(1) = {2, 0, 0, 1.0};
Point(2) = {2, 1, 0, 1.0};
Point(3) = {3, 1, 0, 1.0};
Point(4) = {2, 2, 0, 1.0};
Point(5) = {1, 2, 0, 1.0};
Point(6) = {1, 3, 0, 1.0};
Point(7) = {0, 3, 0, 1.0};
Point(8) = {0, 2, 0, 1.0};
Point(9) = {1, 1, 0, 1.0};
Point(10) = {1.5, 1.5, 0, 1.0};
Point(11) = {2.5, 1.5, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 11};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 1};
Line(10) = {2, 10};
Line(12) = {10, 5};
Line(11) = {5, 8};
Line(13) = {10, 9};
Line(14) = {11, 4};
Line(15) = {11, 10};
Physical Line(16) = {1, 2, 3, 14, 4, 5, 6, 7, 8, 9};
Physical Line(17) = {10, 15, 12, 11, 13};
Line Loop(18) = {2, 3, 15, -10};
Plane Surface(19) = {18};
Line Loop(20) = {-15, 14, 4, -12};
Plane Surface(21) = {20};
Line Loop(22) = {5, 6, 7, -11};
Plane Surface(23) = {22};
Line Loop(24) = {11, 8, -13, 12};
Plane Surface(25) = {24};
Line Loop(26) = {13, 9, 1, 10};
Plane Surface(27) = {26};
Physical Surface(28) = {19, 21, 23, 25, 27};
Transfinite Line {1, 2, 3, 14, 4, 5, 6, 7, 8, 9} = 4 Using Progression 1;
Transfinite Line {10, 15, 12, 11, 13} = 4 Using Progression 1;
Transfinite Surface {19};
Transfinite Surface {21};
Transfinite Surface {23};
Transfinite Surface {25};
Transfinite Surface {27};
