Point(1) = {0.000000, 0.000000, 0, 1.0};
Point(2) = {1.000000, 0.000000, 0, 1.0};
Point(3) = {1.000000, 1.000000, 0, 1.0};
Point(4) = {0.000000, 1.000000, 0, 1.0};
Point(5) = {0.500000, 0.500000, 0, 1.0};
Point(6) = {0.420000, 0.420000, 0, 1.0};
Point(7) = {0.580000, 0.420000, 0, 1.0};
Point(8) = {0.580000, 0.580000, 0, 1.0};
Point(9) = {0.420000, 0.580000, 0, 1.0};
Point(10) = {0.358579, 0.358579, 0, 1.0};
Point(11) = {0.641421, 0.358579, 0, 1.0};
Point(12) = {0.641421, 0.641421, 0, 1.0};
Point(13) = {0.358579, 0.641421, 0, 1.0};
Line(2) = {1, 2};
Transfinite Line {2} = 2 Using Progression 1;
Line(3) = {2, 3};
Transfinite Line {3} = 2 Using Progression 1;
Line(4) = {3, 4};
Transfinite Line {4} = 2 Using Progression 1;
Line(5) = {4, 1};
Transfinite Line {5} = 2 Using Progression 1;
Line(9) = {6, 7};
Transfinite Line {9} = 2 Using Progression 1;
Line(10) = {7, 8};
Transfinite Line {10} = 2 Using Progression 1;
Line(11) = {8, 9};
Transfinite Line {11} = 2 Using Progression 1;
Line(12) = {9, 6};
Transfinite Line {12} = 2 Using Progression 1;
Circle(13) = {10, 5, 11};
Transfinite Line {13} = 2 Using Progression 1;
Circle(14) = {11, 5, 12};
Transfinite Line {14} = 2 Using Progression 1;
Circle(15) = {12, 5, 13};
Transfinite Line {15} = 2 Using Progression 1;
Circle(16) = {13, 5, 10};
Transfinite Line {16} = 2 Using Progression 1;
Line(17) = {6, 10};
Transfinite Line {17} = 2 Using Progression 1;
Line(18) = {7, 11};
Transfinite Line {18} = 2 Using Progression 1;
Line(19) = {8, 12};
Transfinite Line {19} = 2 Using Progression 1;
Line(20) = {9, 13};
Transfinite Line {20} = 2 Using Progression 1;
Line(21) = {10, 1};
Transfinite Line {21} = 2 Using Progression 1;
Line(22) = {11, 2};
Transfinite Line {22} = 2 Using Progression 1;
Line(23) = {12, 3};
Transfinite Line {23} = 2 Using Progression 1;
Line(24) = {13, 4};
Transfinite Line {24} = 2 Using Progression 1;
Line Loop(13) = {9, 10, 11, 12};
Line Loop(14) = {-9, 17, 13, -18};
Line Loop(15) = {-10, 18, 14, -19};
Line Loop(16) = {-11, 19, 15, -20};
Line Loop(17) = {-12, 20, 16, -17};
Line Loop(18) = {-13, 21, 2, -22};
Line Loop(19) = {-14, 22, 3, -23};
Line Loop(20) = {-15, 23, 4, -24};
Line Loop(21) = {-16, 24, 5, -21};
Plane Surface(14) = {13};
Transfinite Surface {14};
Plane Surface(15) = {14};
Transfinite Surface {15};
Plane Surface(16) = {15};
Transfinite Surface {16};
Plane Surface(17) = {16};
Transfinite Surface {17};
Plane Surface(18) = {17};
Transfinite Surface {18};
Plane Surface(19) = {18};
Transfinite Surface {19};
Plane Surface(20) = {19};
Transfinite Surface {20};
Plane Surface(21) = {20};
Transfinite Surface {21};
Plane Surface(22) = {21};
Transfinite Surface {22};
Physical Surface(31) = {14, 15, 16, 17, 18 };
Physical Surface(30) = {19, 20, 21, 22 };
Recombine Surface {14};
Recombine Surface {15};
Recombine Surface {16};
Recombine Surface {17};
Recombine Surface {18};
Recombine Surface {19};
Recombine Surface {20};
Recombine Surface {21};
Recombine Surface {22};
Physical Line(32) = {2};
Physical Line(33) = {3};
Physical Line(34) = {4};
Physical Line(35) = {5};
