// Gmsh project created on Fri Mar 29 13:26:57 2024
SetFactory("OpenCASCADE");
// Fuselage assembly 
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {2.8, 0, 0, 1.0};
//+
Point(3) = {-2.8, 0, 0, 1.0};
//+
Point(4) = {0, -2.8, 0, 1.0};
//+
Point(5) = {0, 2.8, 0, 1.0};
//+
Point(6) = {0, 2.98, 0, 1.0};
//+
Point(7) = {0, -2.98, 0, 1.0};
//+
Point(8) = {-2.98, 0, 0, 1.0};
//+
Point(9) = {2.98, 0, 0, 1.0};
//+
Point(10) = {2.98, -2, 0, 1.0};
//+
Point(11) = {2.98, -2.98, 0, 1.0};
//+
Point(12) = {2, -2.98, 0, 1.0};
//+
Point(13) = {-2, -2.98, 0, 1.0};
//+
Point(14) = {-2.98, -2.98, 0, 1.0};
//+
Point(15) = {-2.98, -2, 0, 1.0};
//+
Circle(1) = {8, 1, 6};
//+
Circle(2) = {3, 1, 5};
//+
Circle(3) = {5, 1, 2};
//+
Circle(4) = {2, 1, 4};
//+
Circle(5) = {4, 1, 3};
//+
Circle(6) = {6, 1, 9};
//+
BSpline(7) = {8, 15, 14, 13, 7};
//+
BSpline(8) = {7, 12, 11, 10, 9};
// Right wing assembly
// Top part
//+ 
Point(16) = {35.35, 2.9, 0, 1.0};
//+
Point(17) = {2.9, -0.7, 0, 1.0};
//+
Point(18) = {2.8, -2.4, 0, 1.0};
//+
Point(19) = {5, -0.6, 0, 1.0};
//+
Point(20) = {9, -0.5, 0, 1.0};
//+
Point(21) = {16, -0.1, 0, 1.0};
//+
Point(22) = {25, 0.6, 0, 1.0};
//+
BSpline(9) = {17, 19, 20, 21, 22, 16};
// Down part
//+
Point(23) = {25, -0.3, 0, 1.0};
//+
Point(24) = {35.35, 2.7, 0, 1.0};
//+
Point(25) = {5, -1.5, 0, 1.0};
//+
Point(26) = {16, -0.8, 0, 1.0};
//+
Point(27) = {9, -1, 0, 1.0};
//+
Bezier(10) = {18, 25, 27, 26, 23, 24};
//+
Spline(11) = {17, 18};
// Right winglet assembly
//+
Point(28) = {35.6, 3, 0, 1.0};
//+
Point(29) = {35.75, 3.55, 0, 1.0};
//+
Point(30) = {35.75, 4.1, 0, 1.0};
//+
Point(31) = {35.8, 4.6, 0, 1.0};
//+
Point(32) = {35.9, 4.8, 0, 1.0};
//+
Point(33) = {35.6, 2.8, 0, 1.0};
//+
Point(34) = {35.8, 3, 0, 1.0};
//+
Point(35) = {35.9, 3.4, 0, 1.0};
//+
Point(36) = {35.9, 3.8, 0, 1.0};
//+
Point(37) = {35.9, 4.2, 0, 1.0};
//+
Point(38) = {35.9, 4.6, 0, 1.0};
//+
BSpline(12) = {24, 33, 34, 35, 36, 37, 38, 32, 31, 30, 29, 28, 16};
// Left wing assembly
// Top part
//+ 
Point(39) = {-35.35, 2.9, 0, 1.0};
//+
Point(40) = {-2.9, -0.7, 0, 1.0};
//+
Point(41) = {-2.8, -2.4, 0, 1.0};
//+
Point(42) = {-5, -0.6, 0, 1.0};
//+
Point(43) = {-9, -0.5, 0, 1.0};
//+
Point(44) = {-16, -0.1, 0, 1.0};
//+
Point(45) = {-25, 0.6, 0, 1.0};
//+
BSpline(13) = {40, 42, 43, 44, 45, 39};
// Down part
//+
Point(46) = {-25, -0.3, 0, 1.0};
//+
Point(47) = {-35.35, 2.7, 0, 1.0};
//+
Point(48) = {-5, -1.5, 0, 1.0};
//+
Point(49) = {-16, -0.8, 0, 1.0};
//+
Point(50) = {-9, -1, 0, 1.0};
//+
Bezier(14) = {41, 48, 50, 49, 46, 47};
//+
Spline(15) = {40, 41};
// Left winglet assembly
//+
Point(51) = {-35.6, 3, 0, 1.0};
//+
Point(52) = {-35.75, 3.55, 0, 1.0};
//+
Point(53) = {-35.75, 4.1, 0, 1.0};
//+
Point(54) = {-35.8, 4.6, 0, 1.0};
//+
Point(55) = {-35.9, 4.8, 0, 1.0};
//+
Point(56) = {-35.6, 2.8, 0, 1.0};
//+
Point(57) = {-35.8, 3, 0, 1.0};
//+
Point(58) = {-35.9, 3.4, 0, 1.0};
//+
Point(59) = {-35.9, 3.8, 0, 1.0};
//+
Point(60) = {-35.9, 4.2, 0, 1.0};
//+
Point(61) = {-35.9, 4.6, 0, 1.0};
//+
BSpline(16) = {47, 56, 57, 58, 59, 60, 61, 55, 54, 53, 52, 51, 39};
// Vertical stab
//+
Point(62) = {0, 11.1, 0, 1.0};
//+
Point(63) = {0.2, 10.5, 0, 1.0};
//+
Point(64) = {-0.2, 10.5, 0, 1.0};
//+
Point(65) = {-0.3, 8, 0, 1.0};
//+
Point(66) = {0.3, 8, 0, 1.0};
//+
Point(67) = {-0.4, 4, 0, 1.0};
//+
Point(68) = {0.4, 4, 0, 1.0};
//+
Point(69) = {-0.5, 2.9, 0, 1.0};
//+
Point(70) = {0.5, 2.9, 0, 1.0};
//+
Spline(17) = {69, 67, 65, 64, 62};
//+
Spline(18) = {62, 63, 66, 68, 70};
//+
Spline(19) = {69, 70};
// Horizontal Stab Right
//+
Point(71) = {2.4, 1.6, 0, 1.0};
//+
Point(72) = {2.1, 2, 0, 1.0};
//+
Point(73) = {4, 2.1, 0, 1.0};
//+
Point(74) = {6.8, 2.3, 0, 1.0};
//+
Point(75) = {4, 1.7, 0, 1.0};
//+
Point(76) = {6.8, 2.1, 0, 1.0};
//+
Point(77) = {8.6, 2.5, 0, 1.0};
//+
Spline(20) = {71, 72};
//+
Spline(21) = {71, 75, 76, 77, 74, 73, 72};
// Horizontal Stab Left
//+
Point(78) = {-2.4, 1.6, 0, 1.0};
//+
Point(79) = {-2.1, 2, 0, 1.0};
//+
Point(80) = {-4, 2.1, 0, 1.0};
//+
Point(81) = {-6.8, 2.3, 0, 1.0};
//+
Point(82) = {-4, 1.7, 0, 1.0};
//+
Point(83) = {-6.8, 2.1, 0, 1.0};
//+
Point(84) = {-8.6, 2.5, 0, 1.0};
//+
Spline(22) = {78, 79};
//+
Spline(23) = {78, 82, 83, 84, 81, 80, 79};
//Creating the surfaces
// Right wing
//+ 
Curve Loop(1) = {9, -12, -10, -11};
//+
Plane Surface(1) = {1};
// Left wing
//+
Curve Loop(2) = {13, -16, -14, -15};
//+
Plane Surface(2) = {2};
// Vertical stab
//+
Curve Loop(3) = {19, -18, -17};
//+
Plane Surface(3) = {3};
// Horizontal stab right
//+
Curve Loop(4) = {20, -21};
//+
Plane Surface(4) = {4};
// Horizontal stab left
//+
Curve Loop(5) = {22, -23};
//+
Plane Surface(5) = {5};
// Fuselage
//+
Curve Loop(6) = {7, 8, -6, -1};
//+
Curve Loop(7) = {5, 2, 3, 4};
//+
Plane Surface(6) = {6, 7};
//+
BooleanUnion{ Surface{1}; Delete; }{ Surface{6}; Delete; }
//+
BooleanUnion{ Surface{2}; Delete; }{ Surface{6}; Delete; }
//+
BooleanUnion{ Surface{3}; Delete; }{ Surface{6}; Delete; }
//+
BooleanUnion{ Surface{4}; Delete; }{ Surface{6}; Delete; }
//+
BooleanUnion{ Surface{5}; Delete; }{ Surface{6}; Delete; }
//+
Coherence;
//+
Curve Loop(16) = {37, -36, -35, 55, 28, 29, -30, 56, 57, 48, 58, -43, -42, 59, -53, 60, 54};
//+
Curve Loop(17) = {61};
//+
Plane Surface(17) = {16, 17};
//+
Recursive Delete {
  Surface{6}; Surface{12}; Surface{16}; Surface{14}; Surface{10}; Surface{8}; Surface{9};   Curve{52}; Curve{50}; Surface{7}; Surface{15}; Surface{11}; Surface{13}; 
}


