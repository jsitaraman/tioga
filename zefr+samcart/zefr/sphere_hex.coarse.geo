size1 = .3;  // .19
size2 = .7;  // .55
R1 = .0625;  // Radius of inner sphere
R2 = .18;  // Radius of outer sphere
NN = 6;  // Number of points in each direction on each spherical surface
NL = 4;  // Number of layers between surfaces
prog1 = 1.2;  // Geometric progression factor for layer width

/* ---- Inner Sphere Surface ---- */

xc = 0;
yc = 0;
zc = 0;
Point(1) = {xc, yc, zc, size1};

r1 = R1 / Sqrt(3);
Point(2) = {xc-r1, yc-r1, zc-r1, size1};
Point(3) = {xc+r1, yc-r1, zc-r1, size1};
Point(4) = {xc+r1, yc+r1, zc-r1, size1};
Point(5) = {xc-r1, yc+r1, zc-r1, size1};

Point(6) = {xc-r1, yc-r1, zc+r1, size1};
Point(7) = {xc+r1, yc-r1, zc+r1, size1};
Point(8) = {xc+r1, yc+r1, zc+r1, size1};
Point(9) = {xc-r1, yc+r1, zc+r1, size1};

// 'Bottom' circles
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

// 'Top' circles
Circle(5) = {6,1,7};
Circle(6) = {7,1,8};
Circle(7) = {8,1,9};
Circle(8) = {9,1,6};

// 'Vertical' circles
Circle(9) = {2,1,6};
Circle(10) = {3,1,7};
Circle(11) = {4,1,8};
Circle(12) = {5,1,9};

Transfinite Line {1:12} = NN+1 Using Progression 1.0;

Line Loop (1) = {1,2,3,4};
Line Loop (2) = {5,6,7,8};
Line Loop (3) = {4,9,-8,-12};
Line Loop (4) = {2,11,-6,-10};
Line Loop (5) = {1,10,-5,-9};
Line Loop (6) = {3,12,-7,-11};

Ruled Surface (1) = {1};
Ruled Surface (2) = {2};
Ruled Surface (3) = {3};
Ruled Surface (4) = {4};
Ruled Surface (5) = {5};
Ruled Surface (6) = {6};

Transfinite Surface {1:6};
Recombine Surface {1:6};


/* ---- Outer Sphere Surface ---- */

Point(11) = {xc, yc, zc, size1};

r2 = R2 / Sqrt(3);
Point(12) = {xc-r2, yc-r2, zc-r2, size1};
Point(13) = {xc+r2, yc-r2, zc-r2, size1};
Point(14) = {xc+r2, yc+r2, zc-r2, size1};
Point(15) = {xc-r2, yc+r2, zc-r2, size1};

Point(16) = {xc-r2, yc-r2, zc+r2, size1};
Point(17) = {xc+r2, yc-r2, zc+r2, size1};
Point(18) = {xc+r2, yc+r2, zc+r2, size1};
Point(19) = {xc-r2, yc+r2, zc+r2, size1};

// 'Bottom' circles
Circle(21) = {12,11,13};
Circle(22) = {13,11,14};
Circle(23) = {14,11,15};
Circle(24) = {15,11,12};

// 'Top' circles
Circle(25) = {16,11,17};
Circle(26) = {17,11,18};
Circle(27) = {18,11,19};
Circle(28) = {19,11,16};

// 'Vertical' circles
Circle(29) = {12,11,16};
Circle(30) = {13,11,17};
Circle(31) = {14,11,18};
Circle(32) = {15,11,19};

Transfinite Line {21:32} = NN+1 Using Progression 1.0;

Line Loop (11) = {21,22,23,24};
Line Loop (12) = {25,26,27,28};
Line Loop (13) = {24,29,-28,-32};
Line Loop (14) = {22,31,-26,-30};
Line Loop (15) = {21,30,-25,-29};
Line Loop (16) = {23,32,-27,-31};

Ruled Surface (11) = {11};
Ruled Surface (12) = {12};
Ruled Surface (13) = {13};
Ruled Surface (14) = {14};
Ruled Surface (15) = {15};
Ruled Surface (16) = {16};

Transfinite Surface {11:16};
Recombine Surface {11:16};


/* ---- Connecting Lines ---- */

Line(41) = {2,12};
Line(42) = {3,13};
Line(43) = {4,14};
Line(44) = {5,15};
Line(45) = {6,16};
Line(46) = {7,17};
Line(47) = {8,18};
Line(48) = {9,19};

Transfinite Line {41:48} = NL+1 Using Progression prog1;

/* ---- Connecting Planes ---- */

Line Loop (21) = {-1,41,21,-42};
Line Loop (22) = {-2,42,22,-43};
Line Loop (23) = {-3,43,23,-44};
Line Loop (24) = {-4,44,24,-41};

Line Loop (25) = {5,46,-25,-45};
Line Loop (26) = {6,47,-26,-46};
Line Loop (27) = {7,48,-27,-47};
Line Loop (28) = {8,45,-28,-48};

Line Loop (29) = {-9,41,29,-45};
Line Loop (30) = {-10,42,30,-46};
Line Loop (31) = {-11,43,31,-47};
Line Loop (32) = {-12,44,32,-48};

Ruled Surface (31) = {21};
Ruled Surface (32) = {22};
Ruled Surface (33) = {23};
Ruled Surface (34) = {24};
Ruled Surface (35) = {25};
Ruled Surface (36) = {26};
Ruled Surface (37) = {27};
Ruled Surface (38) = {28};
Ruled Surface (39) = {29};
Ruled Surface (40) = {30};
Ruled Surface (41) = {31};
Ruled Surface (42) = {32};

Transfinite Surface {31:42};
Recombine Surface {31:42};

/* ---- Construct All Volumes ---- */

Surface Loop (1) = {1,11,31,32,33,34};
Surface Loop (2) = {2,12,35,36,37,38};
Surface Loop (3) = {5,15,31,35,39,40};
Surface Loop (4) = {4,14,32,36,40,41};
Surface Loop (5) = {6,16,33,37,41,42};
Surface Loop (6) = {3,13,34,38,42,39};

Volume (1) = {1};
Volume (2) = {2};
Volume (3) = {3};
Volume (4) = {4};
Volume (5) = {5};
Volume (6) = {6};
Transfinite Volume {1:6};

/* ---- Physical Names (Boundary Names) ---- */

Physical Volume ("FLUID") = {1:6};
Physical Surface ("OVERSET") = {11:16};
Physical Surface ("SPHERE") = {1:6};

Mesh.ElementOrder = 3;
