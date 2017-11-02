R1 = .06;  // Size of inner cube
R2 = .25;  // Length of outer prism
R3 = .25;  // Width of outer prism
NN = 5;  // Number of cells in each direction on inner cube
NL = 6;  // Number of cells going outwards
prog1 = 1.25;  // Geometric progression factor for outward layers

/* ---- Inner Sphere Surface ---- */

Point(1) = {0.0,0.0,0.0,1.};

Point(2) = {-R1, -R1, -R1};
Point(3) = { R1, -R1, -R1};
Point(4) = { R1,  R1, -R1};
Point(5) = {-R1,  R1, -R1};

Point(6) = {-R1, -R1, R1};
Point(7) = { R1, -R1, R1};
Point(8) = { R1,  R1, R1};
Point(9) = {-R1,  R1, R1};

// 'Bottom' circles
Line(1) = {2,3};
Line(2) = {3,4};
Line(3) = {4,5};
Line(4) = {5,2};

// 'Top' circles
Line(5) = {6,7};
Line(6) = {7,8};
Line(7) = {8,9};
Line(8) = {9,6};

// 'Vertical' circles
Line(9) = {2,6};
Line(10) = {3,7};
Line(11) = {4,8};
Line(12) = {5,9};

Transfinite Line {1:12} = NN+1 Using Progression 1.0;

Line Loop (1) = {1,2,3,4};
Line Loop (2) = {5,6,7,8};
Line Loop (3) = {4,9,-8,-12};
Line Loop (4) = {2,11,-6,-10};
Line Loop (5) = {1,10,-5,-9};
Line Loop (6) = {3,12,-7,-11};

Plane Surface (1) = {1};
Plane Surface (2) = {2};
Plane Surface (3) = {3};
Plane Surface (4) = {4};
Plane Surface (5) = {5};
Plane Surface (6) = {6};

Transfinite Surface {1:6};
Recombine Surface {1:6};


/* ---- Outer Sphere Surface ---- */

Point(11) = {0.0,0.0,0.0,1.};

Point(12) = {-R2, -R3, -R2};
Point(13) = { R2, -R3, -R2};
Point(14) = { R2,  R3, -R2};
Point(15) = {-R2,  R3, -R2};

Point(16) = {-R2, -R3, R2};
Point(17) = { R2, -R3, R2};
Point(18) = { R2,  R3, R2};
Point(19) = {-R2,  R3, R2};

// 'Bottom' circles
Line(21) = {12,13};
Line(22) = {13,14};
Line(23) = {14,15};
Line(24) = {15,12};

// 'Top' circles
Line(25) = {16,17};
Line(26) = {17,18};
Line(27) = {18,19};
Line(28) = {19,16};

// 'Vertical' circles
Line(29) = {12,16};
Line(30) = {13,17};
Line(31) = {14,18};
Line(32) = {15,19};

Transfinite Line {21:32} = NN+1 Using Progression 1.0;

Line Loop (11) = {21,22,23,24};
Line Loop (12) = {25,26,27,28};
Line Loop (13) = {24,29,-28,-32};
Line Loop (14) = {22,31,-26,-30};
Line Loop (15) = {21,30,-25,-29};
Line Loop (16) = {23,32,-27,-31};

Plane Surface (11) = {11};
Plane Surface (12) = {12};
Plane Surface (13) = {13};
Plane Surface (14) = {14};
Plane Surface (15) = {15};
Plane Surface (16) = {16};

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

Plane Surface (31) = {21};
Plane Surface (32) = {22};
Plane Surface (33) = {23};
Plane Surface (34) = {24};
Plane Surface (35) = {25};
Plane Surface (36) = {26};
Plane Surface (37) = {27};
Plane Surface (38) = {28};
Plane Surface (39) = {29};
Plane Surface (40) = {30};
Plane Surface (41) = {31};
Plane Surface (42) = {32};

Transfinite Surface {31:42};
Recombine Surface {31:42};

/* ---- Construct All Volumes ---- */

Surface Loop (1) = {1,11,31,32,33,34};
Surface Loop (2) = {2,12,35,36,37,38};
Surface Loop (3) = {5,15,31,35,39,40};
Surface Loop (4) = {4,14,32,36,40,41};
Surface Loop (5) = {6,16,33,37,41,42};
Surface Loop (6) = {3,13,34,38,42,39};
Surface Loop (7) = {1,2,3,4,5,6};

Volume (1) = {1};
Volume (2) = {2};
Volume (3) = {3};
Volume (4) = {4};
Volume (5) = {5};
Volume (6) = {6};
Volume (7) = {7};
Transfinite Volume {1:7};

/* ---- Physical Names (Boundary Names) ---- */

Physical Volume ("FLUID") = {1:7};
Physical Surface ("CHAR") = {11:16};

Mesh.ElementOrder = 1;
