Mesh.MshFileVersion = 2.2;
Mesh.CharacteristicLengthMin = 0.02;
Mesh.CharacteristicLengthFactor = 1;
Mesh.Optimize = 1;

Point(1) = {0,0,0};
Point(2) = {1,0,0};
Point(3) = {0,1,0};
Point(4) = {0,-1,0};
Point(5) = {-1,0,0};
Circle(1) = {2,1,3};
Circle(2) = {3,1,5};
Circle(3) = {5,1,4};
Circle(4) = {4,1,2};
Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};
Extrude {{0,1,0}, {-2,0,0}, 2*Pi/3} {Surface {6}; }
Extrude {{0,1,0}, {-2,0,0}, 2*Pi/3} {Surface {28}; }
Extrude {{0,1,0}, {-2,0,0}, 2*Pi/3} {Surface {50}; }

Physical Volume(1) = {1};
Physical Volume(2) = {2};
Physical Volume(3) = {3};
Physical Surface(11) = {71,67,63,59};
Physical Surface(12) = {49,45,41,37};
Physical Surface(13) = {27,23,19,15};
Physical Surface(14) = {50,28,6};
