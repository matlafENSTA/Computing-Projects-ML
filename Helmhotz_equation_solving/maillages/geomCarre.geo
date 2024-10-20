Mesh.MshFileVersion = 2.2;

// Definition du pas du maillage
h = 0.01;

// Definition des points (en 3D, raison pour laquelle il y a un 0 en z)
Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {1, 1, 0, h};
Point(4) = {0, 1, 0, h};

// Definition des segments qui relient les points
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Definition des contours fermes
Line Loop(1) = {1,2,3,4};

// Definition des surfaces a partir contours fermes
Plane Surface(1) = {1};

// Definition des elements physiques : pour ces elements, nous pourrons recuperer
//	les references
/***** Dirichlet sur tout le bord *****/
/* Decommenter les lignes d'apres pour imposer Dirichlet sur le bord */
Physical Point(1) = {1,2,3,4};
Physical Line(1) = {1,2,3,4};
Physical Surface(1) = {1};
/**************************************/

/****** Fourier sur tout le bord ******/
/* Decommenter les lignes d'apres pour imposer Fourier sur le bord */
// Physical Point(1) = {1,2,3,4};
// Physical Line(2) = {1,2,3,4};
// Physical Surface(1) = {1};
/**************************************/
