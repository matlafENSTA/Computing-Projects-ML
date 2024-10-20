Mesh.MshFileVersion = 2.2;

// Definition du pas du maillage
h = 0.1;

// Definition des points (en 3D, raison pour laquelle il y a un 0 en z)
Point(1) = {0, 0, 0, h};
Point(2) = {4, 0, 0, h};
Point(3) = {4, 2, 0, h};
Point(4) = {0, 2, 0, h};
Point(5) = {4, 1.5, 0, h};
Point(6) = {0., 0.5, 0, h};

// Definition des segments qui relient les points
Line(1) = {1, 2};
Line(3) = {3, 4};
Line(5) = {2, 5};
Line(6) = {5, 6};
Line(7) = {6, 1};
Line(8) = {5, 3};
Line(9) = {4, 6};

// Definition des contours fermes
Line Loop(1) = {1,5,6,7};
Line Loop(2) = {-6,8,3,9};

// Definition des surfaces a partir contours fermes
Plane Surface(1) = {1};
Plane Surface(2) = {2};

// Definition des elements physiques : pour ces elements, nous pourrons recuperer
//								les references
Physical Point(1) = {1,2,3,4};
Physical Surface(1) = {1};
Physical Surface(2) = {2};

/* Ne decommenter qu'UNE des deux lignes ci-dessous pour imposer soit des   */
/* conditions de Dirichlet, soit des conditions de Fourier                  */
Physical Line(1) = {1,5,8,3,9,7};             /* Dirichlet sur tout le bord */
// Physical Line(2) = {1,5,8,3,9,7};            /* Fourier sur tout le bord */
