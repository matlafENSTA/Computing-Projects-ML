/*****************************************************************************/
/* Generation succesive de maillages associes a une geometrie Pacman.        */
/*****************************************************************************/
Mesh.MshFileVersion = 2.2;

/* Déclaration de variables */
DefineConstant[
  R = {1, Name "Rayon"},                         /* Rayon de l'arc exterieur */
  h = {0.5, Name "Pas maillage"}                     /* Pas maillage initial */
  theta = {30, Name "Angle"}                           /* Angle de la bouche */
  nBmaillages = {4, Name "Nombre maillages"}    /* Nombre de maillages (hormis le dernier) */
  prefixeMaillage = {"geomPacman", Name "prefixe maillage"}
];

/*****************/
/* Arc extérieur */
/*****************/
xc = 0;   /* Coordonnées... */
yc = 0;   /* ... du centre  */

/* Points lies a la geometrie */
Point(1) = {xc, yc, 0, h};
Point(2) = {xc + R*Cos(Pi*theta/180), yc + R*Sin(Pi*theta/180), 0, h};
Point(3) = {xc - R, yc, 0, h};
Point(4) = {xc + R*Cos(Pi*theta/180), yc - R*Sin(Pi*theta/180), 0, h};

/* Segments et arcs de cercle */
Line(1) = {1, 2};
Circle(2) = {2,1,3};
Circle(3) = {3,1,4};
Line(4) = {4, 1};

/********/
/* Oeil */
/********/
xo = xc + 0.5*R*Cos(Pi*theta/180) - 0.25*Sqrt(3)*R*Sin(Pi*theta/180);   /* Coordonnées... */
yo = yc + 0.5*R*Sin(Pi*theta/180) + 0.25*Sqrt(3)*R*Cos(Pi*theta/180);   /* ... du centre  */
Ro = 0.125*R;  /* Rayon */

/* Points liés à la géométrie */
Point(5) = {xo, yo, 0, h};
Point(6) = {xo + Ro, yo, 0, h};
Point(7) = {xo, yo + Ro, 0, h};
Point(8) = {xo - Ro, yo, 0, h};
Point(9) = {xo, yo - Ro, 0, h};

/* Segments et arcs de cercle */
Circle(5) = {6,5,7};
Circle(6) = {7,5,8};
Circle(7) = {8,5,9};
Circle(8) = {9,5,6};

/************************************************************/
/* Definition des courbes fermees et des surfaces associees */
/************************************************************/
Curve Loop(1) = {1,2,3,4};      /* Arc extérieur */
Curve Loop(2) = {5,6,7,8};      /* Oeil */
Plane Surface(1) = {1,2};       /* Surface */

/****************************************************/
/* Elements physiques pour recuperer les references */
/****************************************************/
Physical Point(1) = {1,2,4};
Physical Curve(1) = {1,4,5,6,7,8};  /* Bouche et oeil */
Physical Curve(2) = {2,3};           /* Arc exterieur */
Physical Surface(1) = {1,2};             /* Interieur */

/**************************************************************************/
/* Generation de plusieurs maillages obtenus en raffinant successivement  */
/*                                                                        */
/* NOTE : le dernier maillage est celui sur lequel on calcule la solution */
/*        de reference. Il est obtenu en raffinant le precedent deux fois */
/*        de suite.                                                       */
/**************************************************************************/
/* On genere le premier maillage qu'on enregistre */
Mesh 2;
Save Sprintf(StrCat(prefixeMaillage, "-%03g-%01g.msh"), theta, 1);

/* Le premier maillage est ensuite raffine pour obtenir ceux d'apres */
For i In {2:nBmaillages}
  RefineMesh;
  Save Sprintf(StrCat(prefixeMaillage, "-%03g-%01g.msh"), theta, i);
EndFor

/* Le dernier maillage est obtenu en raffinant deux fois */
RefineMesh;
RefineMesh;
Save Sprintf(StrCat(prefixeMaillage, "-%03g-%01g.msh"), theta, 0);
