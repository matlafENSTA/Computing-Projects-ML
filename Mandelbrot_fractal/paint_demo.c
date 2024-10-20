#include <stdio.h>

#include "paint.h"

int main ()
{
  int largeur = 60;
  int hauteur = 40;

  // Creer les pixels pour dessiner
  unsigned char* pixels = create_pixels(largeur,hauteur);

  unsigned char r,g,b;

  // couleur rouge
  r = 255;
  g = 0;
  b = 0;
  // Colorier le pixel (5,10) à (9,10)
  color_pixel(5,10,r,g,b,largeur,hauteur,pixels);
  color_pixel(6,10,r,g,b,largeur,hauteur,pixels);
  color_pixel(7,10,r,g,b,largeur,hauteur,pixels);
  color_pixel(8,10,r,g,b,largeur,hauteur,pixels);
  color_pixel(9,10,r,g,b,largeur,hauteur,pixels);

  /* Ecrire vers un fichier BMP */
  save_BMP("paint_demo.bmp", largeur, hauteur, pixels);

  /* Détruire le canvas pour libérer la mémoire */
  destroy_pixels(pixels);

  return (0) ;
}
