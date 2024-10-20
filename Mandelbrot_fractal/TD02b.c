/* Mathieu LAFITTE
Q1*/
#include "paint.c"
#include <math.h>
#include <stdio.h>

int max_iter = 100;
int compter_iterations(double x, double y, int max_iter) {
  double a = 0;
  double b = 0;
  double c = 0;
  int k = 0;
  while (((a * a + b * b) <= 4) && (k < max_iter)) {
    c = a * b;
    a = a * a - b * b + x;
    b = 2. * c + y;
    k = k + 1;
  }
  return (k);
}
/*Q2
int main() {
  int u = compter_iterations(-1, 0, 100);
  printf("%d\n", u);
  return 0;
}
/*Q3
3
Q4
100
Q5
34
Q6*/
double fx(double x_left, double x_width, double width, int i) {
  double a = x_width / (width);
  double b = x_left;
  return a * i + b;
}

double fy(double y_top, double x_width, double width, double height, int j) {
  double y_height = x_width * height / width;
  double a = -y_height / height;
  double b = y_top;
  return a * j + b;
}
/*Q7
int main() {
  double width = 600;
  double height = 400;
  double x_left = -2;
  double y_top = 1;
  double x_width = 3;
  unsigned char *pixels = create_pixels(width, height);
  unsigned char r, g, b;
  for (int i = 1; i <= width; i = i + 1) {
    for (int j = 0; j <= height; j = j + 1) {
      double x = fx(x_left, x_width, width, i);
      double y = fy(y_top, x_width, width, height, j);

      int nb_iterations = compter_iterations(x, y, max_iter);
      int r = 255 * nb_iterations / max_iter;
      int g = 0;
      int b = 255 - r;
      if (nb_iterations == max_iter) {
        color_pixel(i, j, 0, 0, 0, width, height, pixels);
      } else {
        color_pixel(i, j, r, g, b, width, height, pixels);
      }
    }
  }

  save_BMP("Mandelbrot.bmp", width, height, pixels);

  destroy_pixels(pixels);

  return (0);
}
Q8*/
int main() {
  double width = 600;
  double height = 400;
  double x_left = -2;
  double y_top = 1;
  double x_width = 3;
  unsigned char *pixels = create_pixels(width, height);
  unsigned char r, g, b;
  for (int i = 1; i <= width; i = i + 1) {
    for (int j = 0; j <= height; j = j + 1) {
      double x = fx(x_left, x_width, width, i);
      double y = fy(y_top, x_width, width, height, j);

      int nb_iterations = compter_iterations(x, y, max_iter);
      int r = 150 + floor(sqrt(255 * nb_iterations / max_iter));
      int g = 100 + floor(50 * cos(x) + 50 * sin(y));
      int b = floor(255 - 100 * cos(r));
      if (nb_iterations == max_iter) {
        int rr = floor(0.01 * tan(x));
        int gg = floor(tan(y));
        int bb = floor(150 + 99 * cos(i + j));
        color_pixel(i, j, rr, gg, bb, width, height, pixels);
      } else {
        color_pixel(i, j, r, g, b, width, height, pixels);
      }
    }
  }
  save_BMP("Mandelbrot_colorisé.bmp", width, height, pixels);
  destroy_pixels(pixels);

  return (0);
}
/*Q9
int main() {
  double width = 600;
  double height = 400;
  double x_left = -2;
  double y_top = 1;
  double x_width = 3;
  double dx = 0.01;
  double dy = 0.02;
  unsigned char *pixels = create_pixels(width, height);
  unsigned char r, g, b;
  for (int i = 1; i <= width; i = i + 1) {
    for (int j = 0; j <= height; j = j + 1) {
      double x = fx(x_left, x_width, width, i);
      double y = fy(y_top, x_width, width, height, j);
      x = x + dx;
      y = y + dy;
      int nb_iterations = compter_iterations(x, y, max_iter);
      int r = 150 + floor(sqrt(255 * nb_iterations / max_iter));
      int g = 100 + floor(50 * cos(x) + 50 * sin(y));
      int b = floor(255 - 100 * cos(r));
      if (nb_iterations == max_iter) {
        int rr = floor(0.01 * tan(x));
        int gg = floor(tan(y));
        int bb = floor(150 + 99 * cos(i + j));
        color_pixel(i, j, rr, gg, bb, width, height, pixels);
      } else {
        color_pixel(i, j, r, g, b, width, height, pixels);
      }
    }
  }
  save_BMP("Mandelbrot_translaté.bmp", width, height, pixels);
  destroy_pixels(pixels);

  return (0);
}
/*Q10
int main() {
  double width = 600;
  double height = 400;
  double x_left = -2;
  double y_top = 1;
  double x_width = 3;
  double u = -0.1;
  unsigned char *pixels = create_pixels(width, height);
  unsigned char r, g, b;
  for (int i = 1; i <= width; i = i + 1) {
    for (int j = 0; j <= height; j = j + 1) {
      double x = fx(x_left, x_width, width, i);
      double y = fy(y_top, x_width, width, height, j);
      double xr = x * cos(u) + y * sin(u);
      y = y * cos(u) - x * sin(u);
      x = xr;
      int nb_iterations = compter_iterations(x, y, max_iter);
      int r = 150 + floor(sqrt(255 * nb_iterations / max_iter));
      int g = 100 + floor(50 * cos(x) + 50 * sin(y));
      int b = floor(255 - 100 * cos(r));
      if (nb_iterations == max_iter) {
        int rr = floor(0.01 * tan(x));
        int gg = floor(tan(y));
        int bb = floor(150 + 99 * cos(i + j));
        color_pixel(i, j, rr, gg, bb, width, height, pixels);
      } else {
        color_pixel(i, j, r, g, b, width, height, pixels);
      }
    }
  }
  save_BMP("Mandelbrot_rotation.bmp", width, height, pixels);
  destroy_pixels(pixels);

  return (0);
}
*/
/*
si on fait directement x = x*cos(g)+y*sin(g),
la variable x prend directement sa nouvelle valeur,
alors qu'on a encore besoin de l'ancienne pour définir y
*/