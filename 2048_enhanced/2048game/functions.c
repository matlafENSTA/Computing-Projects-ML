#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "movement.h"
//#include "include/SDL.h"
//#include <SDL2/SDL_ttf.h> 

extern int i,j,k; //arleady defined in movement.c

/*
SDL_Window* window ;
SDL_Renderer* renderer ;
SDL_Color CELL_Colors ;
SDl_Color textColor ;
TTF* font ;

void print_board(int** grid, int size){ //SDL VERSION
    int cell_size = 50 ;
    int edges = 5 ;
    for(i=0;i<size;i++){
        for (j=0;j<size;j++){
            int x = j*(cell_size + edges); 
            int y = i*(cell_size + edges);
            int value = grid[i][j] ;   
                if (value == 0) {
                    color_index = 0 ;  
                }
                else {
                    color_index = log2(value) ;
                }

            // Draw the cell 

            sprintf(value_char, "%d", value) ; 

            SDL_Color cellColor = CELL_Colors[color_index] ;
 
            SDL_Surface* textSurface = TTF_RenderText_Solid(font, value_char, textColor) ; 
            SDL_Texture* textTexture = SDL_CreateTextureFromSurface(renderer, textSurface) ;
             
            SDL_Rect cellRect = { x , y , cell_size , cell_size };

            SDL_SetRenderDrawColor(renderer, cellColor.r, cellColor.g, cellColor.b, 255); 
            SDL_RenderFillRect(renderer, &cellRect);

            // Centering the number inside the cell 

            SDL_Rect text_position = {
                cellRect.x + (cellRect.w - textSurface->w) / 2, 
                cellRect.y + (cellRect.h - textSurface->h) / 2, 
                textSurface->w, 
                textSurface->h 
            }; 

            SDL_RenderCopy(renderer, textTexture, NULL, &text_position); 
        }
    }
    SDL_RenderPresent(renderer) ;
}*/

void colored_cell(int n, int** a, int size){//SHELL VERSION
  if(n == 3){
    printf("\033[7;43;37m%4d   \033[0m|",a[i][j]);
  }
  if(n == 2){
    printf("\033[7;49;39m%4d   \033[0m|",a[i][j]);
  }if(n == 4){
    printf("\033[1;46;37m%4d   \033[0m|",a[i][j]);
  }if(n == 8){
    printf("\033[1;44;37m%4d   \033[0m|",a[i][j]);
  }if(n == 16){
    printf("\033[1;42;37m%4d   \033[0m|",a[i][j]);
  }if(n == 32){
    printf("\033[1;43;37m%4d   \033[0m|",a[i][j]);
  }if(n == 64){
    printf("\033[1;41;36m%4d   \033[0m|",a[i][j]);
  }if(n == 128){
    printf("\033[1;45;95m %4d  \033[0m|",a[i][j]);
  }if(n == 256){
    printf("\033[1;44;32m %4d  \033[0m|",a[i][j]);
  }if(n == 512){
    printf("\033[1;46;31m %4d  \033[0m|",a[i][j]);
  }if(n == 1024){
    printf("\033[1;40;95m  %4d \033[0m|",a[i][j]);
  }if(n == 2048){
    printf("\033[1;47;95m  %4d \033[0m|",a[i][j]);
  }if(n > 2048){
    printf("  %4d ",a[i][j]);
  }
}

void print_board(int size, int** a){ /*just to clear the main*/
  for(k=0 ; k<size ; ++k){
    printf("--------");
  }
  printf("-\n");
  for(i=0;i<size;i++)
    {
    printf("|");
    for(j=0;j<size;j++)
    {
      if(a[i][j]==0)
      printf("       |");
      else
      colored_cell(a[i][j],a,size);
    }
    printf("\n");
    for(k=0 ; k<size ; ++k){
      printf("--------");
    }
    printf("-\n");
  }
}

void my_free(int** grid, int size){ //to free the malloc
  for (i=0 ; i<size ; ++i){
    free(grid[i]);
  }
  free(grid);
}

void copy_board(int** original, int** new, int size){ /*copies the numbers from a grid to another, for the undo button*/
  for(i=0 ; i<size ; ++i){
    for(j=0 ; j<size ; ++j){
      new[i][j] = original[i][j];
    }
  }
}

void score_calculation(int* score, int size, int** a){
  *score = 0;
  for(int i=0 ; i<size ; ++i){
    for(int j=0 ; j<size ; ++j){
      *score = *score + a[i][j];
    }
  }
}

void print_mode(int difficulty){//just to print the mode you're playing in the shell
  switch(difficulty){
    case 0:
      printf("\n         automatic mode");
      break;
    case 1:
      printf("\n            easy mode");
      break;
    case 2:
      printf("\n           regular mode");
      break;
    case 3:
      printf("\n            hard mode");
      break;
    case 4:
      printf("\n           insane mode");
      break;
    default:
      printf("\n           unknown mode");
      break;
  }
  printf("\n");
}

void print_menu(char upmovc, char downmovc, char leftmovc, char rightmovc, char undoc){//print the menu
  system("clear");
  printf(" ~ Menu ~\n- Difficulty\n  1 - automatic mode : the computer is playing lonely using the corner technique\n");
  printf("  2 - easy : high probability to get a 4\n");
  printf("  3 - regular\n  4 - hard\n");
  printf("  5 - insane : when moving left, you divide by 2 instead of multiplying\n");
  printf("- Controls\n  %c - up\n  %c - down\n  %c - left\n  %c - right\n  %c - undo\n  ESC - end the game\n",upmovc, downmovc, leftmovc, rightmovc, undoc);
  printf ("\033[32m press ESC to exit Menu, anywhere else to try again\n \033[0m|");
}

int random_mov(int upmov,int downmov,int leftmov,int rightmov){ //pilot for automatic mode
  srand(time(NULL));
  int r = rand()%10;
  printf("%d\n",r);
  if(r < 1){
    return downmov;
  }else if(r < 2){
    return leftmov;
  }else if(r < 6){
    return upmov;
  }else{
    return rightmov;
  }
}