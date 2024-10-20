#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h> 
#include <unistd.h> //for sleep function while executing automatic mode
#include "movement.h"
#include "functions.h"

//arleady defined in movement.c
extern int i,j,k;

int main(){
  printf("Hello, welcome to 2048 !\n");
  printf("What is your keyboard type ?\nEnter 1 for AZERTY, 2 for QWERTY or anything else if you want to customize your game experience : ");
  int best_score = 0;
  int keyboard = getchar();
  int upmov;
  int leftmov;
  int rightmov;
  int downmov;
  int undo;
  int menu_key;
  if (keyboard == '1'){ /*azerty*/
        upmov = 'z';
        leftmov = 'q';
        rightmov = 'd';
        downmov = 's';
        undo = 'v';
        menu_key = 'm';
  }
  else if (keyboard == '2'){ /*qwerty*/
        upmov = 'w';
        leftmov = 'a';
        rightmov = 'd';
        downmov = 's';
        undo = 'v';
  } 
  else { /*personal config*/
        printf("Your keyboard is not available yet, but you can choose your configuration right now:\n");
        printf("Menu button : ");
        menu_key = getchar();
        printf("\nLeft key : "); /*configuration of all the keys*/
        leftmov = getchar();
        printf(" Right key : ");
        rightmov = getchar();
        printf(" Up key : ");
        upmov = getchar();
        printf(" Down key : ");
        downmov = getchar();
        printf(" Undo key : ");
        undo = getchar();
  }

  printf("\nEnter the size of the grid you wanna play with :  ");
  int size = getchar() - 48;
  while(size<2){ //SIZE ENTRY ERRORS
    printf("\nthe size might be > 1 : ");
    size = getchar() - 48;
    if (size == -21){ //exit
        printf("\n Exit Menu, see you soon !\n");
        return 0;
        }
  }
  int** a = (int **)malloc((size+3) * sizeof(int *));

  for (i=0 ; i<size ; ++i) { // allocation de mémoire pour chaque ligne
      a[i] = (int *)malloc((size+3) * sizeof(int));
  }

  char n;
  char upmovc = upmov;  /*conversion of integers (ASCII codes) in characters to print them*/
  char downmovc = downmov;
  char leftmovc = leftmov;
  char rightmovc = rightmov;
  char undoc = undo;

  while(1){/*first while loop that allows to save the best score of the game*/
    int score = 0;
    printf("\nChoose difficulty : 1 is automatic, 5 is the hardest : ");
    int difficulty = getchar() - 49;
    while(difficulty!=1 && difficulty!=2 && difficulty!=3 && difficulty!=4 && difficulty!=0){ 
      // 0 for automatic mode, see menu for others : DIFFICULTY ERROR CASES
      printf("\npick up an integer between 1 and 4 to choose difficulty : ");
      difficulty = getchar() - 49;
      if (difficulty == -22){ //exit
        printf("\n Exit Menu, see you soon !\n");
        my_free(a,size);
        return 0;
        }
    }

    system("stty -icanon");
    ran(difficulty, size, a);
    
    int** temp = (int **)malloc((size+3)* sizeof(int *)); /*temporary board for the undo button*/
    for (i=0; i<size; ++i) {// allocation de mémoire pour chaque ligne
        temp[i] = (int *)malloc((size+3)* sizeof(int));
    }
    int win = 0;
    while(1){/*while loop for one game*/
      ran(difficulty, size, a);
      
      system("clear");//clear screen
      if (n == 27){//ESC to exit game
          print_menu(upmovc, downmovc, leftmovc, rightmovc, undoc);
          n = getchar();
          if(n == 27){
            printf("\n Exit Menu, see you soon !\n");
            my_free(a,size);
            my_free(temp,size);
            return 0;
          }
      }
      print_mode(difficulty);
      print_board(size, a);
      printf("your score : %5d        best score : %5d\n", score, best_score);
      printf("%c, %c, %c, %c, %c, %c --> up down left right undo Menu; ESC to end the game\n",upmovc, downmovc, leftmovc, rightmovc, undoc, menu_key);
      while(1){ /*the switch loop didn't allow to use variables already defined in the cases so we needed a while, 
        and we replaced "goto lab" by "continue"*/
        if (difficulty == 0){ //for automatic mode
            n = random_mov(upmov,downmov,leftmov,rightmov);
            sleep(1); //to make a pause while printing the board, otherwise it would be too quick
        }else{
            n = getchar();
        }
        score_calculation(&score, size,a);
        if (n == undo){
          copy_board(temp, a, size);
          break;
        }
        else if(n == upmov){ //up
          copy_board(a, temp, size);
          if(0==up(size, a)){
            if(1==fail(size, a)){
              printf("Press ESC to leave\n");
              sleep(2);
              for(i=0 ; i<size ; ++i){
                for(j=0 ; j<size ; ++j){
                  a[i][j] = 0;
                }
              }
              break;}
            else
              
              continue;
            }
            break;
        }
        else if(n == downmov){//down
          copy_board(a, temp, size);
          if(0==down(size,a)){
            if(1==fail(size, a)){
              printf("Press ESC to leave\n");
              sleep(2);
              for(i=0 ; i<size ; ++i){
                for(j=0 ; j<size ; ++j){
                  a[i][j] = 0;
                }
              }
              break;}
            else
              continue;
            }
            break;
        }
        else if(n == leftmov){//left
          copy_board(a, temp, size);
          if(0==left(difficulty, size, a)){
            if(1==fail(size, a)){
              printf("Press ESC to leave\n");
              sleep(2);
              for(i=0 ; i<size ; ++i){
                for(j=0 ; j<size ; ++j){
                  a[i][j] = 0;
                }
              }
              break;}
            else
              continue;
            }
            break;
        }
        else if(n == rightmov){//right
          copy_board(a,temp, size);
          if(0==right(size, a)){
            if(1==fail(size, a)){
              printf("Press ESC to leave\n");
              sleep(2);
              for(i=0 ; i<size ; ++i){
                for(j=0 ; j<size ; ++j){
                  a[i][j] = 0;
                }
              }
              break;}
            else
              continue;
          }
          break;
        }
        else if (n == 27){//ESC
          for(i=0 ; i<size ; ++i){
            for(j=0 ; j<size ; ++j){
              a[i][j] = 0;
            }
          }
          break;
        }
        else if (n == 'm'){
          print_menu(upmovc, downmovc, leftmovc, rightmovc, undoc);
          continue;
        }
        else{
          continue;
        }
      }
      //win or not with possibility to continue
      for(i=0;i<size;i++){
        for(j=0;j<size;j++){
          if(a[i][j] >= 2048 && win<3)//2048 win
          {
            win = win + 1;
          } 
        }
      }
      if(win == 1){
        printf("\nYou win！\npress ESC to end the game, anywhere else to continue\n");
        int try_again = getchar();
        if(try_again == 27){
          printf("Well Done !\n");
          break;
        }
      }
      if (best_score < score){
        best_score = score;
      }
    }
  }
  my_free(a,size);
  return 0;
}
