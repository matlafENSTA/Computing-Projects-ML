# IN104-2048

IMPORTANT :

This version of 2048 allows you to play using the shell and an SDL graphic interface. The game is by default in the shell because you need to install SDL on your computer to play with it. If you wanna try, you can follow the following steps in the "functions.c" file: 
- put under comment the function "colored shell" SHELL version and use the SDL version instead.
- remove the // in front of #include "include/SDL.h" and #include <SDL2/SDL_ttf.h>
- build the executable with the command "make"

2048 RULES :

- The aim : reach the value 2048 in a cell.
- The board : a square grid (you can choose the size at the beginning)
- Use the keyboard to move all the cells to the corresponding side. When two cells with the same number
collide, - they merge in a cell with the number equal to the double of the number in the original cells. Each time a move is made, a cell with a 2 or a 4 randomly appears in the grid.  

HOW TO PLAY :
- in the shell, enter make clean to build the game executable. It's named 2048.x
- enter ./2048.x to start the game, then follow the instructions (you can also enter make run to play instantly)
- you will play using your keyboard and see the grid in the shell : you can
  choose the keys you wanna play with. Don't forget that you can leave the 
  game at any moment using ESC.
- you will be ask to choose your keyboard configuration between AZERTY, QWERTY or with the arrows. 
- then you will be asked to choose the game diffulty : 
    1 : the game is played randomly 
    2 to 4 : the higher the number the lower chance to generate a 4 
    5 : going left will divide all the numbers by 2  
- you can now play using the commands you selected and try to reach the higher score ! 
- if you made a move you did not want to do, you can press v to undo and return to the previous configuration.
 
technical specifications sheet : 
- 2048-template.c is the main file
- movement.c is a file that contains all the basic functions to run the game
- functions.c contains additional code to run new functionalities