#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

int i,j,k;


int move_up(int size, int** a)
{
 int f=0;//have moved then f=1，no movement is 0
 for(j=0;j<size;j++)
 {
  for(i=1;i<size;i++)
  {
   if(a[i][j]!=0)//find the number not 0
   {
    for(k=0;k<i;k++)
    {
     if(a[k][j]==0)//the first 0 upper
     {
      a[k][j]=a[i][j];
      a[i][j]=0;
      f=1;
      break;
     }
    }
   }
  }
 }
 return f;
}
 
int move_down(int size, int** a)
{
 int f=0;
 for(j=0;j<size;j++)
 {
  for (i=0 ; i<size-1 ; i++)
  {
   if(a[i][j]!=0)
   {
    for(k=size-1;k>i;k--)
    {
     if(a[k][j]==0)
     {
      a[k][j]=a[i][j];
      a[i][j]=0;
      f=1;
      break;
     }
    }
   }
  }
 }
 return f;
}
 
int move_left(int size, int** a)
{
 int f=0;
 for(i=0;i<size;i++)
 {
  for(j=1;j<size;j++)
  {
   if(a[i][j]!=0)
   {
    for(k=0;k<j;k++)
    {
     if(a[i][k]==0)
     {
      a[i][k]=a[i][j];
      a[i][j]=0;
      f=1;
      break;
     }
    }
   }
  }
 }
 return f;
}

int move_right(int size, int** a)
{
 int f=0;
 for(i=0;i<size;i++)
 {
  for (j = 0; j < size - 1; j++)
  {
   if(a[i][j]!=0)
   {
    for(k=size-1;k>j;k--)
    {
     if(a[i][k]==0)
     {
      a[i][k]=a[i][j];
      a[i][j]=0;
      f=1;
      break;
     }
    }
   }
  }
 }
 return f;
}

int up(int size, int** a)
{
 int f=0;//have merged，f=1，no :0
 //move
 int f1=move_up(size,a);//have moved f1=1
 //merge
 for(j=0;j<size;j++)
 {
  for(i=0;i<size-1;i++)
  {
   if(a[i][j]==a[i+1][j] && a[i][j]!=0)
   {
    a[i][j]=2 * a[i][j];
    a[i+1][j]=0;
    f=1;
   }
  }
 }
 if(f==1)
  move_up(size,a);
 return (f||f1);
}
 
int down(int size, int** a)
{
 int f=0;
 //move
 int f1=move_down(size,a);
 //merge
 for(j=0;j<size;j++)
 {
  for(i=size-1;i>0;i--)
  {
   if(a[i][j]==a[i-1][j] && a[i][j]!=0)
   {
    a[i][j] = 2*a[i][j];
    a[i-1][j]=0;
    f=1;
   }
  }
 }
 if(f==1)
  move_down(size,a);
 return (f||f1);
}
 
int left(int difficulty, int size, int** a)
{
 int f=0;
 //move
 int f1=move_left(size,a);
 //merge
 for(i=0;i<size;i++)
 {
  for(j=0;j<size-1;j++)
  {
   if(a[i][j]==a[i][j+1]&&a[i][j]!=0)
   {
    if(a[i][j] != 2 && difficulty == 4){
      a[i][j]=a[i][j]/2;
    }
    else{
      a[i][j]=a[i][j]*2;
    }
    a[i][j+1]=0;
    f=1;
   }
  }
 }
 if(f==1)
  move_left(size,a);
 return (f||f1);
}

int right(int size, int** a)
{
 int f=0;
 //move
 int f1=move_right(size,a);
 //merge
 for(i=size-1;i>=0;i--)
 {
  for(j=size-1;j>0;j--)
  {
   if(a[i][j]==a[i][j-1] && a[i][j]!=0)
   {
    a[i][j]=2*a[i][j];
    a[i][j-1]=0;
    f=1;
   }
  }
 }
 if(f==1)
  move_right(size,a);
 return (f||f1);
}

void ran(int difficulty, int size, int** a){
    int z[8];
    switch(difficulty){
      case 0:
        z[0] = 2;
        z[1] = 2;
        z[2] = 2;
        z[3] = 2;
        z[4] = 2;
        z[5] = 2;
        z[6] = 2;
        z[7] = 4;
        break;
      case 1:
        z[0] = 2;
        z[1] = 2;
        z[2] = 2;
        z[3] = 2;
        z[4] = 2;
        z[5] = 2;
        z[6] = 4;
        z[7] = 4; //2 is much more possible
        break;
      case 2: //normal mode
        z[0] = 2;
        z[1] = 2;
        z[2] = 2;
        z[3] = 2;
        z[4] = 2;
        z[5] = 2;
        z[6] = 2;
        z[7] = 4;
        break;
      case 3:
        z[0] = 2;
        z[1] = 2;
        z[2] = 2;
        z[3] = 2;
        z[4] = 2;
        z[5] = 2;
        z[6] = 2;
        z[7] = 2;
        break;
      case 4:
        z[0] = 2;
        z[1] = 2;
        z[2] = 2;
        z[3] = 2;
        z[4] = 2;
        z[5] = 2;
        z[6] = 2;
        z[7] = 4;
        break;
    }
  srand(time(NULL));
  lb:
  i=rand()%size;
  j=rand()%size;
  if(a[i][j]==0)
    a[i][j]=z[rand()%8];
  else
    goto lb;
}

int fail(int size, int** a)//fail
{
 int count=0;
 for(i=0;i<size;i++)
 {
  for(j=0;j<size-1;j++)
  {
   if(a[i][j]==a[i][j+1])//net yet fail
   {
    return 0;
   }
  }
 }
 for(i=0;i<size-1;i++)
 {
  for(j=0;j<size;j++)
  {
   if(a[i][j]==a[i+1][j])//net yet fail
   {
    return 0;
   }
  }
 }
 for(i=0;i<size;i++)
 {
  for(j=0;j<size;j++)
  {
   if(a[i][j]==0)//count the blank squares
   {
    count++;
   }
  }
 }
 if(count==0)//no blank
 {
  printf("\nYou nearly win, unfortunately!\n");
  return 1;
 }
 return 0;
}