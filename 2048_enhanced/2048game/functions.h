void colored_cell(int n);

void print_board(int size, int** a);

void my_free(int** grid, int size);

void copy_board(int** original, int** new, int size);

void score_calculation(int* score, int size, int** a);

void print_mode(int difficulty);

void print_menu(char upmovc, char downmovc, char leftmovc, char rightmovc, char undoc);

int random_mov(int upmov,int downmov,int leftmov,int rightmov);
