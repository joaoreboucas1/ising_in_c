set arg1=%1
gcc -o %arg1% %arg1%.c -IC:\raylib\raylib\src\ -LC:\raylib\raylib\src\ -lraylib -lgdi32 -lwinmm -fopenmp