#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "raylib.h"

#define screenWidth 800
#define screenHeight 600
#define hexagonSize 5
#define Nx 190
#define Ny 130

const float sqrt3 = sqrtf(3.0f);
float T = 1.0f;
float J = 1.0f;
float MB = 0.0f;
size_t states = 3;
#define dT 0.1f
#define dMB 0.2f
#define dJ 0.1f

int random_spin()
{
    if (states%2 != 0) return (rand()%states) - (states-1)/2;
    else return 2*(rand()%states) - states/2 - states/2 + 1;
}

void initialize_spins(int spins[Nx][Ny])
{
    for (size_t i = 0; i < Nx; i++) {
        for (size_t j = 0; j < Ny; j++) {
            spins[i][j] = random_spin();
        }
    }
}

void draw_spins(int spins[Nx][Ny])
{
    for (size_t i = 0; i < Nx; i++) {
        for (size_t j = 0; j < Ny; j++) {
            float x = i*sqrt3*hexagonSize/2.0f;
            float y = (j)*hexagonSize - hexagonSize*sqrt3/2;
            if (i%2 == 1) y -= hexagonSize/2;
            Color color;
            switch (spins[i][j])
            {
            case -1:
                color = RED;
                break;
            case 1:
                color = BLUE;
                break;
            case 0:
                color = BLACK;
                break;
            default:
                break;
            }
            const float oversize = 1.2f; // So there's no padding between the hexagons
            DrawPoly((Vector2) {x, y}, 6, (float) hexagonSize/2.0f*oversize, 0.0f, color);
        }
    }
}

int boundary(int i, int N)
{
    if (i == -1) return N - 1;
    if (i == N) return 0;
    else return i;
}

float deltaE(int spins[Nx][Ny], int i, int j, int original_spin) {
    int neighbors = 0;
    neighbors += spins[boundary(i+1, Nx)][j] + spins[boundary(i-1, Nx)][j];
    neighbors += spins[boundary(i+1, Nx)][boundary(j-1, Ny)] + spins[boundary(i-1, Nx)][boundary(j-1, Ny)] + spins[i][boundary(j-1, Ny)];
    neighbors += spins[i][boundary(j+1, Ny)];
    int delta_spin = spins[i][j] - original_spin;
    return -delta_spin * J * neighbors - MB*delta_spin;
}

void step(int spins[Nx][Ny])
{   
    int i = rand()%Nx;
    int j = rand()%Ny;
    int original_spin = spins[i][j];
    
    spins[i][j] = random_spin();
    
    float delta = deltaE(spins, i, j, original_spin);
    
    if (delta < 0) return;

    float r = (float) rand()/RAND_MAX;
    
    if (r > exp(-delta/T)) spins[i][j] = original_spin;
}

#define N_avg 20
int main(void)
{
    int spins[Nx][Ny];
    size_t counter = 0;
    initialize_spins(spins);
    float E[N_avg], M[N_avg];
    float E_display, M_display;
    
    InitWindow(screenWidth, screenHeight, "Ising Model Simulation");

    SetTargetFPS(60);

    // Main game loop
    while (!WindowShouldClose())
    {
        counter += 1;
        // Update state of the spins
        #pragma omp parallel for
        for (size_t i = 0; i < 2000; i++) step(spins);

        if (IsKeyPressed(KEY_J)) {
            T -= dT;
        }
        if (IsKeyPressed(KEY_K)) {
            T += dT;
        }
        if (IsKeyPressed(KEY_M)) {
            MB += dMB;
        }
        if (IsKeyPressed(KEY_N)) {
            MB -= dMB;
        }
        if (IsKeyPressed(KEY_O)) {
            J += dJ;
        }
        if (IsKeyPressed(KEY_I)) {
            J -= dJ;
        }
        if (IsKeyPressed(KEY_TWO)) {
            states = 2;
            initialize_spins(spins);
        }
        if (IsKeyPressed(KEY_THREE)) {
            states = 3;
            initialize_spins(spins);
        }

        // Draw new state
        BeginDrawing();
            ClearBackground(BLACK);
            draw_spins(spins);
            DrawText(TextFormat("Temperature: %f", T), 30, 30, 20, RAYWHITE);
            DrawText(TextFormat("Magnetic field: %f", MB), 30, 60, 20, RAYWHITE);
            DrawText(TextFormat("J: %f", J), 30, 90, 20, RAYWHITE);
        EndDrawing();
    }

    CloseWindow();

    return 0;
}