#ifndef MAIN_H
#define MAIN_H

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <GL/glew.h>
#include <SDL2/SDL.h>
#include "structs.h"

//screen properties
const int SCREEN_FPS = 60;
#define screenWidth 640
#define screenHeight 480
#define viewWidth 10.0f
#define viewHeight (screenHeight*viewWidth/screenWidth)

// general Particle properties
#define yield 0.08f
#define stiffness 0.08f
#define nearStiffness 0.1f
#define linearViscocity 0.5f
#define quadraticViscocity 1.f
#define maxParticles 6000

//Particle size
#define radius 0.05f
#define particleHeight (6*radius)

//Particle rendering properties
#define drawRatio 3.0f
#define velocityFactor 0.1f

#define frameRate 20
#define timeStep 7
#define Pi 3.14159265f
#define deltaTime ((1.0f/frameRate) / timeStep)

#define restDensity 75.0f
#define surfaceTension 0.0006f
#define multipleParticleTypes true

SDL_Window* screen = NULL;
SDL_GLContext gContext;

//world boundaries
Boundary boundaries[4] =
{
	Boundary(1, 0, 0),
	Boundary(0, 1, 0),
	Boundary(-1, 0, -viewWidth),
	Boundary(0, -1, -viewHeight)
};

//types of materials, rgb orders lightest to heaviest, and the sources that will use them
#define red ParticleType(Vec4(0.75f, 0.1f, 0.1f, 1), 1.0f)
#define green ParticleType(Vec4(0.1f, 0.75f, 0.1f, 1), 1.2f)
#define blue ParticleType(Vec4(0.1f, 0.1f, 0.75f, 1), 1.4f)
Source sources[3] =
{
	Source(red, Vec2(0.05f*viewWidth, 0.8f*viewHeight), Vec2(4, 1), 5),
	Source(blue, Vec2(viewWidth - 0.05f*viewWidth, 0.8f*viewHeight), Vec2(-4, 1), 5),
	Source(green, Vec2(0, 0), Vec2(0, -1), 7),
};
bool activeSpout = false;

//Particle references for different functions
int totalParticles = 0;
Particle particles[maxParticles];
Springs neighbours[maxParticles];
Vec2 prevPos[maxParticles];
ParticleType particleTypes[maxParticles];
Vec4 savedParticleColors[maxParticles];

//the map of the scene
const int mapWidth = (int)(viewWidth / particleHeight);
const int mapHeight = (int)(viewHeight / particleHeight);
Particle* map[mapWidth][mapHeight];
int mapCoords[maxParticles][2];

Vec2 gravity(0.f, -9.8f);

void handleEvents(SDL_Event e);

#endif