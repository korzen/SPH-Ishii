#include "SPH.h"

//Resets the map of the scene, re-adding every particle based on where it is at this moment
void updateMap()
{
	memset(map, 0, mapWidth*mapHeight*sizeof(Particle*));

	for (int i = 0; i<totalParticles; i++)
	{
		Particle& p = particles[i];
		int x = p.posX / particleHeight;
		int y = p.posY / particleHeight;

		if (x < 1) x = 1;
		else if (x > mapWidth - 2) x = mapWidth - 2;

		if (y < 1)
			y = 1;
		else if (y > mapHeight - 2)
			y = mapHeight - 2;

		//this handles the linked list between particles on the same square
		p.next = map[x][y];
		map[x][y] = &p;

		mapCoords[i][0] = x;
		mapCoords[i][1] = y;
	}
}

//saves neighbors for lookup in other functions
void storeNeighbors(){
	for (int i = 0; i<totalParticles; i++){
		Particle& p = particles[i];
		int pX = mapCoords[i][0];
		int pY = mapCoords[i][1];

		neighbours[i].count = 0;

		//iterate over the nine squares on the grid around p
		for (int mapX = pX - 1; mapX <= pX + 1; mapX++){
			for (int mapY = pY - 1; mapY <= pY + 1; mapY++){
				//go through the current square's linked list of overlapping values, if there is one
				for (Particle* nextDoor = map[mapX][mapY]; nextDoor != NULL; nextDoor = nextDoor->next){
					const Particle& pj = *nextDoor;

					float diffX = pj.posX - p.posX;
					float diffY = pj.posY - p.posY;
					float r2 = diffX*diffX + diffY*diffY;
					float r = sqrt(r2);
					float q = r / particleHeight;

					//save this neighbor
					if (q < 1 && q > 0.0000000000001f){
						if (neighbours[i].count < maxSprings){
							neighbours[i].particles[neighbours[i].count] = &pj;
							neighbours[i].r[neighbours[i].count] = r;
							neighbours[i].Lij[neighbours[i].count] = particleHeight;
							neighbours[i].count++;
						}
					}
				}
			}
		}
	}
}

//2-Dimensional gravity for player input
void applyGravity(){
	for (int i = 0; i<totalParticles; i++){
		Particle& p = particles[i];
		p.velY += gravity.y*deltaTime;
		p.velX += gravity.x*deltaTime;
	}
}

//Advances particle along its velocity
void advance(){
	for (int i = 0; i<totalParticles; i++){
		Particle& p = particles[i];

		prevPos[i].x = p.posX;
		prevPos[i].y = p.posY;

		p.posX += deltaTime * p.velX;
		p.posY += deltaTime * p.velY;
	}
}

//applies viscosity impulses to particles
void applyViscosity(){
	//cycle through all particles
	for (int i = 0; i<totalParticles; i++){
		Particle& p = particles[i];

		for (int j = 0; j < neighbours[i].count; j++){
			const Particle& pNear = *neighbours[i].particles[j];

			float diffX = pNear.posX - p.posX;
			float diffY = pNear.posY - p.posY;

			float r2 = diffX*diffX + diffY*diffY;
			float r = sqrt(r2);

			float q = r / particleHeight;

			if (q>1) continue;

			float diffVelX = p.velX - pNear.velX;
			float diffVelY = p.velY - pNear.velY;
			float u = diffVelX*diffX + diffVelY*diffY;

			if (u > 0){
				float a = 1 - q;
				diffX /= r;
				diffY /= r;
				u /= r;

				float I = 0.5f * deltaTime * a * (linearViscocity*u + quadraticViscocity*u*u);

				particles[i].velX -= I * diffX;
				particles[i].velY -= I * diffY;
			}
		}
	}
}

//This combines both algorithms for spring adjustment and displacement.
void adjustSprings(){
	//iterate through all particles
	for (int i = 0; i < totalParticles; i++)	{
		const Particle& p = particles[i];
		//iterate through that particles neighbors
		for (int j = 0; j < neighbours[i].count; j++){
			const Particle& pNear = *neighbours[i].particles[j];

			float r = neighbours[i].r[j];
			float q = r / particleHeight;

			if (q < 1 && q > 0.0000000000001f){
				float d = yield*neighbours[i].Lij[j];

				//calculate spring values
				if (r>particleHeight + d){
					neighbours[i].Lij[j] += deltaTime*alpha*(r - particleHeight - d);
				}
				else if (r<particleHeight - d){
					neighbours[i].Lij[j] -= deltaTime*alpha*(particleHeight - d - r);
				}

				//apply those changes to the particle
				float Lij = neighbours[i].Lij[j];
				float diffX = pNear.posX - p.posX;
				float diffY = pNear.posY - p.posY;
				float displaceX = deltaTime*deltaTime*kSpring*(1 - Lij / particleHeight)*(Lij - r)*diffX;
				float displaceY = deltaTime*deltaTime*kSpring*(1 - Lij / particleHeight)*(Lij - r)*diffY;
				particles[i].posX -= 0.5f*displaceX;
				particles[i].posY -= 0.5f*displaceY;
			}
		}
	}
}

//This maps pretty closely to the outline in the paper. Find density and pressure for all particles, 
//then apply a displacement based on that. There is an added if statement to handle surface tension for multiple weights of particles
void doubleDensityRelaxation(){
	for (int i = 0; i < totalParticles; i++)	{
		Particle& p = particles[i];

		float density = 0;
		float nearDensity = 0;
		
		for (int j = 0; j < neighbours[i].count; j++){
			const Particle& pNear = *neighbours[i].particles[j];

			float r = neighbours[i].r[j];
			float q = r / particleHeight;
			if (q>1) continue;
			float a = 1 - q;

			density += a*a * pNear.m * 20;
			nearDensity += a*a*a * pNear.m * 30;
		}
		p.pressure = stiffness * (density - p.m*restDensity);
		p.nearPressure = nearStiffness * nearDensity;
		float dx = 0, dy = 0;

		for (int j = 0; j < neighbours[i].count; j++){
			const Particle& pNear = *neighbours[i].particles[j];

			float diffX = pNear.posX - p.posX;
			float diffY = pNear.posY - p.posY;

			float r = neighbours[i].r[j];
			float q = r / particleHeight;
			if (q>1) continue;
			float a = 1 - q;
			float d = (deltaTime*deltaTime) * ((p.nearPressure + pNear.nearPressure)*a*a*a*53 + (p.pressure + pNear.pressure)*a*a*35) / 2;

			// weight is added to the denominator to reduce the change in dx based on its weight
			dx -= d * diffX / (r*p.m);
			dy -= d * diffY / (r*p.m);

			//surface tension is mapped with one type of particle, 
			//this allows multiple weights of particles to behave appropriately
			if (p.m == pNear.m && multipleParticleTypes == true){
				dx += surfaceTension * diffX;
				dy += surfaceTension * diffY;
			}
		}

		p.posX += dx;
		p.posY += dy;
	}
}

void computeNextVelocity(){
	for (int i = 0; i<totalParticles; ++i){
		Particle& p = particles[i];
		p.velX = (p.posX - prevPos[i].x) / deltaTime;
		p.velY = (p.posY - prevPos[i].y) / deltaTime;
	}
}

//Only checks if particles have left window, and push them back if so
void resolveCollisions(){
	for (int i = 0; i<totalParticles; i++){
		Particle& p = particles[i];

		for (int j = 0; j<4; j++){
			const Boundary& boundary = boundaries[j];
			float distance = boundary.x*p.posX + boundary.y*p.posY - boundary.c;

			if (distance < radius){
				if (distance < 0)
					distance = 0;
				p.velX += (radius - distance) * boundary.x / deltaTime;
				p.velY += (radius - distance) * boundary.y / deltaTime;
			}

			//The resolve collisions tends to overestimate the needed counter velocity, this limits that
			if (p.velX > 10) p.velX = 10;
			if (p.velY > 10) p.velY = 10;
			if (p.velX < -10) p.velX = -10;
			if (p.velY < -10) p.velY = -10;
		}
	}
}

//Iterates through every particle and multiplies its RGB values based on speed. 
//speed^2 is just used to make the difference in speeds more noticeable.
void adjustColor(){
	for (int i = 0; i<totalParticles; i++){
		const Particle& p = particles[i];
		const ParticleType& pt = particleTypes[i];

		float speed2 = p.velX*p.velX + p.velY*p.velY;

		Vec4& color = savedParticleColors[i];
		color = pt.color;
		color.r *= 1.0f + velocityFactor*speed2;
		color.g *= 1.0f + velocityFactor*speed2;
		color.b *= 1.0f + velocityFactor*speed2;
	}
}

//OpenGL operations. Because I'm just coloring in points, I don't need to deal with fragment shaders
void render(){
	glClearColor(0.5f, 0.5f, 0.5f, 1);
	glClear(GL_COLOR_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, viewWidth, 0, viewHeight, 0, 1);

	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	adjustColor();

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	glPointSize(drawRatio*radius*screenWidth / viewWidth);

	glColorPointer(4, GL_FLOAT, sizeof(Vec4), savedParticleColors);
	glVertexPointer(2, GL_FLOAT, sizeof(Particle), particles);
	glDrawArrays(GL_POINTS, 0, totalParticles);

	glDisableClientState(GL_COLOR_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);
}

//Source particle generation
int delay = 0;
void generateParticles(){
	if (totalParticles == maxParticles)
		return;

	if (delay++ < 2)
		return;

	for (int turn = 0; turn<2; turn++){
		Source& source = sources[turn];
		if (source.count >= maxParticles / 4) continue;

		for (int i = 0; i <= 2 && totalParticles<maxParticles; i++){
			Particle& p = particles[totalParticles];
			ParticleType& pt = particleTypes[totalParticles];
			totalParticles++;

			source.count++;

			//for an even distribution of particles
			float offset = float(i) / 1.5f;
			offset *= 0.2f;
			p.posX = source.position.x - offset*source.direction.y;
			p.posY = source.position.y + offset*source.direction.x;
			p.velX = source.speed *source.direction.x;
			p.velY = source.speed *source.direction.y;
			p.m = source.pt.mass;

			pt = source.pt;
		}
	}
	delay = 0;
}

//Mouse particle generator. generates less particles, but does it more consistently. 
//It can also only be run every frame, while generate can run 7 times a frame
void addParticles(int x, int y){
	if (totalParticles == maxParticles)
		return;

	float localX = (viewWidth * (float(x) / screenWidth)),
		localY = (viewHeight * (float(screenHeight - y) / screenHeight));
	sources[2].position = Vec2(localX, localY);

	for (int i = 0; i <= 2 && totalParticles<maxParticles; i++){
		Particle& p = particles[totalParticles];
		ParticleType& pt = particleTypes[totalParticles];
		totalParticles++;

		sources[2].count++;

		p.posX = sources[2].position.x + (particleHeight / 2)*i;
		p.posY = sources[2].position.y;
		p.velX = sources[2].speed * sources[2].direction.x;
		p.velY = sources[2].speed * sources[2].direction.y;
		p.m = sources[2].pt.mass;

		pt = sources[2].pt;
	}
}

//Runs through all of the logic 7 times a frame
void update(){
	for (int step = 0; step<timeStep; step++){
		generateParticles();
		applyGravity();	
		applyViscosity();
		advance();
		adjustSprings();
		updateMap();
		storeNeighbors();
		doubleDensityRelaxation();
		computeNextVelocity();
		resolveCollisions();
	}
}

int main(int argc, char** argv){
	SDL_Init(SDL_INIT_VIDEO);
	screen = SDL_CreateWindow("SPH", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, screenWidth, screenHeight, SDL_WINDOW_OPENGL);
	gContext = SDL_GL_CreateContext(screen);

	SDL_GL_SetSwapInterval(1);

	//Main loop flag
	bool quit = false;
	SDL_Event e;
	SDL_StartTextInput();
	memset(particles, 0, maxParticles*sizeof(Particle));
	updateMap();

	while (!quit){
		while (SDL_PollEvent(&e) != 0){
			if (e.type == SDL_QUIT)
				quit = true;
			else{
				handleEvents(e);
			}
		}
		if (activeSpout){
			int x = 0;
			int y = 0;
			SDL_GetMouseState(&x, &y);
			addParticles(x, y);
		}
		update(); // update logic
		render(); // update buffer

		SDL_GL_SwapWindow(screen); //update window
	}

	//close program, return true
	SDL_StopTextInput();
	SDL_DestroyWindow(screen);
	screen = NULL;
	SDL_Quit();
	return 0;
}

//this is where i handle all of the input from the user
void handleEvents(SDL_Event e){
	if (e.type == SDL_KEYDOWN){
		if (e.key.keysym.sym == SDLK_DOWN){
			gravity.x = 0.f;
			gravity.y = -9.8f;
		}
		else if (e.key.keysym.sym == SDLK_UP){
			gravity.x = 0.f;
			gravity.y = 9.8f;
		}
		else if (e.key.keysym.sym == SDLK_LEFT){
			gravity.x = -9.8f;
			gravity.y = 0.f;
		}
		else if (e.key.keysym.sym == SDLK_RIGHT){
			gravity.x = 9.8f;
			gravity.y = 0.f;
		}
		else if (e.key.keysym.sym == SDLK_SPACE){
			gravity.x = 0.f;
			gravity.y = 0.f;
		}
		else if (e.key.keysym.sym == SDLK_KP_PLUS){
			if (sources[2].pt == green){
				printf("1");
				sources[2].pt = blue;
			}
			else if (sources[2].pt == red){
				printf("2");
				sources[2].pt = green;
			}
		}
		else if (e.key.keysym.sym == SDLK_KP_MINUS){
			if (sources[2].pt == green){
				printf("3");
				sources[2].pt = red;
			}
			else if (sources[2].pt == blue){
				printf("4");
				sources[2].pt = green;
			}
		}
	}
	if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT){
		activeSpout = true;
	}
	else if (e.type == SDL_MOUSEBUTTONUP && e.button.button == SDL_BUTTON_LEFT){
		activeSpout = false;
	}
}