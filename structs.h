struct Particle{
	float posX;
	float posY;

	float velX;
	float velY;

	float m; //mass

	float density;
	float nearDensity;
	float pressure;
	float nearPressure;

	Particle* next;
};

//Spring properties
#define maxSprings 64
#define alpha 0.3f
#define kSpring 0.3f
struct Springs{
	const Particle* particles[maxSprings];
	float r[maxSprings];
	float Lij[maxSprings];
	int count;
};

struct Vec2{
	Vec2() { }
	Vec2(float x, float y) : x(x), y(y) { }
	float x;
	float y;
	float magnitude(){
		return sqrt(x*x+y*y);
	}
	Vec2 operator + (Vec2 v){
		return Vec2(x+v.x, y+v.y);
	}
	Vec2 operator - (Vec2 v){
		return Vec2(x - v.x, y - v.y);
	}
};

struct Boundary{
	Boundary() { }
	Boundary(float x, float y, float c) : x(x), y(y), c(c) { }
	float x;
	float y;
	float c;
};

struct Vec4{
	Vec4() { }
	Vec4(float r, float g, float b, float a) : r(r), g(g), b(b), a(a) { }
	float r, g, b, a;
	bool operator == (Vec4 v){
		return r==v.r && g==v.g && b==v.b && a==v.a;
	}
};

struct ParticleType{
	ParticleType() { }
	ParticleType(const Vec4& color, float mass) : color(color), mass(mass) { }
	Vec4 color;
	float mass;
	bool operator == (ParticleType pt){
		return mass == pt.mass && color == pt.color;
	}
};


struct Source{
	Source(const ParticleType& pt, const Vec2& position, const Vec2& direction, float speed)
		: pt(pt), position(position), direction(direction), speed(speed), count(0)
	{
		float len = sqrt(direction.x*direction.x + direction.y*direction.y);
		this->direction.x /= len;
		this->direction.y /= len;
	}
	ParticleType pt;
	Vec2 position;
	Vec2 direction;
	float speed;
	int count;
};