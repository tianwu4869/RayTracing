#pragma once
#include "image.h"


class coordinates {
public:
	float x, y, z;
	coordinates() {}
	coordinates(float m, float n, float p) {
		x = m; y = n; z = p;
	}
};

class color {
	public:
		float r, g, b;
		color() {}
		color(float o, float p, float q) {
			r = o; g = p; b = q;
		}
};

class vector {
public:
	float norm2;
	float x, y, z;
	vector() {}
	vector(float m, float n, float p) {
		x = m; y = n; z = p;
		norm2 = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
	}
	vector(coordinates A, coordinates B) {
		x = B.x - A.x;
		y = B.y - A.y;
		z = B.z - A.z;
		norm2 = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
	}
	void normalization() {
		x = x / norm2;
		y = y / norm2;
		z = z / norm2;
	}
	void computeModule() {
		norm2 = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
	}
};

float dotProduct(vector A, vector B);

vector crossProduct(vector A, vector B);

class shape {

};

class spheres : public shape {
	public:
		spheres* next;
		float r;
		float ar, ag, ab, dr, dg, db, sr, sg, sb, ns, tr, tg, tb, ior;
		coordinates center;
};

class triangle : public shape {
public:
	int vertex1, vertex2, vertex3;
	float ar, ag, ab, dr, dg, db, sr, sg, sb, ns, tr, tg, tb, ior;
	triangle* next;
};

class light {

};

class pointlights : public light {
public:
	pointlights *next;
	coordinates center;
	color colors;
};

class directionalLight : public light {
public:
	vector direction;
	color colors;
};

class spotlight : public light {
public:
	coordinates center;
	color colors;
	vector centralline;
	float angle1, angle2;
	spotlight* next;
};

class camera {
	public:
		coordinates p;
		vector d, u;
		float ha;
		camera(float px, float py, float pz, float dx, float dy, float dz, float ux, float uy, float uz, float h) {
			p.x = px; p.y = py; p.z = pz;
			d.x = dx; d.y = dy; d.z = dz;
			u.x = ux; u.y = uy; u.z = uz;
			ha = h;
			d.computeModule();
			u.computeModule();
		}
};

class origin {
	virtual void f() {

	}
};

class sphereAndt : public origin {
public:
	spheres* nearestPoint;
	float t;
};

class triangleAndt : public origin {
public:
	triangle* nearestPoint;
	vector normalVector;
	float t;
};

class ray {
	public:
		coordinates startPoint;
		vector direction;
		ray(coordinates c, vector v) {
			startPoint.x = c.x; startPoint.y = c.y; startPoint.z = c.z;
			direction.x = v.x; direction.y = v.y; direction.z = v.z;
			direction.computeModule();
		}
};

float rgbThreshold(float color);

void setBackgroundColor(float r, float g, float b, Image* img);

void setRelativeCoordinate();

void setImageNull(Image* img);

//node_spheres* sphereListing(node_spheres* sphere_header, int sphere_count);

int shadow_pointlight(coordinates intersect, pointlights *pl, coordinates* vertex);

int shadow_directionallight(coordinates intersect, directionalLight *pl, coordinates* vertex);

int shadow_spotlight(coordinates intersect, spotlight *pl, coordinates* vertex);

Image* rayTracing(Image* img, coordinates* vertex);

triangleAndt triangleRayIntersection(coordinates* vertex, ray R);

sphereAndt sphereRayIntersection(ray R);

color recursiveRay(int depth, ray R, coordinates* vertex);
