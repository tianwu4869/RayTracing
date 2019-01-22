// hw2Raytracing.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "image.h"
#include <math.h>
#include <iostream>
#include <algorithm>
#include "hw2Raytracing.h"
#include <omp.h>
#include <ctime>
using namespace std;

extern camera cam;
extern directionalLight* dl_header;
extern color ambientLight;
extern color background;
extern triangle *triangle_header;
extern int triangle_count;

extern spheres *sphere_header;
extern int spheres_count;

extern pointlights *pointlight_header;
extern int pointlight_count;
vector U, V;
extern int max_depth;
extern spotlight* spotlight_header;
extern int spotlight_count;
extern spotlight* spotlight_temp;

float dotProduct(vector A, vector B) {
	float result;
	result = A.x * B.x + A.y * B.y + A.z * B.z;
	return result;
}

vector crossProduct(vector A, vector B) {
	float x, y, z;
	x = A.y * B.z - A.z * B.y;
	y = A.z * B.x - A.x * B.z;
	z = A.x * B.y - A.y * B.x;
	vector C(x, y, z);
	return C;
}

float rgbThreshold(float color) {
	if (color > 255) {
		color = 255;
	}
	else if (color < 0) {
		color = 0;
	}
	return color;
}

void setBackgroundColor(float r, float g, float b, Image* img) {
	int h = img->Height();
	int w = img->Width();
	for (int i = 0; i < w; i++) {
		for (int j = 0; j < h; j++) {
			Pixel p;
			p.r = 255 * r;
			p.g = 255 * g;
			p.b = 255 * b;
			img->SetPixel(i, j, p);	
		}
	}
}

void setRelativeCoordinate() {
	U = crossProduct(cam.d, cam.u);
	V = crossProduct(U, cam.d);
	U.normalization();
	V.normalization();
	cam.d.normalization();
}

void setImageNull(Image* img) {
	int h = img->Height();
	int w = img->Width();
	for (int i = 0; i < w; i++) {
		for (int j = 0; j < h; j++) {
			Pixel p;
			p.r = NULL;
			p.g = NULL;
			p.b = NULL;
			img->SetPixel(i, j, p);
		}
	}
}

//node_spheres* sphereListing(node_spheres* sphere_header, int sphere_count) {
//	for (int i = 0; i < sphere_count - 1; i++) {
//		node_spheres* current, *former;
//		current = sphere_header;
//		former = sphere_header;
//		
//		for (int j = 0; j < sphere_count - 1; j++) {
//			float distance1, distance2;
//			distance1 = sqrt(pow(current->x - px, 2) + pow(current->y - py, 2) + pow(current->z - pz, 2));
//			distance2 = sqrt(pow(current->next->x - px, 2) + pow(current->next->y - py, 2) + pow(current->next->z - pz, 2));
//			if (distance1 < distance2) {
//				if (j == 0) {
//					sphere_header = sphere_header->next;
//					current->next = sphere_header->next;
//					sphere_header->next = current;
//					former = sphere_header;
//				}
//				else {
//					former->next = current->next;
//					current->next = current->next->next;
//					former = former->next;
//					former->next = current;
//				}
//			}
//			else {
//				if (j == 0) {
//					current = current->next;
//				}
//				else{
//					current = current->next;
//					former = former->next;
//				}
//			}
//		}
//	}
//
//	return sphere_header;
//}

int shadow_pointlight(coordinates intersect, pointlights *pl, coordinates* vertex) {
	int shadow;
	float jittered1, jittered2, mod = 10;
	pointlights current, *c = NULL;
	c = &current;
	c->center = pl->center;
	c->colors = pl->colors;
	c->next = pl->next;
	jittered1 = rand() % 10 / mod - 0.5;
	jittered2 = rand() % 10 / mod - 0.5;
	vector ll(intersect, c->center);
	ll.normalization();
	vector a(-1 / ll.x, 2 / ll.y, -1 / ll.z);
	vector b(2 / ll.x, -1 / ll.y, -1 / ll.z);
	a.normalization();
	b.normalization();

	c->center.x = c->center.x + jittered1 * a.x + jittered2 * b.x;
	c->center.y = c->center.y + jittered1 * a.y + jittered2 * b.y;
	c->center.z = c->center.z + jittered1 * a.z + jittered2 * b.z;

	spheres* s = sphere_header;
	shadow = 1;
	for (int m = 0; m < spheres_count; m++) {
		//if (s->center.x != sph->center.x || s->center.y != sph->center.y || s->center.z != sph->center.z || s->r != sph->r) {
			vector l(intersect, c->center);
			float a, b, c, delta, t, t_module;
			a = pow(l.x, 2) + pow(l.y, 2) + pow(l.z, 2);
			b = 2 * (intersect.x - s->center.x) * l.x + 2 * (intersect.y - s->center.y) * l.y + 2 * (intersect.z - s->center.z) * l.z;
			c = pow((intersect.x - s->center.x), 2) + pow((intersect.y - s->center.y), 2) + pow((intersect.z - s->center.z), 2) - pow(s->r, 2);
			delta = pow(b, 2) - 4 * a * c;

			if (delta > 0) {
				t = (-b - sqrt(delta)) / (2 * a);
				t_module = sqrt(pow(t * l.x, 2) + pow(t * l.y, 2) + pow(t * l.z, 2));
				if (t > 0 && t_module < l.norm2) {
					shadow = 0;
					break;
				}
			}
		//}
		s = s->next;
	}

	triangle* t = triangle_header;
	if (shadow == 1) {
		for (int k = 0; k < triangle_count; k++) {
			//if (t->vertex1 != tri->vertex1 || t->vertex2 != tri->vertex2 || t->vertex3 != tri->vertex3) {
				coordinates v1(vertex[t->vertex1].x, vertex[t->vertex1].y, vertex[t->vertex1].z);
				coordinates v2(vertex[t->vertex2].x, vertex[t->vertex2].y, vertex[t->vertex2].z);
				coordinates v3(vertex[t->vertex3].x, vertex[t->vertex3].y, vertex[t->vertex3].z);

				vector v1v2(v1, v2);
				vector v2v3(v2, v3);
				vector v3v1(v3, v1);

				vector normalVector;// normal vector of the triangle plain 
				normalVector = crossProduct(v1v2, v2v3);

				vector shadowed(intersect, c->center);
				ray R(intersect, shadowed);
				float m;
				coordinates p;// point where ray intersects with the triangle plain
				m = (normalVector.x * (v1.x - R.startPoint.x) + normalVector.y * (v1.y - R.startPoint.y) + normalVector.z * (v1.z - R.startPoint.z)) / (normalVector.x * R.direction.x + normalVector.y * R.direction.y + normalVector.z * R.direction.z);
				if (m > 0) {
					p.x = R.startPoint.x + R.direction.x * m;
					p.y = R.startPoint.y + R.direction.y * m;
					p.z = R.startPoint.z + R.direction.z * m;

					int inside = 1;
					vector v1p(v1, p);
					vector v1v3(v1, v3);
					if (dotProduct(crossProduct(v1p, v1v2), crossProduct(v1v3, v1v2)) < 0) {// see if p is inside triangle
						inside = 0;
					}
					vector v2p(v2, p);
					vector v2v1(v2, v1);
					if (dotProduct(crossProduct(v2p, v2v3), crossProduct(v2v1, v2v3)) < 0) {
						inside = 0;
					}
					vector v3p(v3, p);
					vector v3v2(v3, v2);
					if (dotProduct(crossProduct(v3p, v3v1), crossProduct(v3v2, v3v1)) < 0) {
						inside = 0;
					}

					if (inside == 1) {
						shadow = 0;
						break;
					}
				}
			//}
			t = t->next;
		}
	}
	return shadow;
}

int shadow_directionallight(coordinates intersect, directionalLight *pl, coordinates* vertex) {
	int shadow;
	float jittered1, mod = 10;
	directionalLight current, *c = NULL;
	c = &current;
	c->direction = pl->direction;
	c->colors = pl->colors;
	//jittered1 = rand() % 10 / mod - 0.5; soft shadow is not applied to directional light
	jittered1 = 0;

	c->direction.x = c->direction.x + jittered1;
	c->direction.y = c->direction.y + jittered1;
	c->direction.z = c->direction.z + jittered1;
	
	spheres* s = sphere_header;
	shadow = 1;
	for (int m = 0; m < spheres_count; m++) {
		//if (s->center.x != sph->center.x || s->center.y != sph->center.y || s->center.z != sph->center.z || s->r != sph->r) {
			vector l;
			l.x = c->direction.x * (-1);
			l.y = c->direction.y * (-1);
			l.z = c->direction.z * (-1);
			l.computeModule();
			l.normalization();
			float a, b, c, delta, t, t_module;
			a = pow(l.x, 2) + pow(l.y, 2) + pow(l.z, 2);
			b = 2 * (intersect.x - s->center.x) * l.x + 2 * (intersect.y - s->center.y) * l.y + 2 * (intersect.z - s->center.z) * l.z;
			c = pow((intersect.x - s->center.x), 2) + pow((intersect.y - s->center.y), 2) + pow((intersect.z - s->center.z), 2) - pow(s->r, 2);
			delta = pow(b, 2) - 4 * a * c;

			if (delta > 0) {
				t = (-b - sqrt(delta)) / (2 * a);
				t_module = sqrt(pow(t * l.x, 2) + pow(t * l.y, 2) + pow(t * l.z, 2));
				if (t > 0 && t_module < l.norm2) {
					shadow = 0;
					break;
				}
			}
		//}
		s = s->next;
	}

	triangle* t = triangle_header;
	vector l;
	l.x = c->direction.x * (-1);
	l.y = c->direction.y * (-1);
	l.z = c->direction.z * (-1);
	l.computeModule();
	l.normalization();
	ray R(intersect, l);

	if (shadow == 1) {
		for (int k = 0; k < triangle_count; k++) {
			//if (tri->vertex1 != t->vertex1 || tri->vertex2 != t->vertex2 || tri->vertex3 != t->vertex3) {
				coordinates v1(vertex[t->vertex1].x, vertex[t->vertex1].y, vertex[t->vertex1].z);
				coordinates v2(vertex[t->vertex2].x, vertex[t->vertex2].y, vertex[t->vertex2].z);
				coordinates v3(vertex[t->vertex3].x, vertex[t->vertex3].y, vertex[t->vertex3].z);

				vector v1v2(v1, v2);
				vector v2v3(v2, v3);
				vector v3v1(v3, v1);

				vector normalVector;// normal vector of the triangle plain 
				normalVector = crossProduct(v1v2, v2v3);
				normalVector.normalization();

				float m;
				coordinates p;// point where ray intersects with the triangle plain
				m = (normalVector.x * (v1.x - R.startPoint.x) + normalVector.y * (v1.y - R.startPoint.y) + normalVector.z * (v1.z - R.startPoint.z)) / (normalVector.x * R.direction.x + normalVector.y * R.direction.y + normalVector.z * R.direction.z);
				p.x = R.startPoint.x + R.direction.x * m;
				p.y = R.startPoint.y + R.direction.y * m;
				p.z = R.startPoint.z + R.direction.z * m;
				vector intersect2p(intersect, p);
				if ((m - 0.005) > 0) {

					int inside = 1;
					vector v1p(v1, p);
					vector v1v3(v1, v3);
					if (dotProduct(crossProduct(v1p, v1v2), crossProduct(v1v3, v1v2)) < 0) {// see if p is inside triangle
						inside = 0;
					}
					vector v2p(v2, p);
					vector v2v1(v2, v1);
					if (dotProduct(crossProduct(v2p, v2v3), crossProduct(v2v1, v2v3)) < 0) {
						inside = 0;
					}
					vector v3p(v3, p);
					vector v3v2(v3, v2);
					if (dotProduct(crossProduct(v3p, v3v1), crossProduct(v3v2, v3v1)) < 0) {
						inside = 0;
					}

					if (inside == 1) {
						shadow = 0;
						break;
					}
				}
			//}
			t = t->next;
		}
	}
	return shadow;
}

int shadow_spotlight(coordinates intersect, spotlight *pl, coordinates* vertex) {
	int shadow;
	float jittered1, jittered2, mod = 10;
	spotlight current, *c = NULL;
	c = &current;
	c->center = pl->center;
	c->colors = pl->colors;
	jittered1 = rand() % 10 / mod - 0.5;
	jittered2 = rand() % 10 / mod - 0.5;
	vector ll(intersect, c->center);
	ll.normalization();
	vector a(-1 / ll.x, 2 / ll.y, -1 / ll.z);
	vector b(2 / ll.x, -1 / ll.y, -1 / ll.z);
	a.normalization();
	b.normalization();

	c->center.x = c->center.x + jittered1 * a.x + jittered2 * b.x;
	c->center.y = c->center.y + jittered1 * a.y + jittered2 * b.y;
	c->center.z = c->center.z + jittered1 * a.z + jittered2 * b.z;

	spheres* s = sphere_header;
	shadow = 1;
	for (int m = 0; m < spheres_count; m++) {
		//if (s->center.x != sph->center.x || s->center.y != sph->center.y || s->center.z != sph->center.z || s->r != sph->r) {
		vector l(intersect, c->center);
		float a, b, c, delta, t, t_module;
		a = pow(l.x, 2) + pow(l.y, 2) + pow(l.z, 2);
		b = 2 * (intersect.x - s->center.x) * l.x + 2 * (intersect.y - s->center.y) * l.y + 2 * (intersect.z - s->center.z) * l.z;
		c = pow((intersect.x - s->center.x), 2) + pow((intersect.y - s->center.y), 2) + pow((intersect.z - s->center.z), 2) - pow(s->r, 2);
		delta = pow(b, 2) - 4 * a * c;

		if (delta > 0) {
			t = (-b - sqrt(delta)) / (2 * a);
			t_module = sqrt(pow(t * l.x, 2) + pow(t * l.y, 2) + pow(t * l.z, 2));
			if (t > 0 && t_module < l.norm2) {
				shadow = 0;
				break;
			}
		}
		//}
		s = s->next;
	}

	triangle* t = triangle_header;
	vector shadowed(intersect, c->center);
	ray R(intersect, shadowed);

	if (shadow == 1) {
		for (int k = 0; k < triangle_count; k++) {
			//if (tri->vertex1 != t->vertex1 || tri->vertex2 != t->vertex2 || tri->vertex3 != t->vertex3) {
			coordinates v1(vertex[t->vertex1].x, vertex[t->vertex1].y, vertex[t->vertex1].z);
			coordinates v2(vertex[t->vertex2].x, vertex[t->vertex2].y, vertex[t->vertex2].z);
			coordinates v3(vertex[t->vertex3].x, vertex[t->vertex3].y, vertex[t->vertex3].z);

			vector v1v2(v1, v2);
			vector v2v3(v2, v3);
			vector v3v1(v3, v1);

			vector normalVector;// normal vector of the triangle plain 
			normalVector = crossProduct(v1v2, v2v3);
			normalVector.normalization();

			float m;
			coordinates p;// point where ray intersects with the triangle plain
			m = (normalVector.x * (v1.x - R.startPoint.x) + normalVector.y * (v1.y - R.startPoint.y) + normalVector.z * (v1.z - R.startPoint.z)) / (normalVector.x * R.direction.x + normalVector.y * R.direction.y + normalVector.z * R.direction.z);
			p.x = R.startPoint.x + R.direction.x * m;
			p.y = R.startPoint.y + R.direction.y * m;
			p.z = R.startPoint.z + R.direction.z * m;
			vector intersect2p(intersect, p);
			if (m > 0) {

				int inside = 1;
				vector v1p(v1, p);
				vector v1v3(v1, v3);
				if (dotProduct(crossProduct(v1p, v1v2), crossProduct(v1v3, v1v2)) < 0) {// see if p is inside triangle
					inside = 0;
				}
				vector v2p(v2, p);
				vector v2v1(v2, v1);
				if (dotProduct(crossProduct(v2p, v2v3), crossProduct(v2v1, v2v3)) < 0) {
					inside = 0;
				}
				vector v3p(v3, p);
				vector v3v2(v3, v2);
				if (dotProduct(crossProduct(v3p, v3v1), crossProduct(v3v2, v3v1)) < 0) {
					inside = 0;
				}

				if (inside == 1) {
					shadow = 0;
					break;
				}
			}
			//}
			t = t->next;
		}
	}
	return shadow;
}

color recursiveRay(int depth, ray R, coordinates* vertex) {
	color I(0, 0, 0);
	I.r = background.r * 255;
	I.g = background.g * 255;
	I.b = background.b * 255;

	sphereAndt sph;
	sph = sphereRayIntersection(R);
	triangleAndt tri;
	tri = triangleRayIntersection(vertex, R);

	if ((sph.t < tri.t || tri.t < 0) && sph.t > 0) {// intersect with sphere
		I.r = ambientLight.r * 255 * sph.nearestPoint->ar;
		I.g = ambientLight.g * 255 * sph.nearestPoint->ag;
		I.b = ambientLight.b * 255 * sph.nearestPoint->ab;
		coordinates intersect;
		intersect.x = R.startPoint.x + R.direction.x * sph.t;
		intersect.y = R.startPoint.y + R.direction.y * sph.t;
		intersect.z = R.startPoint.z + R.direction.z * sph.t;

		vector n(sph.nearestPoint->center, intersect);
		n.normalization();
		vector v(intersect, R.startPoint);
		v.normalization();

		pointlights *pl = pointlight_header;
		for (int pointlight = 0; pointlight < pointlight_count; pointlight++) {// Influence from point light
			srand(time(NULL));
			color temp(0, 0, 0);
			//#pragma omp parallel for
			for (int k = 0; k < 32; k++) {// This loop is used for sampling if soft shadow is wanted.
				int shadow = shadow_pointlight(intersect, pl, vertex);
				vector l(intersect, pl->center);
				l.normalization();
				
				vector h;
				h.x = l.x + v.x;
				h.y = l.y + v.y;
				h.z = l.z + v.z;
				h.computeModule();
				h.normalization();
				float nl_dotproduct, zero, nh_dotproduct;
				nl_dotproduct = dotProduct(n, l);
				nh_dotproduct = dotProduct(n, h);
				zero = 0;
				temp.r = temp.r + sph.nearestPoint->dr * pl->colors.r * 255 * max(zero, nl_dotproduct) / pow(l.norm2, 2) * shadow + sph.nearestPoint->sr * pl->colors.r * 255 * pow(max(zero, nh_dotproduct), sph.nearestPoint->ns) / pow(l.norm2, 2) * shadow;
				temp.g = temp.g + sph.nearestPoint->dg * pl->colors.g * 255 * max(zero, nl_dotproduct) / pow(l.norm2, 2) * shadow + sph.nearestPoint->sg * pl->colors.g * 255 * pow(max(zero, nh_dotproduct), sph.nearestPoint->ns) / pow(l.norm2, 2) * shadow;
				temp.b = temp.b + sph.nearestPoint->db * pl->colors.b * 255 * max(zero, nl_dotproduct) / pow(l.norm2, 2) * shadow + sph.nearestPoint->sb * pl->colors.b * 255 * pow(max(zero, nh_dotproduct), sph.nearestPoint->ns) / pow(l.norm2, 2) * shadow;
			}
			I.r = I.r + temp.r / 32;
			I.g = I.g + temp.g / 32;
			I.b = I.b + temp.b / 32;
			pl = pl->next;
		}

		spotlight* sl = spotlight_header;
		for (int spotlight = 0; spotlight < spotlight_count; spotlight++) {// influence from spot light
			srand(time(NULL));
			color temp(0, 0, 0);
			//#pragma omp parallel for
			for (int k = 0; k < 32; k++) {// This loop is used for sampling if soft shadow is wanted, the number of k depends on the number of samples.
				vector spot2intersect(sl->center, intersect);
				spot2intersect.normalization();
				sl->centralline.computeModule();
				sl->centralline.normalization();
				float angle = acos(dotProduct(spot2intersect, sl->centralline)) * 180 / 3.1416;
				float alpha;// coefficient about spotlight decline with different angles
				if (angle > sl->angle2) {
					continue;
				}
				else if (angle <= sl->angle2 && angle > sl->angle1) {
					alpha = (angle - sl->angle1) / (sl->angle2 - sl->angle1);
				}
				else if (0 <= angle && angle <= sl->angle1) {
					alpha = 1;
				}
				//cout << alpha;
				int shadow = shadow_spotlight(intersect, sl, vertex);

				vector l(intersect, sl->center);
				l.normalization();
				vector h;
				h.x = l.x + v.x;
				h.y = l.y + v.y;
				h.z = l.z + v.z;
				h.computeModule();
				h.normalization();
				float nl_dotproduct, zero, nh_dotproduct;
				nl_dotproduct = dotProduct(n, l);
				nh_dotproduct = dotProduct(n, h);
				zero = 0;
				temp.r = temp.r + sph.nearestPoint->dr * sl->colors.r * 255 * max(zero, nl_dotproduct) / pow(l.norm2, 2) * shadow * alpha + sph.nearestPoint->sr * sl->colors.r * 255 * pow(max(zero, nh_dotproduct), sph.nearestPoint->ns) / pow(l.norm2, 2) * shadow * alpha;
				temp.g = temp.g + sph.nearestPoint->dg * sl->colors.g * 255 * max(zero, nl_dotproduct) / pow(l.norm2, 2) * shadow * alpha + sph.nearestPoint->sg * sl->colors.g * 255 * pow(max(zero, nh_dotproduct), sph.nearestPoint->ns) / pow(l.norm2, 2) * shadow * alpha;
				temp.b = temp.b + sph.nearestPoint->db * sl->colors.b * 255 * max(zero, nl_dotproduct) / pow(l.norm2, 2) * shadow * alpha + sph.nearestPoint->sb * sl->colors.b * 255 * pow(max(zero, nh_dotproduct), sph.nearestPoint->ns) / pow(l.norm2, 2) * shadow * alpha;
			}
			I.r = I.r + temp.r / 32;
			I.g = I.g + temp.g / 32;
			I.b = I.b + temp.b / 32;
			sl = sl->next;
		}

		if (dl_header != NULL) {// influene from directional light
			srand(time(NULL));
			color temp(0, 0, 0);
			for (int k = 0; k < 32; k++) {
				int shadow = shadow_directionallight(intersect, dl_header, vertex);
				vector l;
				l.x = dl_header->direction.x * (-1);
				l.y = dl_header->direction.y * (-1);
				l.z = dl_header->direction.z * (-1);
				l.computeModule();
				l.normalization();
				vector h;
				h.x = l.x + v.x;
				h.y = l.y + v.y;
				h.z = l.z + v.z;
				h.computeModule();
				h.normalization();
				float nl_dotproduct, zero, nh_dotproduct;
				nl_dotproduct = dotProduct(n, l);
				nh_dotproduct = dotProduct(n, h);
				zero = 0;
				temp.r = temp.r + sph.nearestPoint->dr * dl_header->colors.r * 255 * max(zero, nl_dotproduct) * shadow + sph.nearestPoint->sr * dl_header->colors.r * 255 * pow(max(zero, nh_dotproduct), sph.nearestPoint->ns) * shadow;
				temp.g = temp.g + sph.nearestPoint->dg * dl_header->colors.g * 255 * max(zero, nl_dotproduct) * shadow + sph.nearestPoint->sg * dl_header->colors.g * 255 * pow(max(zero, nh_dotproduct), sph.nearestPoint->ns) * shadow;
				temp.b = temp.b + sph.nearestPoint->db * dl_header->colors.b * 255 * max(zero, nl_dotproduct) * shadow + sph.nearestPoint->sb * dl_header->colors.b * 255 * pow(max(zero, nh_dotproduct), sph.nearestPoint->ns) * shadow;
			}
			I.r = I.r + temp.r / 32;
			I.g = I.g + temp.g / 32;
			I.b = I.b + temp.b / 32;
		}

		if (depth < max_depth) {// recursion
			float nv_dotproduct = dotProduct(n, v);
			vector r;
			r.x = -v.x + 2 * nv_dotproduct * n.x;
			r.y = -v.y + 2 * nv_dotproduct * n.y;
			r.z = -v.z + 2 * nv_dotproduct * n.z;
			r.computeModule();
			depth++;
			R.direction = r;
			R.startPoint = intersect;
			color temp_I = recursiveRay(depth, R, vertex);// reflection
			I.r = I.r + temp_I.r * sph.nearestPoint->sr;
			I.g = I.g + temp_I.g * sph.nearestPoint->sg;
			I.b = I.b + temp_I.b * sph.nearestPoint->sb;

			vector T;
			float cosRefractionAngle = 1 - (1 - pow(nv_dotproduct, 2)) / pow(sph.nearestPoint->ior, 2);
			//cout << cosRefractionAngle;
			if (cosRefractionAngle >= 0) {// refraction
				cosRefractionAngle = sqrt(cosRefractionAngle);
				T.x = (nv_dotproduct / sph.nearestPoint->ior - cosRefractionAngle) * n.x - v.x / sph.nearestPoint->ior;
				T.y = (nv_dotproduct / sph.nearestPoint->ior - cosRefractionAngle) * n.y - v.y / sph.nearestPoint->ior;
				T.z = (nv_dotproduct / sph.nearestPoint->ior - cosRefractionAngle) * n.z - v.z / sph.nearestPoint->ior;

				float a, b, c, delta, t;
				a = pow(T.x, 2) + pow(T.y, 2) + pow(T.z, 2);
				b = 2 * (intersect.x - sph.nearestPoint->center.x) * T.x + 2 * (intersect.y - sph.nearestPoint->center.y) * T.y + 2 * (intersect.z - sph.nearestPoint->center.z) * T.z;
				c = pow((intersect.x - sph.nearestPoint->center.x), 2) + pow((intersect.y - sph.nearestPoint->center.y), 2) + pow((intersect.z - sph.nearestPoint->center.z), 2) - pow(sph.nearestPoint->r, 2);
				delta = pow(b, 2) - 4 * a * c;
				t = (-b + sqrt(delta)) / (2 * a);
				if (t > 0) {

					coordinates exitPoint;
					exitPoint.x = intersect.x + T.x * t;
					exitPoint.y = intersect.y + T.y * t;
					exitPoint.z = intersect.z + T.z * t;

					vector v2(-v.x, -v.y, -v.z);
					v2.normalization();
					R.direction = v2;
					R.startPoint = exitPoint;
					color temp_I = recursiveRay(depth, R, vertex);
					I.r = I.r + temp_I.r * sph.nearestPoint->tr;
					I.g = I.g + temp_I.g * sph.nearestPoint->tg;
					I.b = I.b + temp_I.b * sph.nearestPoint->tb;
				}
			}
		}
	}
	else if ((tri.t < sph.t || sph.t <0) && tri.t >0) {// intersect with triangle
		I.r = ambientLight.r * 255 * tri.nearestPoint->ar;
		I.g = ambientLight.g * 255 * tri.nearestPoint->ag;
		I.b = ambientLight.b * 255 * tri.nearestPoint->ab;

		coordinates intersect;
		intersect.x = R.startPoint.x + R.direction.x * tri.t;
		intersect.y = R.startPoint.y + R.direction.y * tri.t;
		intersect.z = R.startPoint.z + R.direction.z * tri.t;

		vector v(intersect, R.startPoint);
		v.normalization();

		pointlights *pl = pointlight_header;
		for (int pointlight = 0; pointlight < pointlight_count; pointlight++) {// influence from point light
			srand(time(NULL));
			color temp(0, 0, 0);
			for (int k = 0; k < 32; k++) {
				int shadow = shadow_pointlight(intersect, pl, vertex);
				vector l(intersect, pl->center);
				l.normalization();

				vector h;
				h.x = l.x + v.x;
				h.y = l.y + v.y;
				h.z = l.z + v.z;
				h.computeModule();
				h.normalization();
				float nl_dotproduct, zero, nh_dotproduct;
				if (dotProduct(v, tri.normalVector) < 0) {
					tri.normalVector.x = (-1) * tri.normalVector.x;
					tri.normalVector.y = (-1) * tri.normalVector.y;
					tri.normalVector.z = (-1) * tri.normalVector.z;
				}
				nl_dotproduct = dotProduct(tri.normalVector, l);
				nh_dotproduct = dotProduct(tri.normalVector, h);
				zero = 0;
				temp.r = temp.r + tri.nearestPoint->dr * pl->colors.r * 255 * max(zero, nl_dotproduct) / pow(l.norm2, 2) * shadow + tri.nearestPoint->sr * pl->colors.r * 255 * pow(max(zero, nh_dotproduct), tri.nearestPoint->ns) / pow(l.norm2, 2) * shadow;
				temp.g = temp.g + tri.nearestPoint->dg * pl->colors.g * 255 * max(zero, nl_dotproduct) / pow(l.norm2, 2) * shadow + tri.nearestPoint->sg * pl->colors.g * 255 * pow(max(zero, nh_dotproduct), tri.nearestPoint->ns) / pow(l.norm2, 2) * shadow;
				temp.b = temp.b + tri.nearestPoint->db * pl->colors.b * 255 * max(zero, nl_dotproduct) / pow(l.norm2, 2) * shadow + tri.nearestPoint->sb * pl->colors.b * 255 * pow(max(zero, nh_dotproduct), tri.nearestPoint->ns) / pow(l.norm2, 2) * shadow;
			}
			I.r = I.r + temp.r / 32;
			I.g = I.g + temp.g / 32;
			I.b = I.b + temp.b / 32;
			pl = pl->next;
		}

		spotlight* sl = spotlight_header;
		for (int spotlight = 0; spotlight < spotlight_count; spotlight++) {// influence from spot light
			srand(time(NULL));
			color temp(0, 0, 0);
			//#pragma omp parallel for
			for (int k = 0; k < 32; k++) { // This loop is used for sampling if soft shadow is wanted, the number of k depends on the number of samples.
				vector spot2intersect(sl->center, intersect);
				spot2intersect.normalization();
				sl->centralline.computeModule();
				sl->centralline.normalization();
				float angle = acos(dotProduct(spot2intersect, sl->centralline)) * 180 / 3.1416;
				float alpha;// coefficient about spotlight decline with different angles
				if (angle > sl->angle2) {
					continue;
				}
				else if (angle <= sl->angle2 && angle > sl->angle1) {
					alpha = (angle - sl->angle1) / (sl->angle2 - sl->angle1);
				}
				else if (0 <= angle && angle <= sl->angle1) {
					alpha = 1;
				}

				int shadow = shadow_spotlight(intersect, sl, vertex);

				vector l(intersect, sl->center);
				l.normalization();
				vector h;
				h.x = l.x + v.x;
				h.y = l.y + v.y;
				h.z = l.z + v.z;
				h.computeModule();
				h.normalization();
				float nl_dotproduct, zero, nh_dotproduct;
				if (dotProduct(v, tri.normalVector) < 0) {
					tri.normalVector.x = (-1) * tri.normalVector.x;
					tri.normalVector.y = (-1) * tri.normalVector.y;
					tri.normalVector.z = (-1) * tri.normalVector.z;
				}
				nl_dotproduct = dotProduct(tri.normalVector, l);
				nh_dotproduct = dotProduct(tri.normalVector, h);
				zero = 0;
				temp.r = temp.r + tri.nearestPoint->dr * sl->colors.r * 255 * max(zero, nl_dotproduct) / pow(l.norm2, 2) * shadow * alpha + tri.nearestPoint->sr * sl->colors.r * 255 * pow(max(zero, nh_dotproduct), tri.nearestPoint->ns) / pow(l.norm2, 2) * shadow * alpha;
				temp.g = temp.g + tri.nearestPoint->dg * sl->colors.g * 255 * max(zero, nl_dotproduct) / pow(l.norm2, 2) * shadow * alpha + tri.nearestPoint->sg * sl->colors.g * 255 * pow(max(zero, nh_dotproduct), tri.nearestPoint->ns) / pow(l.norm2, 2) * shadow * alpha;
				temp.b = temp.b + tri.nearestPoint->db * sl->colors.b * 255 * max(zero, nl_dotproduct) / pow(l.norm2, 2) * shadow * alpha + tri.nearestPoint->sb * sl->colors.b * 255 * pow(max(zero, nh_dotproduct), tri.nearestPoint->ns) / pow(l.norm2, 2) * shadow * alpha;
			}
			I.r = I.r + temp.r / 32;
			I.g = I.g + temp.g / 32;
			I.b = I.b + temp.b / 32;
			sl = sl->next;
		}

		if (dl_header != NULL) {// influence from directional light
			srand(time(NULL));
			color temp(0, 0, 0);
			for (int k = 0; k < 32; k++) {
				int shadow = shadow_directionallight(intersect, dl_header, vertex);
				vector l;
				l.x = dl_header->direction.x * (-1);
				l.y = dl_header->direction.y * (-1);
				l.z = dl_header->direction.z * (-1);
				l.computeModule();
				l.normalization();

				vector h;
				h.x = l.x + v.x;
				h.y = l.y + v.y;
				h.z = l.z + v.z;
				h.computeModule();
				h.normalization();

				if (dotProduct(v, tri.normalVector) < 0) {
					tri.normalVector.x = (-1) * tri.normalVector.x;
					tri.normalVector.y = (-1) * tri.normalVector.y;
					tri.normalVector.z = (-1) * tri.normalVector.z;
				}

				float nl_dotproduct, zero, nh_dotproduct;
				nl_dotproduct = dotProduct(tri.normalVector, l);
				nh_dotproduct = dotProduct(tri.normalVector, h);

				zero = 0;
				temp.r = temp.r + tri.nearestPoint->dr * dl_header->colors.r * 255 * max(zero, nl_dotproduct) * shadow + tri.nearestPoint->sr * dl_header->colors.r * 255 * pow(max(zero, nh_dotproduct), tri.nearestPoint->ns) * shadow;
				temp.g = temp.g + tri.nearestPoint->dg * dl_header->colors.g * 255 * max(zero, nl_dotproduct) * shadow + tri.nearestPoint->sg * dl_header->colors.g * 255 * pow(max(zero, nh_dotproduct), tri.nearestPoint->ns) * shadow;
				temp.b = temp.b + tri.nearestPoint->db * dl_header->colors.b * 255 * max(zero, nl_dotproduct) * shadow + tri.nearestPoint->sb * dl_header->colors.b * 255 * pow(max(zero, nh_dotproduct), tri.nearestPoint->ns) * shadow;
			}
			I.r = I.r + temp.r / 32;
			I.g = I.g + temp.g / 32;
			I.b = I.b + temp.b / 32;
		}

		if (depth < max_depth) {// recursion
			float nv_dotproduct = dotProduct(tri.normalVector, v);
			vector r;
			r.x = -v.x + 2 * nv_dotproduct * tri.normalVector.x;
			r.y = -v.y + 2 * nv_dotproduct * tri.normalVector.y;
			r.z = -v.z + 2 * nv_dotproduct * tri.normalVector.z;
			r.computeModule();
			depth++;
			R.direction = r;
			R.startPoint = intersect;
			color temp_I = recursiveRay(depth, R, vertex);// reflection
			I.r = I.r + temp_I.r * tri.nearestPoint->sr;
			I.g = I.g + temp_I.g * tri.nearestPoint->sg;
			I.b = I.b + temp_I.b * tri.nearestPoint->sb;

			vector T;
			float cosRefractionAngle = 1 - (1 - pow(nv_dotproduct, 2)) / pow(tri.nearestPoint->ior, 2);
			if (cosRefractionAngle >= 0) {
				cosRefractionAngle = sqrt(cosRefractionAngle);
				T.x = (nv_dotproduct / tri.nearestPoint->ior - cosRefractionAngle) * tri.normalVector.x - v.x / tri.nearestPoint->ior;
				T.y = (nv_dotproduct / tri.nearestPoint->ior - cosRefractionAngle) * tri.normalVector.y - v.y / tri.nearestPoint->ior;
				T.z = (nv_dotproduct / tri.nearestPoint->ior - cosRefractionAngle) * tri.normalVector.z - v.z / tri.nearestPoint->ior;
				T.computeModule();

				R.direction = T;
				R.startPoint = intersect;
				color temp_I = recursiveRay(depth, R, vertex);// refraction
				I.r = I.r + temp_I.r * tri.nearestPoint->tr;
				I.g = I.g + temp_I.g * tri.nearestPoint->tg;
				I.b = I.b + temp_I.b * tri.nearestPoint->tb;
			}
		}
	}

	return I;
}

Image* rayTracing(Image* img, coordinates* vertex) {
	int h = img->Height();
	int w = img->Width();
	float radian = cam.ha * 3.1416 / 180;
	float focalLength = 0.5 * h / tan(radian);
	srand(time(NULL));
	#pragma omp parallel for
	for (int i = 0; i < w; i++) {
		for (int j = 0; j < h; j++) {
			color totalcolor(0, 0, 0);
			int count = 0;
			for (int k = 0; k < 1; k++) {// loop is used for jittered sampling. The number of k depends on the number of samples
				float jittered_i, jittered_j, mod = 10;
				jittered_i = rand() % 10 / mod - 0.5;
				jittered_j = rand() % 10 / mod - 0.5;
				vector pixel;
				pixel.x = cam.d.x * focalLength + V.x * (h / 2 - j + 0) + U.x * (0 + i - w / 2);// pixel at assumed direction (0, 0, 1)
				pixel.y = cam.d.y * focalLength + V.y * (h / 2 - j + 0) + U.y * (0 + i - w / 2);
				pixel.z = cam.d.z * focalLength + V.z * (h / 2 - j + 0) + U.z * (0 + i - w / 2);
				ray R(cam.p, pixel);
				color temp;
				int depth = 0;
				temp = recursiveRay(depth, R, vertex);
				totalcolor.r = temp.r + totalcolor.r;
				totalcolor.g = temp.g + totalcolor.g;
				totalcolor.b = temp.b + totalcolor.b;
				count++;
			}
			totalcolor.r = totalcolor.r / count;
			totalcolor.g = totalcolor.g / count;
			totalcolor.b = totalcolor.b / count;
			Pixel p;
			p.r = rgbThreshold(totalcolor.r);
			p.g = rgbThreshold(totalcolor.g);
			p.b = rgbThreshold(totalcolor.b);
			img->SetPixel(i, j, p);
		}
	}
	return img;
}

sphereAndt sphereRayIntersection(ray R) {
	float a, b, c, delta, t, min;// quadratic equation
	min = NULL;
	
	sphereAndt sph;
	sph.t = -20;
	spheres* s = sphere_header;
	for (int spheres = 0; spheres < spheres_count; spheres++) {
		a = pow(R.direction.x, 2) + pow(R.direction.y, 2) + pow(R.direction.z, 2);
		b = 2 * (R.startPoint.x - s->center.x) * R.direction.x + 2 * (R.startPoint.y - s->center.y) * R.direction.y + 2 * (R.startPoint.z - s->center.z) * R.direction.z;
		c = pow((R.startPoint.x - s->center.x), 2) + pow((R.startPoint.y - s->center.y), 2) + pow((R.startPoint.z - s->center.z), 2) - pow(s->r, 2);
		delta = pow(b, 2) - 4 * a * c;
		t = (-b - sqrt(delta)) / (2 * a);
		if (min == NULL && t > 0) {
			min = t;
			sph.nearestPoint = s;
			sph.t = t;
		}
		else if (t < min && t > 0) {
			min = t;
			sph.nearestPoint = s;
			sph.t = t;
		}
		s = s->next;
	}

	return sph;
}

triangleAndt triangleRayIntersection(coordinates* vertex, ray R) {
	
	triangleAndt tri;
	tri.t = -20;
	triangle* t = triangle_header;
	float min = NULL;
	int min_check = 0;
	vector normalVector;// normal vector of the triangle plain 
	for (int k = 0; k < triangle_count; k++) {
		coordinates v1(vertex[t->vertex1].x, vertex[t->vertex1].y, vertex[t->vertex1].z);
		coordinates v2(vertex[t->vertex2].x, vertex[t->vertex2].y, vertex[t->vertex2].z);
		coordinates v3(vertex[t->vertex3].x, vertex[t->vertex3].y, vertex[t->vertex3].z);

		vector v1v2(v1, v2);
		vector v2v3(v2, v3);
		vector v3v1(v3, v1);

		normalVector = crossProduct(v1v2, v2v3);
		normalVector.normalization();
		
		float m = 0;
		coordinates p;// point where ray intersects with the triangle plain
		m = (normalVector.x * (v1.x - R.startPoint.x) + normalVector.y * (v1.y - R.startPoint.y) + normalVector.z * (v1.z - R.startPoint.z)) / (normalVector.x * R.direction.x + normalVector.y * R.direction.y + normalVector.z * R.direction.z);
		if (m > 0) {
			p.x = R.startPoint.x + R.direction.x * m;
			p.y = R.startPoint.y + R.direction.y * m;
			p.z = R.startPoint.z + R.direction.z * m;

			int inside = 1;
			vector v1p(v1, p);
			vector v1v3(v1, v3);
			if (dotProduct(crossProduct(v1p, v1v2), crossProduct(v1v3, v1v2)) < 0) {// see if p is inside triangle
				inside = 0;
			}
			vector v2p(v2, p);
			vector v2v1(v2, v1);
			if (dotProduct(crossProduct(v2p, v2v3), crossProduct(v2v1, v2v3)) < 0) {
				inside = 0;
			}
			vector v3p(v3, p);
			vector v3v2(v3, v2);
			if (dotProduct(crossProduct(v3p, v3v1), crossProduct(v3v2, v3v1)) < 0) {
				inside = 0;
			}
			/*cout << normalVector.x << ' ' << normalVector.y << ' ' << normalVector.z << endl;
			cout << m;*/
			
			if (inside == 1) {// shade
				//cout << normalVector.x << ' ' << normalVector.y << ' ' << normalVector.z << endl;
				if (min_check == 0) {
					min_check = 1;
					min = m;
					tri.nearestPoint = t;
					tri.t = m;
					tri.normalVector = normalVector;
				}
				else if (m < min && min_check == 1) {
					min = m;
					tri.nearestPoint = t;
					tri.t = m;
					//cout << normalVector.x << ' ' << normalVector.y << ' ' << normalVector.z << endl;
					tri.normalVector.x = normalVector.x;
					tri.normalVector.y = normalVector.y;
					tri.normalVector.z = normalVector.z;
					//cout << tri.normalVector.x << ' ' << tri.normalVector.y << ' ' << tri.normalVector.z << endl;
				}
			}
		}
		t = t->next;
	}
	//cout << tri.normalVector.x << ' ' << tri.normalVector.y << ' ' << tri.normalVector.z << endl;
	tri.t = tri.t * 1.0;
	return tri;
}