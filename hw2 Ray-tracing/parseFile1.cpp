//File parsing example

#include "stdafx.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include "hw2Raytracing.h"
#include "image.h"

#define STB_IMAGE_IMPLEMENTATION //only place once in one .cpp file
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION //only place once in one .cpp files
#include "stb_image_write.h"

using namespace std;

// camera default setting 
//float px = 0, py = 0, pz = 0, dx = 0, dy = 0, dz = 1, ux = 0, uy = 1, uz = 0, ha = 45;
camera cam(0, 0, 0, 0, 0, 1, 0, 1, 0, 45);
// film resolution default setting
int width = 640, height = 480;
// background color default setting
color background(0, 0, 0);
// material default setting
float ar = 0, ag = 0, ab = 0, dr = 1, dg = 1, db = 1, sr = 0, sg = 0, sb = 0, ns = 5, tr = 0, tg = 0, tb = 0, ior = 1;
// filename default setting
char outFile[1024] = "raytraced.bmp";
// ambient light default setting
color ambientLight(0, 0, 0);
// set default maximum recursion depth
int max_depth = 5;

int max_vertices = NULL;
int max_normals = NULL;
coordinates* vertex = NULL;
coordinates* normal = NULL;
int vertices_count = 0;

triangle *triangle_header = NULL;
int triangle_count = 0;
triangle *triangle_temp = NULL;

spheres *sphere_header = NULL;
int spheres_count = 0;
spheres *sphere_temp = NULL;

pointlights *pointlight_header = NULL;
int pointlight_count = 0;
pointlights *pointlight_temp = NULL;

directionalLight* dl_header = NULL;

spotlight* spotlight_header = NULL;
int spotlight_count = 0;
spotlight* spotlight_temp = NULL;

int main(){
  FILE *fp;
  long length;
  char line[1024]; //Assumes no line is longer than 1024 characters!

  string fileName = "dragon.scn";

  // open the file containing the scene description
  fp = fopen(fileName.c_str(), "r");

  // check for errors in opening the file
  if (fp == NULL) {
	  printf("Can't open file '%s'\n", fileName.c_str());
	  return 0;  //Exit
  }
	
  // determine the file size (this is optional -- feel free to delete the 4 lines below)
  fseek(fp, 0, SEEK_END); // move position indicator to the end of the file;
  length = ftell(fp);  // return the value of the current position
  printf("File '%s' is %ld bytes long.\n\n",fileName.c_str(),length);
  fseek(fp, 0, SEEK_SET);  // move position indicator to the start of the file
  
  Image* img = NULL;
  //Loop through reading each line
  while( fgets(line,1024,fp) ) { //Assumes no line is longer than 1024 characters!
    if (line[0] == '#'){
		printf("Skipping comment: %s", line);
		continue;
    }
    
    char command[100];
    int fieldsRead = sscanf(line,"%s ",command); //Read first word in the line (i.e., the command type)
    
    if (fieldsRead < 1){ //No command read
		//Blank line
		continue;
    }
    
	if (strcmp(command, "triangle") == 0) { //If the command is a triangle command
		triangle *p = new triangle;
		if (triangle_count == 0) {
			triangle_header = p;
		}
		else {
			triangle_temp->next = p;
		}
		triangle_count++;
		triangle_temp = p;
		p->next = NULL;
		sscanf(line, "triangle %d %d %d", &p->vertex1, &p->vertex2, &p->vertex3);
		printf("Triangle with vertices %d %d %d\n", p->vertex1, p->vertex2, p->vertex3);
		p->ar = ar; p->ag = ag; p->ab = ab;
		p->dr = dr; p->dg = dg; p->db = db;
		p->sr = sr; p->sg = sg; p->sb = sb; p->ns = ns;
		p->tr = tr; p->tg = tg; p->tb = tb; p->ior = ior;
	}
    else if (strcmp(command, "sphere")==0){ //If the command is a sphere command
		spheres *p = new spheres;
		if (spheres_count == 0) {
			sphere_header = p;
		}
		else {
			sphere_temp->next = p;
		}
		spheres_count++;
		sphere_temp = p;
		p->next = NULL;
		sscanf(line, "sphere %f %f %f %f", &p->center.x, &p->center.y, &p->center.z, &p->r);
		printf("Sphere as position (%f,%f,%f) with radius %f\n", p->center.x, p->center.y, p->center.z, p->r);
		p->ar = ar; p->ag = ag; p->ab = ab;
		p->dr = dr; p->dg = dg; p->db = db;
		p->sr = sr; p->sg = sg; p->sb = sb; p->ns = ns;
		p->tr = tr; p->tg = tg; p->tb = tb; p->ior = ior;
    }
	else if (strcmp(command, "material") == 0) { //If the command is a material command
		sscanf(line, "material %f %f %f %f %f %f %f %f %f %f %f %f %f %f", &ar, &ag, &ab, &dr, &dg, &db, &sr, &sg, &sb, &ns, &tr, &tg, &tb, &ior);
		printf("Material: %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", ar, ag, ab, dr, dg, db, sr, sg, sb, ns, tr, tg, tb, ior);
	}
    else if (strcmp(command, "background")==0){ //If the command is a background command
		sscanf(line, "background %f %f %f", &background.r, &background.g, &background.b);
		printf("Background color of (%f,%f,%f)\n", background.r, background.g, background.b);
    }
	else if (strcmp(command, "ambient_light") == 0) { //If the command is an ambient_light command
		sscanf(line, "ambient_light %f %f %f", &ambientLight.r, &ambientLight.g, &ambientLight.b);
		printf("Ambient_light color of (%f,%f,%f)\n", ambientLight.r, ambientLight.g, ambientLight.b);
	}
    else if (strcmp(command, "output_image")==0){ //If the command is an output_image command
		sscanf(line,"output_image %s", outFile);
		printf("Render to file named: %s\n", outFile);
    }
	else if (strcmp(command, "film_resolution") == 0) { //If the command is a film_resolution command
		sscanf(line, "film_resolution %d %d", &width, &height);
		printf("Film_resolution of (%d,%d)\n", width, height);
		//setImageNull(img);
	}
	else if (strcmp(command, "camera") == 0) { //If the command is a camera command
		sscanf(line, "camera %f %f %f %f %f %f %f %f %f %f", &cam.p.x, &cam.p.y, &cam.p.z, &cam.d.x, &cam.d.y, &cam.d.z, &cam.u.x, &cam.u.y, &cam.u.z, &cam.ha);
		printf("Camera: %f %f %f %f %f %f %f %f %f %f\n", cam.p.x, cam.p.y, cam.p.z, cam.d.x, cam.d.y, cam.d.z, cam.u.x, cam.u.y, cam.u.z, cam.ha);
		cam.d.computeModule();
		cam.u.computeModule();
		setRelativeCoordinate();
	}
	else if (strcmp(command, "point_light") == 0) { //If the command is a point_light command
		pointlights *p = new pointlights;
		if (pointlight_count == 0) {
			pointlight_header = p;
		}
		else {
			pointlight_temp->next = p;
		}
		pointlight_count++;
		pointlight_temp = p;
		p->next = NULL;
		sscanf(line, "point_light %f %f %f %f %f %f ", &p->colors.r, &p->colors.g, &p->colors.b, &p->center.x, &p->center.y, &p->center.z);
		printf("Point light: %f %f %f %f %f %f \n", p->colors.r, p->colors.g, p->colors.b, p->center.x, p->center.y, p->center.z);
	}
	else if (strcmp(command, "spot_light") == 0) { //If the command is a spot_light command
		spotlight *p = new spotlight;
		if (spotlight_count == 0) {
			spotlight_header = p;
		}
		else {
			spotlight_temp->next = p;
		}
		spotlight_count++;
		spotlight_temp = p;
		p->next = NULL;
		sscanf(line, "spot_light %f %f %f %f %f %f %f %f %f %f %f", &p->colors.r, &p->colors.g, &p->colors.b, &p->center.x, &p->center.y, &p->center.z, &p->centralline.x, &p->centralline.y, &p->centralline.z, &p->angle1, &p->angle2);
		printf("Spot light: %f %f %f %f %f %f %f %f %f %f %f\n", p->colors.r, p->colors.g, p->colors.b, p->center.x, p->center.y, p->center.z, p->centralline.x, p->centralline.y, p->centralline.z, p->angle1, p->angle2);
	}
	else if (strcmp(command, "directional_light") == 0) { //If the command is a directional_light command
		dl_header = new directionalLight;
		sscanf(line, "directional_light %f %f %f %f %f %f ", &dl_header->colors.r, &dl_header->colors.g, &dl_header->colors.b, &dl_header->direction.x, &dl_header->direction.y, &dl_header->direction.z);
		printf("Point light: %f %f %f %f %f %f \n", dl_header->colors.r, dl_header->colors.g, dl_header->colors.b, dl_header->direction.x, dl_header->direction.y, dl_header->direction.z);
	}
	else if (strcmp(command, "max_vertices") == 0) { //If the command is a max_vertices command
		sscanf(line, "max_vertices %d", &max_vertices);
		printf("Max_vertices: %d\n", max_vertices);
		vertex = new coordinates[max_vertices];
	}
	else if (strcmp(command, "vertex") == 0) { //If the command is a vertex command
		float vertex_x, vertex_y, vertex_z;
		sscanf(line, "vertex %f %f %f", &vertex_x, &vertex_y, &vertex_z);
		printf("Vertex: %f %f %f\n", vertex_x, vertex_y, vertex_z);
		vertex[vertices_count].x = vertex_x;
		vertex[vertices_count].y = vertex_y;
		vertex[vertices_count].z = vertex_z;
		vertices_count++;
	}
	else if (strcmp(command, "max_depth") == 0) { //If the command is a max_depth command
		sscanf(line, "max_depth %d", &max_depth);
		printf("Max_depth: %d\n", max_depth);
	}
    else {
		printf("WARNING. Do not know command: %s\n", command);
    }
  }
  setRelativeCoordinate();
  img = new Image(width, height);// create film
  //setBackgroundColor(r, g, b, img);// set background color

  //sphere_header = sphereListing(sphere_header, spheres_count);
  rayTracing(img, vertex);
  img->Write(outFile);

  return 0;
}
