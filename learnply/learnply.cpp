#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <iostream>
#include <algorithm>

#include "glError.h"
#include "gl/glew.h"
#include "gl/freeglut.h"
#include "ply.h"
#include "icVector.H"
#include "icMatrix.H"
#include "polyhedron.h"
#include "trackball.h"
#include "tmatrix.h"

#include "drawUtil.h"

Polyhedron* poly;
std::vector<PolyLine> streamlines;
std::vector<icVector3> points;
std::vector<icVector3> sources;
std::vector<icVector3> sinks;
std::vector<icVector3> saddles;
std::vector<icVector3> higher_order;
const double STEP = .1; // You should experiment to find the optimal step size.
const int STEP_MAX = 1000; // Upper limit of steps to take for tracing each streamline.
std::vector<PolyLine> lines; // Used for storing streamlines.

/*scene related variables*/
const float zoomspeed = 0.9;
int win_width = 1024;
int win_height = 1024;
float aspectRatio = win_width / win_height;
const int view_mode = 0;		// 0 = othogonal, 1=perspective
const double radius_factor = 0.9;

/*
Use keys 1 to 0 to switch among different display modes.
Each display mode can be designed to show one type 
visualization result.

Predefined ones: 
display mode 1: solid rendering
display mode 2: show wireframes
display mode 3: render each quad with colors of vertices
*/
int display_mode = 1;

/*User Interaction related variabes*/
float s_old, t_old;
float rotmat[4][4];
double zoom = 1.0;
double translation[2] = { 0, 0 };
int mouse_mode = -2;	// -1 = no action, 1 = tranlate y, 2 = rotate

// IBFV related variables (Van Wijk 2002)
//https://www.win.tue.nl/~vanwijk/ibfv/
#define NPN		64
#define SCALE	4.0
#define ALPHA	8
float tmax = win_width / (SCALE * NPN);
float dmax = SCALE / win_width;
unsigned char* pixels;

/******************************************************************************
Forward declaration of functions
******************************************************************************/

void init(void);
void initIBFV();
void find_singularities();
void classify_singularities();
bool onEdge(icVector2, Quad*);
Quad* streamline_step(icVector2&, icVector2&, Quad*, bool);
PolyLine build_streamline(double, double);

/*glut attaching functions*/
void keyboard(unsigned char key, int x, int y);
void motion(int x, int y);
void displayIBFV();
void display(void);
void mouse(int button, int state, int x, int y);
void mousewheel(int wheel, int direction, int x, int y);
void reshape(int width, int height);

/*functions for element picking*/
void display_vertices(GLenum mode, Polyhedron* poly);
void display_quads(GLenum mode, Polyhedron* poly);
void display_selected_vertex(Polyhedron* poly);
void display_selected_quad(Polyhedron* poly);

void display3DIBFV();
void displayIBFVSlices(Polyhedron**);
void drawQuad();

/*display vis results*/
void display_polyhedron(Polyhedron* poly);

/******************************************************************************
Main program.
******************************************************************************/

int main(int argc, char* argv[])
{
	/*load mesh from ply file*/
	FILE* this_file = fopen("../data/3D/rotstrat4096.ply", "r");
	poly = new Polyhedron(this_file);
	fclose(this_file);
	
	/*initialize the mesh*/
	poly->initialize(); // initialize the mesh
	poly->write_info();


	/*init glut and create window*/
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowPosition(20, 20);
	glutInitWindowSize(win_width, win_height);
	glutCreateWindow("Scientific Visualization");


	/*initialize openGL*/
	init();
	
	/*the render function and callback registration*/
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutMotionFunc(motion);
	glutMouseFunc(mouse);
	glutMouseWheelFunc(mousewheel);
	
	/*event processing loop*/
	glutMainLoop();
	
	/*clear memory before exit*/
	poly->finalize();	// finalize everything
	free(pixels);
	return 0;
}

/******************************************************************************
Set projection mode
******************************************************************************/

void set_view(GLenum mode)
{
	GLfloat light_ambient0[] = { 0.3, 0.3, 0.3, 1.0 };
	GLfloat light_diffuse0[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat light_specular0[] = { 0.0, 0.0, 0.0, 1.0 };

	GLfloat light_ambient1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_diffuse1[] = { 0.5, 0.5, 0.5, 1.0 };
	GLfloat light_specular1[] = { 0.0, 0.0, 0.0, 1.0 };

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);

	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);


	glMatrixMode(GL_PROJECTION);
	if (mode == GL_RENDER)
		glLoadIdentity();

	if (aspectRatio >= 1.0) {
		if (view_mode == 0)
			glOrtho(-radius_factor * zoom * aspectRatio, radius_factor * zoom * aspectRatio, -radius_factor * zoom, radius_factor * zoom, -1000, 1000);
		else
			glFrustum(-radius_factor * zoom * aspectRatio, radius_factor * zoom * aspectRatio, -radius_factor * zoom, radius_factor * zoom, 0.1, 1000);
	}
	else {
		if (view_mode == 0)
			glOrtho(-radius_factor * zoom, radius_factor * zoom, -radius_factor * zoom / aspectRatio, radius_factor * zoom / aspectRatio, -1000, 1000);
		else
			glFrustum(-radius_factor * zoom, radius_factor * zoom, -radius_factor * zoom / aspectRatio, radius_factor * zoom / aspectRatio, 0.1, 1000);
	}

	GLfloat light_position[3];
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	light_position[0] = 5.5;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	light_position[0] = -0.1;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT2, GL_POSITION, light_position);
}

/******************************************************************************
Update the scene
******************************************************************************/

void set_scene(GLenum mode, Polyhedron* poly)
{
	glTranslatef(translation[0], translation[1], -3.0);

	/*multiply rotmat to current mat*/
	{
		int i, j, index = 0;

		GLfloat mat[16];

		for (i = 0; i < 4; i++)
			for (j = 0; j < 4; j++)
				mat[index++] = rotmat[i][j];

		glMultMatrixf(mat);
	}

	glScalef(0.9 / poly->radius, 0.9 / poly->radius, 0.9 / poly->radius);
	glTranslatef(-poly->center.entry[0], -poly->center.entry[1], -poly->center.entry[2]);
}

/******************************************************************************
Init scene
******************************************************************************/

void init(void) {

	mat_ident(rotmat);

	/* select clearing color */
	glClearColor(0.0, 0.0, 0.0, 0.0);  // background
	glShadeModel(GL_FLAT);
	glPolygonMode(GL_FRONT, GL_FILL);

	glDisable(GL_DITHER);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	
	//set pixel storage modes
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	
	glEnable(GL_NORMALIZE);
	if (poly->orientation == 0)
		glFrontFace(GL_CW);
	else
		glFrontFace(GL_CCW);
}

/******************************************************************************
Initialize IBFV patterns
******************************************************************************/

void initIBFV()
{
	pixels = (unsigned char*)malloc(sizeof(unsigned char) * win_width * win_height * 3);
	memset(pixels, 255, sizeof(unsigned char) * win_width * win_height * 3);

	tmax = win_width / (SCALE * NPN);
	dmax = SCALE / win_width;

	int lut[256];
	int phase[NPN][NPN];
	GLubyte pat[NPN][NPN][4];
	int i, j, k;

	for (i = 0; i < 256; i++) lut[i] = i < 127 ? 0 : 255;
	for (i = 0; i < NPN; i++)
		for (j = 0; j < NPN; j++) phase[i][j] = rand() % 256;

	for (i = 0; i < NPN; i++)
	{
		for (j = 0; j < NPN; j++)
		{
			pat[i][j][0] =
				pat[i][j][1] =
				pat[i][j][2] = lut[(phase[i][j]) % 255];
			pat[i][j][3] = ALPHA;
		}
	}
	
	glNewList(1, GL_COMPILE);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0, GL_RGBA, GL_UNSIGNED_BYTE, pat);
	glEndList();
}

/******************************************************************************
Pick objects from the scene
******************************************************************************/

int processHits(GLint hits, GLuint buffer[])
{
	unsigned int i, j;
	GLuint names, * ptr;
	double smallest_depth = 1.0e+20, current_depth;
	int seed_id = -1;
	unsigned char need_to_update;

	ptr = (GLuint*)buffer;
	for (i = 0; i < hits; i++) {  /* for each hit  */
		need_to_update = 0;
		names = *ptr;
		ptr++;

		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		for (j = 0; j < names; j++) {  /* for each name */
			if (need_to_update == 1)
				seed_id = *ptr - 1;
			ptr++;
		}
	}
	return seed_id;
}

/******************************************************************************
Diaplay all quads for selection
******************************************************************************/

void display_quads(GLenum mode, Polyhedron* this_poly)
{
	unsigned int i, j;
	GLfloat mat_diffuse[4];

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glDisable(GL_LIGHTING);

	for (i = 0; i < this_poly->nquads; i++) {
		if (mode == GL_SELECT)
			glLoadName(i + 1);

		Quad* temp_q = this_poly->qlist[i];
		
		glBegin(GL_POLYGON);
		for (j = 0; j < 4; j++) {
			Vertex* temp_v = temp_q->verts[j];
			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
		}
		glEnd();
	}
}

/******************************************************************************
Diaplay all vertices for selection
******************************************************************************/

void display_vertices(GLenum mode, Polyhedron* this_poly)
{
	for (int i = 0; i < this_poly->nverts; i++) {
		if (mode == GL_SELECT)
			glLoadName(i + 1);

		CHECK_GL_ERROR();

		Vertex* temp_v = this_poly->vlist[i];
		drawDot(temp_v->x, temp_v->y, temp_v->z, 0.15);
	}
	CHECK_GL_ERROR();
}

/******************************************************************************
Diaplay selected quad
******************************************************************************/

void display_selected_quad(Polyhedron* this_poly)
{
	if (this_poly->selected_quad == -1)
	{
		return;
	}

	unsigned int i, j;

	glDisable(GL_POLYGON_OFFSET_FILL);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glDisable(GL_LIGHTING);

	Quad* temp_q = this_poly->qlist[this_poly->selected_quad];

	glBegin(GL_POLYGON);
	for (j = 0; j < 4; j++) {
		Vertex* temp_v = temp_q->verts[j];
		glColor3f(1.0, 0.0, 0.0);
		glVertex3d(temp_v->x, temp_v->y, temp_v->z);
	}
	glEnd();
}

/******************************************************************************
Diaplay selected vertex
******************************************************************************/

void display_selected_vertex(Polyhedron* this_poly)
{
	if (this_poly->selected_vertex == -1)
	{
		return;
	}

	Vertex* temp_v = this_poly->vlist[this_poly->selected_vertex];
	drawDot(temp_v->x, temp_v->y, temp_v->z, 0.15, 1.0, 0.0,0.0);

	CHECK_GL_ERROR();
}

/******************************************************************************
Process a keyboard action.  In particular, exit the program when an
"escape" is pressed in the window.
******************************************************************************/


void display3DIBFV() {
	//n = number of slices in S_i?
	int N = 16;
	//std::vector<Polyhedron*> S;
	Polyhedron* S[16];
	for (int i = 0; i < 16; i++) {
		//S[i]->vlist = poly->vlist;
		//Polyhedron* S_temp = new Polyhedron(poly->vlist, 256, i*256);
		//S.push_back(S_temp);
		S[i] = new Polyhedron(poly->vlist, 256, i * 256);
		//delete(S_temp);
		//std::cout << S[i]->vlist[0]->vz << std::endl;
	}

	//GL_RGBA A[16];

	GLsizei size = N;
	GLuint *A = new GLuint[N];
	glGenTextures(size, A);

	for (int k = 0; k < N; k++) {
		for (int i = 1; i < S[k]->nverts; i++) {
			if (k > 0) {

				float vzk = std::max(S[k - 1]->vlist[i]->vz, 0.0);
				//Do 1D Z-axis advection from Si-1 to Si
				glDisable(GL_BLEND);
				glBindTexture(GL_TEXTURE_2D, A[k - 1]); drawQuad();

				glEnable(GL_BLEND);
				glBlendFunc(GL_ZERO, GL_SRC_COLOR);
				glBindTexture(GL_TEXTURE_2D, vzk ); drawQuad();

				GLuint temp;
				glGenTextures(1, &temp);
				glBindTexture(GL_TEXTURE_2D, temp);
				glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 0, 0, N, N, 0);

				glDisable(GL_BLEND);
				glBindTexture(GL_TEXTURE_2D, A[k]); drawQuad();

				glEnable(GL_BLEND);
				glBlendFunc(GL_ZERO, GL_ONE_MINUS_SRC_COLOR);
				glBindTexture(GL_TEXTURE_2D, vzk ); drawQuad();

				glBlendFunc(GL_ONE, GL_ONE);
				glBindTexture(GL_TEXTURE_2D, temp); drawQuad();

				glBindTexture(GL_TEXTURE_2D, A[k]);
				glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 0, 0, N, N, 0);

			}
			if (k < N - 1) {
				//Do 1D Z-axis advection from Si+1 to Si

				float vzk = std::max(-1 * S[k + 1]->vlist[i]->vz, 0.0);
				glDisable(GL_BLEND);
				glBindTexture(GL_TEXTURE_2D, A[k + 1]); drawQuad();

				glEnable(GL_BLEND);
				glBlendFunc(GL_ZERO, GL_SRC_COLOR);
				glBindTexture(GL_TEXTURE_2D, vzk); drawQuad();

				GLuint temp;
				glGenTextures(1, &temp);

				glBindTexture(GL_TEXTURE_2D, temp);
				glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 0, 0, N, N, 0);

				glDisable(GL_BLEND);
				glBindTexture(GL_TEXTURE_2D, A[k]); drawQuad();

				glEnable(GL_BLEND);
				glBlendFunc(GL_ZERO, GL_ONE_MINUS_SRC_COLOR);
				glBindTexture(GL_TEXTURE_2D, vzk ); drawQuad();

				glBlendFunc(GL_ONE, GL_ONE);
				glBindTexture(GL_TEXTURE_2D, temp); drawQuad();

				glBindTexture(GL_TEXTURE_2D, A[k]);
				glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 0, 0, N, N, 0);
				
			}
			//Do 2d IBFV-based advection in the slice Si
		}
	}
	displayIBFVSlices(S);
	

	for (int i = 0; i < 16; i++) {
		delete(S[i]);
	}

}

void drawQuad() {
	//This function does not seem to actually do anything rihgt now
	GLenum mode;
	//draws a single textured quadrilateral that covers the whole image.
	unsigned int i, j;
	GLfloat mat_diffuse[4];

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glDisable(GL_LIGHTING);

	glBegin(GL_POLYGON);
	glVertex3d(10, 10, 0);
	glVertex3d(10, -10, 0);
	glVertex3d(-10, -10, 0);
	glVertex3d(-10, 10, 0);
	glEnd();
}


void displayIBFVSlices(Polyhedron** S) {
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_BLEND);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

	glEnable(GL_TEXTURE_2D);
	glShadeModel(GL_FLAT);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glClearColor(0.5, 0.5, 0.5, 1.0);  // background for rendering color coding and lighting
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// draw the mesh using pixels and use vector field to advect texture coordinates
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, win_width, win_height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels);

	double modelview_matrix[16], projection_matrix[16];
	int viewport[4];
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix);
	glGetDoublev(GL_PROJECTION_MATRIX, projection_matrix);
	glGetIntegerv(GL_VIEWPORT, viewport);
	for (int k = 0; k < 16; k++) {
		for (int i = 0; i < S[k]->nquads; i++)
		{
			Quad* qtemp = S[k]->qlist[i];

			glBegin(GL_QUADS);
			for (int j = 0; j < 4; j++)
			{
				Vertex* vtemp = qtemp->verts[j];

				double tx, ty, dummy;
				gluProject((GLdouble)vtemp->x, (GLdouble)vtemp->y, (GLdouble)vtemp->z,
					modelview_matrix, projection_matrix, viewport, &tx, &ty, &dummy);

				tx = tx / win_width;
				ty = ty / win_height;

				icVector2 dp = icVector2(vtemp->vx, vtemp->vy);
				normalize(dp);
				dp *= dmax;

				double dx = -dp.x;
				double dy = -dp.y;

				float px = tx + dx;
				float py = ty + dy;

				glTexCoord2f(px, py);
				glVertex3d(vtemp->x, vtemp->y, vtemp->z);
			}
			glEnd();
		}
	}

	glEnable(GL_BLEND);

	// blend in noise pattern
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glTranslatef(-1.0, -1.0, 0.0);
	glScalef(2.0, 2.0, 1.0);

	glCallList(1);

	glBegin(GL_QUAD_STRIP);

	glTexCoord2f(0.0, 0.0);  glVertex2f(0.0, 0.0);
	glTexCoord2f(0.0, tmax); glVertex2f(0.0, 1.0);
	glTexCoord2f(tmax, 0.0);  glVertex2f(1.0, 0.0);
	glTexCoord2f(tmax, tmax); glVertex2f(1.0, 1.0);
	glEnd();
	glDisable(GL_BLEND);

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glReadPixels(0, 0, win_width, win_height, GL_RGB, GL_UNSIGNED_BYTE, pixels);

	// draw the mesh using pixels without advecting texture coords
	glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, win_width, win_height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels);
	for (int k = 0; k < 16; k++) {
		for (int i = 0; i < S[k]->nquads; i++)
		{
			Quad* qtemp = S[k]->qlist[i];
			glBegin(GL_QUADS);
			for (int j = 0; j < 4; j++)
			{
				Vertex* vtemp = qtemp->verts[j];
				double tx, ty, dummy;
				gluProject((GLdouble)vtemp->x, (GLdouble)vtemp->y, (GLdouble)vtemp->z,
					modelview_matrix, projection_matrix, viewport, &tx, &ty, &dummy);
				tx = tx / win_width;
				ty = ty / win_height;
				glTexCoord2f(tx, ty);
				glVertex3d(vtemp->x, vtemp->y, vtemp->z);
			}
			glEnd();
		}
	}
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);
	glShadeModel(GL_SMOOTH);
}

void keyboard(unsigned char key, int x, int y) {
	int i;

	// clear out lines and points
	lines.clear();
	points.clear();

	switch (key) {
	case 27:	// set excape key to exit program
		poly->finalize();  // finalize_everything
		exit(0);
		break;

	case '1':	// solid color display with lighting
		display_mode = 1;
		glutPostRedisplay();
		break;

	case '2':	// wireframe display
		display_mode = 2;
		glutPostRedisplay();
		break;

	case '3':	// checkerboard display
	{
		display_mode = 3;

		double L = (poly->radius * 2) / 30;
		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			for (int j = 0; j < 4; j++) {

				Vertex* temp_v = temp_q->verts[j];

				temp_v->R = int(temp_v->x / L) % 2 == 0 ? 1 : 0;
				temp_v->G = int(temp_v->y / L) % 2 == 0 ? 1 : 0;
				temp_v->B = 0.0;
			}
		}
		glutPostRedisplay();
	}
	break;

	case '4':	// Drawing points and lines created by the dots_and_lines_example() function
		display_mode = 4;
		dots_and_lines_example(&points, &lines);
		glutPostRedisplay();
		break;

	case '5':	// IBFV vector field display
		display_mode = 5;
		glutPostRedisplay();
		break;

	case '6':	// add your own display mode
		display_mode = 6;
		{

			find_singularities();
			classify_singularities();
			
		}
		glutPostRedisplay();
		break;

	case '7':	// add your own display mode
		display_mode = 7;
		{
			double wanted_x = 0, wanted_y = 0;
			std::cout << "Enter an X and Y coordinate for the seed point of the streamline" << std::endl;
			std::cout << "X: ";
			std::cin >> wanted_x;
			std::cout << "Y: ";
			std::cin >> wanted_y;
			streamlines.push_back(build_streamline(wanted_x, wanted_y));

		}
		glutPostRedisplay();
		break;
	
	case '8':	// add your own display mode
		display_mode = 7;
		{
			double wanted_x = 0, wanted_y = 0;
			while (abs(wanted_x) <= 10 && abs(wanted_y) <= 10) {
				std::cout << "Enter an X and Y coordinate for the seed point of the streamline. Enter a point outside the bounds to stop entering new points" << std::endl;
				std::cout << "X: ";
				std::cin >> wanted_x;
				std::cout << "Y: ";
				std::cin >> wanted_y;
				streamlines.push_back(build_streamline(wanted_x, wanted_y));
			}

		}
		glutPostRedisplay();
		break;

	case '9':	// solid color display with lighting
		display_mode = 9;
		glutPostRedisplay();
		break;
	case '0':	// solid color display with lighting
		display_mode = 0;
		glutPostRedisplay();
		break;



	case 'r':	// reset rotation and transformation
		mat_ident(rotmat);
		translation[0] = 0;
		translation[1] = 0;
		zoom = 1.0;
		glutPostRedisplay();
		break;
		streamlines.clear();
	}
}


/******************************************************************************
Function for finding singularities of a vector set
******************************************************************************/
void find_singularities()
{
	for (int i = 0; i < poly->nquads; i++)
	{
		Quad* quad = poly->qlist[i];
		// 1. Delete any old singularity points
		
		delete quad->singularity;

		// 2. Find x1, x2, y1, y2, fx1y1, fx2y1, fx1y2, fx2y2, gx1y1, gx2y1, gx1y2, gx2y2
		// x1 is the smallest x coordinate of the 4 vertices
		// x2 is the largest x coordinate of the 4 vertices
		// y1 is the smallest y coordinate of the 4 vertices
		// y2 is the largest y coordinate of the 4 vertices
		// fx1y1 is the vector value x component at the vertex with coordinates (x1,y1)
		// gx1y1 is the vector value y component at the vertex with coordinates (x1,y1)
		// ...
		// to access the vector components of a Vertex* object, use vert->vx and vert->vy
		
		double smallestX = quad->verts[0]->x;
		double largestX = quad->verts[0]->x;
		double smallestY = quad->verts[0]->y;
		double largestY = quad->verts[0]->y;


		for (int j = 0; j < 4; j++) {
			if (smallestX > quad->verts[j]->x) {
				smallestX = quad->verts[j]->x;
			}
			if (largestX < quad->verts[j]->x) {
				largestX = quad->verts[j]->x;
			}
			if (smallestY > quad->verts[j]->y) {
				smallestY = quad->verts[j]->y;
			}
			if (largestY < quad->verts[j]->y) {
				largestY = quad->verts[j]->y;
			}
		}

		double fx1y1 = 0, fx1y2 = 0, fx2y1 = 0, fx2y2 = 0, gx1y1 = 0, gx1y2 = 0, gx2y1 = 0, gx2y2 = 0;

		for (int j = 0; j < 4; j++) {
			if (quad->verts[j]->x == smallestX && quad->verts[j]->y == smallestY) {
				fx1y1 = quad->verts[j]->vx;
				gx1y1 = quad->verts[j]->vy;
			}
			if (quad->verts[j]->x == smallestX && quad->verts[j]->y == largestY) {
				fx1y2 = quad->verts[j]->vx;
				gx1y2 = quad->verts[j]->vy;
			}
			if (quad->verts[j]->x == largestX && quad->verts[j]->y == smallestY) {
				fx2y1 = quad->verts[j]->vx;
				gx2y1 = quad->verts[j]->vy;
			}
			if (quad->verts[j]->x == largestX && quad->verts[j]->y == largestY) {
				fx2y2 = quad->verts[j]->vx;
				gx2y2 = quad->verts[j]->vy;
			}
		}


		// 3. compute the coefficients for solving the quadratic equation
		double a00 = fx1y1;
		double a10 = fx2y1 - fx1y1;
		double a01 = fx1y2 - fx1y1;
		double a11 = fx1y1 - fx2y1 - fx1y2 + fx2y2;
		double b00 = gx1y1;
		double b10 = gx2y1 - gx1y1;
		double b01 = gx1y2 - gx1y1;
		double b11 = gx1y1 - gx2y1 - gx1y2 + gx2y2;
		double c00 = a11 * b00 - a00 * b11;
		double c10 = a11 * b10 - a10 * b11;
		double c01 = a11 * b01 - a01 * b11;
		
		// 4. Compute the coefficients of the quadratic equation about s:
		// (-a11*c10)s2 + (-a11*c00 - a01*c10 + a10*c01)s + (a00*c01 - a01*c00) = 0.
		
		double a = (-a11 * c10);
		double b = (-a11 * c00 - a01 * c10 + a10 * c01);
		double c = (a00 * c01 - a01 * c00);
		
		// 5. Use the quadratic formula to solve for the s.
		// (check beforehand for complex values or dividing by zero)
		// You will get two values for s because of the ±.
		double sPositive = 0;
		double sNegative = 0;
		if (a != 0) {
			double insideSqrt = b * b - 4 * a * c;
			if (insideSqrt >= 0) {
				sPositive = ((-1 * b) + sqrt(insideSqrt)) / (2 * a);
				sNegative = ((-1 * b) - sqrt(insideSqrt)) / (2 * a);
			}
		}
		
		// 6. Use both values of s to get two corresponding values for t:
		// t = -(c00/c01) - (c10/c01)s
		 
		double  tPositive = -(c00 / c01) - (c10 / c01) * sPositive;
		double  tNegative = -(c00 / c01) - (c10 / c01) * sNegative;

		// 7. For each (s,t) pair, check that both values are between 0 and 1.
		// Either one or none of these pairs will satisfy this condition.
		// If one of the pairs has both components between 0 and 1,
		// then it corresponds to a singularity inside the quad.
		
		double s = -1;
		double t = -1;

		if (sPositive > 0 && sPositive < 1 && tPositive > 0 && tPositive < 1) {
			s = sPositive;
			t = tPositive;
		}

		if (sNegative > 0 && sNegative < 1 && tNegative > 0 && tNegative < 1) {
			s = sNegative;
			t = tNegative;
		}

		if (s != -1 && t != -1) {
			// 8. Compute the coordinates of the singularity inside the quad using s and t
			// use s to interpolate between x1 and x2 (s tells you how far inbetween
			// x1 and x2 the x coordinate is).
			// use t to interpolate between y1 and y2 (t tells you how far inbetween
			// y1 and y2 the y coordinate is).
			 
			double x = smallestX + s * (largestX - smallestX);
			double y = smallestY + t * (largestY - smallestY);

			// 9. Insert the singularity into the quad data structure
			// quad->singularity = new icVector2(x,y);
			// You will need to create a new field in the Quad data structure to store
			// the singularity point. The Quad class definition is located in the 
			// polyhedron.h file.
			
			quad->singularity = new icVector2(x, y);
		}
	}
}


/******************************************************************************
Function for classifying singularities that have already been found in a vector data set
******************************************************************************/
void classify_singularities()
{
	sources.clear();
	saddles.clear();
	higher_order.clear();
	for (int i = 0; i < poly->nquads; i++)
	{
		Quad* quad = poly->qlist[i];
		// 1. Check if the quad has a non-null singularity pointer. If not, continue.
		
		if (quad->singularity == NULL) {
			continue;
		}

		// 2. Find x1, x2, y1, y2, fx1y1, fx2y1, fx1y2, fx2y2, gx1y1, gx2y1, gx1y2, gx2y2
		// This is the same as in the find_singularities function
		
		double smallestX = quad->verts[0]->x;
		double largestX = quad->verts[0]->x;
		double smallestY = quad->verts[0]->y;
		double largestY = quad->verts[0]->y;

		for (int j = 0; j < 4; j++) {
			if (smallestX > quad->verts[j]->x) {
				smallestX = quad->verts[j]->x;
			}
			if (largestX < quad->verts[j]->x) {
				largestX = quad->verts[j]->x;
			}
			if (smallestY > quad->verts[j]->y) {
				smallestY = quad->verts[j]->y;
			}
			if (largestY < quad->verts[j]->y) {
				largestY = quad->verts[j]->y;
			}
		}

		double fx1y1 = 0, fx1y2 = 0, fx2y1 = 0, fx2y2 = 0, gx1y1 = 0, gx1y2 = 0, gx2y1 = 0, gx2y2 = 0;

		for (int j = 0; j < 4; j++) {
			if (quad->verts[j]->x == smallestX && quad->verts[j]->y == smallestY) {
				fx1y1 = quad->verts[j]->vx;
				gx1y1 = quad->verts[j]->vy;
			}
			if (quad->verts[j]->x == smallestX && quad->verts[j]->y == largestY) {
				fx1y2 = quad->verts[j]->vx;
				gx1y2 = quad->verts[j]->vy;
			}
			if (quad->verts[j]->x == largestX && quad->verts[j]->y == smallestY) {
				fx2y1 = quad->verts[j]->vx;
				gx2y1 = quad->verts[j]->vy;
			}
			if (quad->verts[j]->x == largestX && quad->verts[j]->y == largestY) {
				fx2y2 = quad->verts[j]->vx;
				gx2y2 = quad->verts[j]->vy;
			}
		}

		// 3. Get the coordinates of the singularity
		
		double x = quad->singularity->x;
		double y = quad->singularity->y;
		
		// 4. calculate the jacobian values dfdx, dfdy, dgdx, dgdy
		// dfdx = (-(y2 - y0) * fx1y1 + (y2 - y0) * fx2y1 - (y0 - y1) * fx1y2 + (y0 - y1) * fx2y2) / (x2 - x1) * (y2 - y1);
		// dfdy = (-(x2 - x0) * fx1y1 - (x0 - x1) * fx2y1 + (x2 - x0) * fx1y2 + (x0 - x1) * fx2y2) / (x2 - x1) * (y2 - y1);
		// dgdx = (-(y2 - y0) * gx1y1 + (y2 - y0) * gx2y1 - (y0 - y1) * gx1y2 + (y0 - y1) * gx2y2) / (x2 - x1) * (y2 - y1);
		// dgdy = (-(x2 - x0) * gx1y1 - (x0 - x1) * gx2y1 + (x2 - x0) * gx1y2 + (x0 - x1) * gx2y2) / (x2 - x1) * (y2 - y1);
		// Use the icMatrix2x2 class to store these values in a 2x2 matrix
		

		// x/y0 are the singularity points
		// x/y1 are the mins of the verts
		// x/y2 are the maxs of the verts
		double x1 = smallestX;
		double x2 = largestX;
		double y1 = smallestY;
		double y2 = largestY;

		double dfdx = (-(y2 - y) * fx1y1 + (y2 - y) * fx2y1 - (y - y1) * fx1y2 + (y - y1) * fx2y2) / (x2 - x1) * (y2 - y1);
		double dfdy = (-(x2 - x) * fx1y1 - (x - x1) * fx2y1 + (x2 - x) * fx1y2 + (x - x1) * fx2y2) / (x2 - x1) * (y2 - y1);
		double dgdx = (-(y2 - y) * gx1y1 + (y2 - y) * gx2y1 - (y - y1) * gx1y2 + (y - y1) * gx2y2) / (x2 - x1) * (y2 - y1);
		double dgdy = (-(x2 - x) * gx1y1 - (x - x1) * gx2y1 + (x2 - x) * gx1y2 + (x - x1) * gx2y2) / (x2 - x1) * (y2 - y1);


		icMatrix2x2 jacobian;
		jacobian.entry[0][0] = (dfdx);
		jacobian.entry[0][1] = (dfdy);
		jacobian.entry[1][0] = (dgdx);
		jacobian.entry[1][1] = (dgdy);

		// 5. Use the determinant() function in the icMatrix.H file to compute the determinant of the jacobian.
		// Or you can compute the determinant manually.
		
		double determ = determinant(jacobian);
		
		// 5.5 (553 students only) Find the eigenvalues of the Jacobian matrix 
		// by solving the quadratic equation P(λ)=det(Jacobian - λI)
		// where λ is the variable being solved for and I is the identity matrix.
		// 5. Use the Jacobian determinant to classify the singularity:
		// If the determinant is greater than 0, the singularity is a source or a sink.
		// (553 students only) Use the eigenvalues to differentiate between sources and sinks.
		// sources.push_back(icVector3(x,y,0));
		// If the determinant is exactly equal to zero, the singularity is higher order.
		// higher_order.push_back(icVector3(x,y,0));
		// If the determinant is less than zero, the singularity is a saddle.
		// saddles.push_back(icVector3(x,y,0));
		if (determ > 0) {
			//We have a source/sink
			//I think I have to find if it is a source or sink but idk
			//Right now just call them all sources
			sources.push_back(icVector3(x, y, 0));
		}
		if (determ < 0) {
			//We have a saddel
			saddles.push_back(icVector3(x, y, 0));
		}
		if (determ == 0) {
			//We have a higher order(center/focus?)
			higher_order.push_back(icVector3(x, y, 0));
		}

	}
}



// find the distance from the input position ((x,y) coordinates)
// to the nearest singularity
double sing_prox(icVector2 pos)
{
	double prox = DBL_MAX;
	for (int i = 0; i < sources.size(); i++)
	{
		// get the (x,y) coordinates of the
		// singularity in an icVector2 object.
		icVector2 spos;
		spos.x = sources[i].x;
		spos.y = sources[i].y;
		double dist = length(pos - spos);
		if (dist < prox)
			prox = dist;
	}
	for (int i = 0; i < saddles.size(); i++)
	{
		// get the (x,y) coordinates of the
		// singularity in an icVector2 object.
		icVector2 spos;
		spos.x = saddles[i].x;
		spos.y = saddles[i].y;
		double dist = length(pos - spos);
		if (dist < prox)
			prox = dist;
	}
	for (int i = 0; i < higher_order.size(); i++)
	{
		// get the (x,y) coordinates of the
		// singularity in an icVector2 object.
		icVector2 spos;
		spos.x = higher_order[i].x;
		spos.y = higher_order[i].y;
		double dist = length(pos - spos);
		if (dist < prox)
			prox = dist;
	}
	return prox;
}


// takes x and y coordinates of a seed point and draws a streamline through that point
PolyLine build_streamline(double x, double y)
{
	// 1. Initialize current quad, current position, new position, step counter, and PolyLine variables
	// Quad* cquad = /*use the find_quad() method of the polyhedron class*/
	// icVector2 cpos = ???
	// icVector2 npos;
	// int step_counter = 0;
	// PolyLine pline;
	
	Quad* cquad = poly->find_quad(x, y);
	icVector2 cpos(x,y);
	icVector2 npos;
	int step_counter = 0;
	PolyLine pline;

	// 2. Trace the streamline forward until you hit a singularity, the edge of the domain, or the max step limit
	while (cquad != NULL && step_counter < STEP_MAX)
	{
		// take a step forward using streamline_step()
		cquad = streamline_step(cpos, npos, cquad, TRUE);
		// create a LineSegment object between cpos and npos
		LineSegment line(cpos.x, cpos.y, 0, npos.x, npos.y, 0);
		// push the LineSegment object onto the PolyLine pline
		pline.push_back(line);
		// set cpos to npos for next loop
		cpos = npos;
		// increment counter
		step_counter++;
	}
	// 3. Reset the cquad, cpos, and step_counter variables
	cquad = poly->find_quad(x, y);
	cpos.x = x;
	cpos.y = y;
	step_counter = 0;
	// 4. Trace the streamline backward until you hit a singularity, the edge of the domain, or the max step limit.
	// This should have the exact same structure as step 2 but with a slight change to the streamline_step function call
	while (cquad != NULL && step_counter < STEP_MAX)
	{
		// take a step forward using streamline_step()
		cquad = streamline_step(cpos, npos, cquad, FALSE);
		// create a LineSegment object between cpos and npos
		LineSegment line(cpos.x, cpos.y, 0, npos.x, npos.y, 0);
		// push the LineSegment object onto the PolyLine pline
		pline.push_back(line);
		// set cpos to npos for next loop
		cpos = npos;
		// increment counter
		step_counter++;
	}
	// 5. Return the PolyLine object
	return(pline);
}


Quad* streamline_step(icVector2& cpos, icVector2& npos, Quad* cquad, bool forward)
{
	double x1, y1, x2, y2, f11, f12, f21, f22, g11, g21, g12, g22;
	Vertex* v11, * v12, * v21, * v22;
	// 1. Find x1, x2, y1, y2, f11, f21, f12, f22, g11, g21, g12, g22, v11, v21, v12, v22
	// x1 is the smallest x coordinate of the 4 vertices
	// x2 is the largest x coordinate of the 4 vertices
	// y1 is the smallest y coordinate of the 4 vertices
	// y2 is the largest y coordinate of the 4 vertices
	// v11 is the vertex with coordinates (x1,y1)
	// f11 is the vector x component at vertex v11
	// g11 is the vector y component at vertex v11
	
	x1 = x2 = cquad->verts[0]->x;
	y1 = y2 = cquad->verts[0]->y;

	for (int j = 0; j < 4; j++) {
		if (x1 > cquad->verts[j]->x) {
			x1 = cquad->verts[j]->x;
		}
		if (x2 < cquad->verts[j]->x) {
			x2 = cquad->verts[j]->x;
		}
		if (y1 > cquad->verts[j]->y) {
			y1 = cquad->verts[j]->y;
		}
		if (y2 < cquad->verts[j]->y) {
			y2 = cquad->verts[j]->y;
		}
	}

	for (int j = 0; j < 4; j++) {
		if (cquad->verts[j]->x == x1 && cquad->verts[j]->y == y1) {
			v11 = cquad->verts[j];
			f11 = cquad->verts[j]->vx;
			g11 = cquad->verts[j]->vy;
		}
		if (cquad->verts[j]->x == x1 && cquad->verts[j]->y == y2) {
			v12 = cquad->verts[j];
			f12 = cquad->verts[j]->vx;
			g12 = cquad->verts[j]->vy;
		}
		if (cquad->verts[j]->x == x2 && cquad->verts[j]->y == y1) {
			v21 = cquad->verts[j];
			f21 = cquad->verts[j]->vx;
			g21 = cquad->verts[j]->vy;
		}
		if (cquad->verts[j]->x == x2 && cquad->verts[j]->y == y2) {
			v22 = cquad->verts[j];
			f22 = cquad->verts[j]->vx;
			g22 = cquad->verts[j]->vy;
		}
	}

	// 2. Get the coordinates (x0,y0) and find the vector components at cpos.
	// icVector2 vect;
	// vect.x = (use bilinear interpolation)
	// vect.y = (use bilinear interpolation)
	// normalize the vector
	
	//Trying to find the vector components at cpos, but only know the surrounding 4 vertexes vector components
	icVector2 vect;

	vect.x = (x2 - cpos.x) / (x2 - x1) * (y2 - cpos.y) / (y2 - y1) * (f11)
		+ (cpos.x - x1) / (x2 - x1) * (y2 - cpos.y) / (y2 - y1) * f21
		+ (x2 - cpos.x) / (x2 - x1) * (cpos.y - y1) / (y2 - y1) * f12
		+ (cpos.x - x1) / (x2 - x1) * (cpos.y - y1) / (y2 - y1) * f22;

	vect.y = (x2 - cpos.x) / (x2 - x1) * (y2 - cpos.y) / (y2 - y1) * g11
		+ (cpos.x - x1) / (x2 - x1) * (y2 - cpos.y) / (y2 - y1) * g21
		+ (x2 - cpos.x) / (x2 - x1) * (cpos.y - y1) / (y2 - y1) * g12
		+ (cpos.x - x1) / (x2 - x1) * (cpos.y - y1) / (y2 - y1) * g22;



	normalize(vect);

	// 3. Check the direction and calculate the new position
	// if (!forward) {vect *= -1.0;}
	// calculate new position npos using cpos, STEP, and vect
	
	if (!forward) {
		vect *= -1.0;
	}
	npos = cpos + STEP * vect;

	// 4. If npos is outside the current quad cquad, then we need to find the crossing
	// point where the streamline leaves the cquad, and determine what quad we will use
	// for the next step.
	Quad* nquad = cquad; //guess that the next quad will be the same
	if (npos.x < x1 || npos.x > x2 || npos.y < y1 || npos.y > y2)
	{
		// set up local variables
		Edge* cross_edge;
		// use parametric forms of lines to find the intersection along
		// each edge. (this method is described in the lecture lecture)
		
		double vx = npos.x - cpos.x;
		double vy = npos.y - cpos.y;
		
		double ex = 0, ey = 0;

		icVector2 cross;

		if (vx > 0) {
			ex = x2;
		}
		else {
			ex = x1;
		}
		if (vy > 0) {
			ey = y2;
		}
		else {
			ey = y1;
		}

		if (vx == 0) {
			cross.x = cpos.x;
			cross.y = ey;
		}
		else if (vy == 0) {
			cross.x = ex;
			cross.y = cpos.y;
		}
		else {
			double tx = (ex - cpos.x) / vx;
			double ty = (ey - cpos.y) / vy;
			if (tx <= ty) {
				cross.x = ex;
				cross.y = cpos.y + tx * vy;
			}
			else {
				cross.x = cpos.x + ty * vx;
				cross.y = ey;
			}
		}



		// once you have found each crossing point, check whether the vector at cpos
		// is pointing toward each crossing point.
		double dprod = dot(vect, cross - cpos);

		// Find the crossing point that is both on the border of cquad and
		// has a dprod value greater than 0.
		if (onEdge(cross, cquad) && dprod > 0)
		{
			// set npos to the crossing point
			npos = cross;
			Vertex* vert1 = NULL;
			Vertex* vert2 = NULL;
			// get cross_edge using the poly->find_edge() function
			for (int j = 0; j < 4; j++) {
				if (vert1 == NULL){
					if (cross.x == cquad->verts[j]->x || cross.y == cquad->verts[j]->y) {
						vert1 = cquad->verts[j];
					}
				}
				else if (vert2 == NULL) {
					if (cross.x == cquad->verts[j]->x || cross.y == cquad->verts[j]->y) {
						vert2 = cquad->verts[j];
					}
				}

			}
			cross_edge = poly->find_edge(vert1, vert2);
			// set nquad using the poly->other_quad() function
			nquad = poly->other_quad(cross_edge, cquad);
		}
		else {
			nquad = poly->find_quad(npos.x, npos.y);
		}
		// if none of the crossing points meet these conditions, use the poly->find_quad()
		// function with npos to get the appropriate nquad
	}
	// 5. Check the current singularity proximity to see if we should stop tracing the streamline.
	// You can do this using the sing_prox() function.
	// If the proximity is less than the step size, we don’t want to take another step

	if (sing_prox(npos) < STEP) {
		nquad = NULL;
	}

	// 6. Return nquad to be used in the next call to streamline_step().
	return nquad;
}


bool onEdge(icVector2 point, Quad* quad) {
	for (int i = 0; i < 4; i++) {
		if (point.x == quad->verts[i]->x || point.y == quad->verts[i]->y) {
			return(TRUE);
		}
	}
	return(FALSE);
}



/******************************************************************************
Callback function for dragging mouse
******************************************************************************/

void motion(int x, int y) {
	float r[4];
	float s, t;

	s = (2.0 * x - win_width) / win_width;
	t = (2.0 * (win_height - y) - win_height) / win_height;

	if ((s == s_old) && (t == t_old))
		return;

	switch (mouse_mode) {
	case 2:

		Quaternion rvec;

		mat_to_quat(rotmat, rvec);
		trackball(r, s_old, t_old, s, t);
		add_quats(r, rvec, rvec);
		quat_to_mat(rvec, rotmat);

		s_old = s;
		t_old = t;

		display();
		break;

	case 1:

		translation[0] += (s - s_old);
		translation[1] += (t - t_old);

		s_old = s;
		t_old = t;

		display();
		break;
	}
}

/******************************************************************************
Callback function for mouse clicks
******************************************************************************/

void mouse(int button, int state, int x, int y) {

	int key = glutGetModifiers();

	if (button == GLUT_LEFT_BUTTON || button == GLUT_RIGHT_BUTTON) {
		
		if (state == GLUT_DOWN) {
			float xsize = (float)win_width;
			float ysize = (float)win_height;

			float s = (2.0 * x - win_width) / win_width;
			float t = (2.0 * (win_height - y) - win_height) / win_height;

			s_old = s;
			t_old = t;

			/*translate*/
			if (button == GLUT_LEFT_BUTTON)
			{
				mouse_mode = 1;
			}

			/*rotate*/
			if (button == GLUT_RIGHT_BUTTON)
			{
				mouse_mode = 2;
			}
		}
		else if (state == GLUT_UP) {

			if (button == GLUT_LEFT_BUTTON && key == GLUT_ACTIVE_SHIFT) {  // build up the selection feedback mode

				/*select face*/

				GLuint selectBuf[512];
				GLint hits;
				GLint viewport[4];

				glGetIntegerv(GL_VIEWPORT, viewport);

				glSelectBuffer(win_width, selectBuf);
				(void)glRenderMode(GL_SELECT);

				glInitNames();
				glPushName(0);

				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();

				/*create 5x5 pixel picking region near cursor location */
				gluPickMatrix((GLdouble)x, (GLdouble)(viewport[3] - y), 1.0, 1.0, viewport);

				set_view(GL_SELECT);
				set_scene(GL_SELECT, poly);
				display_quads(GL_SELECT, poly);

				glMatrixMode(GL_PROJECTION);
				glPopMatrix();
				glFlush();

				glMatrixMode(GL_MODELVIEW);

				hits = glRenderMode(GL_RENDER);
				poly->selected_quad = processHits(hits, selectBuf);
				printf("Selected quad id = %d\n", poly->selected_quad);
				glutPostRedisplay();

				CHECK_GL_ERROR();

			}
			else if (button == GLUT_LEFT_BUTTON && key == GLUT_ACTIVE_CTRL)
			{
				/*select vertex*/

				GLuint selectBuf[512];
				GLint hits;
				GLint viewport[4];

				glGetIntegerv(GL_VIEWPORT, viewport);

				glSelectBuffer(win_width, selectBuf);
				(void)glRenderMode(GL_SELECT);

				glInitNames();
				glPushName(0);

				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();

				/*  create 5x5 pixel picking region near cursor location */
				gluPickMatrix((GLdouble)x, (GLdouble)(viewport[3] - y), 1.0, 1.0, viewport);

				set_view(GL_SELECT);
				set_scene(GL_SELECT, poly);
				display_vertices(GL_SELECT, poly);

				glMatrixMode(GL_PROJECTION);
				glPopMatrix();
				glFlush();

				glMatrixMode(GL_MODELVIEW);

				hits = glRenderMode(GL_RENDER);
				poly->selected_vertex = processHits(hits, selectBuf);
				printf("Selected vert id = %d\n", poly->selected_vertex);
				for (int j = 0; j < poly->nverts; j++) {
					if (j == poly->selected_vertex) {
						printf("Selected vert x: %f, y: %f\n", poly->vlist[j]->x, poly->vlist[j]->y);
					}
				}
				glutPostRedisplay();

			}

			mouse_mode = -1;
		}
	}
}

/******************************************************************************
Callback function for mouse wheel scroll
******************************************************************************/

void mousewheel(int wheel, int direction, int x, int y) {
	if (direction == 1) {
		zoom *= zoomspeed;
		glutPostRedisplay();
	}
	else if (direction == -1) {
		zoom /= zoomspeed;
		glutPostRedisplay();
	}
}

/******************************************************************************
Callback function for window reshaping
******************************************************************************/

void reshape(int width, int height)
{
	win_width = width;
	win_height = height;

	aspectRatio = (float)width / (float)height;

	glViewport(0, 0, width, height);

	set_view(GL_RENDER);

	// reset IBFV pixels buffer
	free(pixels);
	initIBFV();
}

/******************************************************************************
Display IBFV vector field visualization (used for Project 3)
******************************************************************************/

void displayIBFV()
{
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_BLEND);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

	glEnable(GL_TEXTURE_2D);
	glShadeModel(GL_FLAT);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glClearColor(0.5, 0.5, 0.5, 1.0);  // background for rendering color coding and lighting
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// draw the mesh using pixels and use vector field to advect texture coordinates
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, win_width, win_height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels);

	double modelview_matrix[16], projection_matrix[16];
	int viewport[4];
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix);
	glGetDoublev(GL_PROJECTION_MATRIX, projection_matrix);
	glGetIntegerv(GL_VIEWPORT, viewport);

	for (int i = 0; i < poly->nquads; i++)
	{
		Quad* qtemp = poly->qlist[i];

		glBegin(GL_QUADS);
		for (int j = 0; j < 4; j++)
		{
			Vertex* vtemp = qtemp->verts[j];

			double tx, ty, dummy;
			gluProject((GLdouble)vtemp->x, (GLdouble)vtemp->y, (GLdouble)vtemp->z,
				modelview_matrix, projection_matrix, viewport, &tx, &ty, &dummy);

			tx = tx / win_width;
			ty = ty / win_height;
			
			icVector2 dp = icVector2(vtemp->vx, vtemp->vy);
			normalize(dp);
			dp *= dmax;

			double dx = -dp.x;
			double dy = -dp.y;

			float px = tx + dx;
			float py = ty + dy;

			glTexCoord2f(px, py);
			glVertex3d(vtemp->x, vtemp->y, vtemp->z);
		}
		glEnd();
	}

	glEnable(GL_BLEND);

	// blend in noise pattern
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glTranslatef(-1.0, -1.0, 0.0);
	glScalef(2.0, 2.0, 1.0);

	glCallList(1);

	glBegin(GL_QUAD_STRIP);

	glTexCoord2f(0.0, 0.0);  glVertex2f(0.0, 0.0);
	glTexCoord2f(0.0, tmax); glVertex2f(0.0, 1.0);
	glTexCoord2f(tmax, 0.0);  glVertex2f(1.0, 0.0);
	glTexCoord2f(tmax, tmax); glVertex2f(1.0, 1.0);
	glEnd();
	glDisable(GL_BLEND);

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glReadPixels(0, 0, win_width, win_height, GL_RGB, GL_UNSIGNED_BYTE, pixels);

	// draw the mesh using pixels without advecting texture coords
	glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, win_width, win_height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels);
	for (int i = 0; i < poly->nquads; i++)
	{
		Quad* qtemp = poly->qlist[i];
		glBegin(GL_QUADS);
		for (int j = 0; j < 4; j++)
		{
			Vertex* vtemp = qtemp->verts[j];
			double tx, ty, dummy;
			gluProject((GLdouble)vtemp->x, (GLdouble)vtemp->y, (GLdouble)vtemp->z,
				modelview_matrix, projection_matrix, viewport, &tx, &ty, &dummy);
			tx = tx / win_width;
			ty = ty / win_height;
			glTexCoord2f(tx, ty);
			glVertex3d(vtemp->x, vtemp->y, vtemp->z);
		}
		glEnd();
	}
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);
	glShadeModel(GL_SMOOTH);
}

/******************************************************************************
Callback function for scene display
******************************************************************************/

void display(void)
{
	glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	set_view(GL_RENDER);
	set_scene(GL_RENDER, poly);

	/*display the mesh*/
	display_polyhedron(poly);

	/*display selected elements*/
	display_selected_vertex(poly);
	display_selected_quad(poly);


	glFlush();
	glutSwapBuffers();
	glFinish();

	CHECK_GL_ERROR();
}

/******************************************************************************
Diaplay the polygon with visualization results
******************************************************************************/

void display_polyhedron(Polyhedron* poly)
{
	unsigned int i, j;

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);

	switch (display_mode)
	{
	case 1:	// solid color display with lighting
	{
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHT1);

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		GLfloat mat_diffuse[4] = { 0.24, 0.4, 0.47, 0.0 };
		GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialf(GL_FRONT, GL_SHININESS, 50.0);

		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}
	}
	break;

	case 2:	// wireframe display
	{
		glDisable(GL_LIGHTING);
		glEnable(GL_LINE_SMOOTH);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glLineWidth(1.0);

		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];

			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_q->normal.entry[0], temp_q->normal.entry[1], temp_q->normal.entry[2]);
				glColor3f(0.0, 0.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}

		glDisable(GL_BLEND);
	}
	break;

	case 3:	// checkerboard pattern display
	{
		glDisable(GL_LIGHTING);
		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glColor3f(temp_v->R, temp_v->G, temp_v->B);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}
	}
	break;

	case 4: // points and lines drawing example
	{
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHT1);

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		GLfloat mat_diffuse[4] = { 0.24, 0.4, 0.47, 0.0 };
		GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialf(GL_FRONT, GL_SHININESS, 50.0);

		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}

		// draw lines
		for (int k = 0; k < lines.size(); k++)
		{
			drawPolyLine(lines[k], 1.0, 1.0, 0.0, 0.0);
		}

		// draw points
		for (int k = 0; k < points.size(); k++)
		{
			icVector3 point = points[k];
			drawDot(point.x, point.y, point.z);
		}

		break;
	}
	break;

	case 5:	// IBFV vector field display
	{
		displayIBFV();
		glutPostRedisplay();
	}
	break;

	case 6: // add your own display mode
	{	
		displayIBFV();
		glutPostRedisplay();

		for (int j = 0; j < poly->nverts; j++) {
			drawDot(poly->vlist[j]->x, poly->vlist[j]->y, poly->vlist[j]->z, .05f, 1, 0 ,0);
		}


	}
	break;

	case 7:	// Streamlines
	{
		displayIBFV();
		glutPostRedisplay();

		for (int k = 0; k < streamlines.size(); k++)
		{
			drawPolyLine(streamlines[k], 1.0, 1.0, 0.0, 0.0);
		}
	}
	break;


	case 9:	// With vertexes colored
	{
		display3DIBFV();
		for (int j = 0; j < poly->nverts; j++) {
			drawDot(poly->vlist[j]->x, poly->vlist[j]->y, poly->vlist[j]->z, .02f, 1, 0, 0);
		}
	}
	break;

	case 0:	// Without
	{
		display3DIBFV();
		//for (int j = 0; j < poly->nverts; j++) {
			//drawDot(poly->vlist[j]->x, poly->vlist[j]->y, poly->vlist[j]->z, .05f, 1, 0, 0);
		//}
	}
	break;

	default:
	{
		// don't draw anything
	}

	}
}
