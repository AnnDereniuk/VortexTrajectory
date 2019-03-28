#include<GL/freeglut.h>
#include<math.h>
#include<vector>

#define pi 3.14159265359
#define e 2.71828182846
#define y_0 0.2
#define x_0 1.2
#define intensity 1.
#define steptime 0.1


std::vector <double> dotsEulX = {};
std::vector <double> dotsEulY = {};

std::vector <double> dotsAdamsX = {};
std::vector <double> dotsAdamsY = {};


const int n = 5000;

float coef = intensity / (2. * pi);
double vortXarr[]{ -1., -1., 1., 1. };	//vortices coordinates
double vortYarr[]{ -1., 1., 1., -1. };	//vortices coordinates

//calculating u(x,y), v(x,y) -  (from 1 to 4) for each vortex 

double nextEulX(double x_k, double y_k)
{
	float u = 0.;
	for (int i = 0; i < 4; i++) {
		u = u + coef*(vortYarr[i] - y_k)/(pow((x_k - vortXarr[i]), 2) + pow((y_k - vortYarr[i]), 2));
	}
	u = x_k + steptime * u;
	return u;
}

double nextEulY(double x_k, double y_k)
{
	float v = 0.;
	for (int i = 0; i < 4; i++) {
		v = v + coef * (-vortXarr[i] + x_k)/(pow((x_k - vortXarr[i]), 2) + pow((y_k - vortYarr[i]), 2));
		//(((x_k - vortXarr[i])*(x_k - vortXarr[i])) + ((y_k - vortYarr[i])*(y_k - vortYarr[i]))));
	}
	v = y_k + steptime * v;
	return v;
}
//-------------------------------------

double nextAdamsX(double x_1, double y_1, double x_2, double y_2) {
	float u_1 = 0.;
	float u_2 = 0.;
	float res = 0.;

	for (int i = 0; i < 4; i++) {
		u_1 = u_1 + coef * (vortYarr[i] - y_1) / (pow((x_1 - vortXarr[i]), 2) + pow((y_1 - vortYarr[i]), 2));
		u_2 = u_2 + coef * (vortYarr[i] - y_2) / (pow((x_2 - vortXarr[i]), 2) + pow((y_2 - vortYarr[i]), 2));
	}
	res = x_2 + (3 * u_2 - u_1)*(steptime / 2);
	return res;
}

double nextAdamsY(double x_1, double y_1, double x_2, double y_2) {
	float v_1 = 0.;
	float v_2 = 0.;
	float res = 0.;
	for (int i = 0; i < 4; i++) {
		v_1 = v_1 + coef * (-vortXarr[i] + x_1) / (pow((x_1 - vortXarr[i]), 2) + pow((y_1 - vortYarr[i]), 2));
		v_2 = v_2 + coef * (-vortXarr[i] + x_2) / (pow((x_2 - vortXarr[i]), 2) + pow((y_2 - vortYarr[i]), 2));
	}
	res = y_2 + (3 * v_2 - v_1)*(steptime / 2);
	return res;
}

void Reshape(int width, int height)
{
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(-5, 5, -5, 5);
	glMatrixMode(GL_MODELVIEW);
}

void Draw(void)
{
	double nextEulDotX, nextEulDotY;
	double nextAdamsDotX, nextAdamsDotY;

	glClear(GL_COLOR_BUFFER_BIT);
	glColor3d(0.0, 0.0, 0.0);
	glLineWidth(1);
	//*****************************
	glBegin(GL_LINES);
	glVertex2f(0., -10.0f); //y
	glVertex2f(0., 10.0f);
	glVertex2f(-10.0f, 0.);//x
	glVertex2f(10.0f, 0.);
	glEnd();
	//*****************************
	glColor3d(0.0, 0.0, 0.0);
	glPointSize(3.f);
	glBegin(GL_POINTS);
	glVertex2f(-1.0f, -1.0f);
	glVertex2f(1.0f, 1.0f);
	glVertex2f(1.0f, -1.0f);
	glVertex2f(-1.0f, 1.0f);
	glEnd();
	//*****************************
	//glBegin(GL_LINE_STRIP);
	dotsEulX.push_back(x_0);
	dotsEulY.push_back(y_0);

	for (int j = 0; j < n; j++) {
		nextEulDotX = nextEulX(dotsEulX[dotsEulX.size()-1], dotsEulY[dotsEulY.size()-1]);
		nextEulDotY = nextEulY(dotsEulX[dotsEulX.size()-1], dotsEulY[dotsEulY.size()-1]);
		dotsEulX.push_back(nextEulDotX);
		dotsEulY.push_back(nextEulDotY);
	}


	glColor3d(0.0, 0.0, 1.0);
	glBegin(GL_POINTS);
	for (int i = 0; i < n-1;i++) {
		glVertex2d(dotsEulX[i], dotsEulY[i]);
	}
	glEnd();
	//glutSwapBuffers();
/*
	dotsAdamsX.push_back(dotsEulX[0]);
	dotsAdamsY.push_back(dotsEulY[0]);
	dotsAdamsX.push_back(dotsEulX[1]);
	dotsAdamsY.push_back(dotsEulY[1]);
	for (int j = 0; j < n; j++) {
		nextAdamsDotX = nextAdamsX(dotsAdamsX[dotsAdamsX.size()-1], dotsAdamsY[dotsAdamsY.size() - 1], dotsAdamsX[dotsAdamsX.size() - 2], dotsAdamsY[dotsAdamsY.size() - 2]);
		nextAdamsDotY = nextAdamsY(dotsAdamsX[dotsAdamsX.size() - 1], dotsAdamsY[dotsAdamsY.size() - 1], dotsAdamsX[dotsAdamsX.size() - 2], dotsAdamsY[dotsAdamsY.size() - 2]);
		dotsAdamsX.push_back(nextAdamsDotX);
		dotsAdamsY.push_back(nextAdamsDotY);
	}
	glColor3d(1.0, 0.0, 0.0);

	glBegin(GL_POINTS);
	for (int i = 0; i <(n/10.);i++) {
		glVertex2d(dotsAdamsX[i], dotsAdamsY[i]);
	}
	glEnd();*/

	glFlush();
}

int main(int argc, char *argv[])
{
	glutInit(&argc, argv);
	glutInitWindowSize(1000, 700);
	glutInitWindowPosition(100, 100);

	glutInitDisplayMode(GLUT_SINGLE|GLUT_RGB);
	glutCreateWindow("Trajectory");

	glutReshapeFunc(Reshape);
	glutDisplayFunc(Draw);
	glClearColor(1., 1., 1., 1.);

	glutMainLoop();
	return 0;
}