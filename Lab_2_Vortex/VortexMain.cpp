#include<GL/freeglut.h>
#include<math.h>
#define pi 3.14159265359
#define e 2.71828182846
#define y_0 3
#define x_0 1
#define intensity 1
#define steptime 0.1

double vortexCoordU(double x_k, double y_k) {
	return (intensity / (2 * pi))*((y_0 - y_k) / pow((x_0 - x_k),2)+pow((y_0 - y_k),2));
}

double vortexCoordV(double x_k, double y_k) {
	return (intensity / (2 * pi))*((x_0 - x_k) / pow((x_0 - x_k), 2) + pow((y_0 - y_k), 2));
}

double nextX(double x_k, double y_k) {
	return x_k + steptime*vortexCoordU(x_k, y_k);
}

double nextY(double x_k, double y_k) {
	return y_k + steptime*vortexCoordV(x_k, y_k);

}


void display(void) {

	const int n = 100;
	double arrX[n];
	double arrY[n];
	arrX[0] = x_0;
	arrY[0] = y_0;
	/*glColor3f(0.0, 0.0, 0.0);
	glBegin(GL_LINES);

	glVertex2d(-10, -10); //õ
	glVertex2d(10, -10);

	glVertex2d(-10, 10); //ó
	glVertex2d(-10, -10);

	/*for (int i = 1; i <= n; i++) {
		arrX[i] = nextX(arrX[i - 1], arrY[i - 1]);
		arrY[i] = nextY(arrX[i - 1], arrY[i - 1]);

	}
	glEnd();*/
	glLineWidth(2);
	glColor3f(1, 0, 0); 
	glBegin(GL_LINE_STRIP);
	glVertex2d(arrX[0], arrY[0]);
	for (int i = 1; i < n; i++) {
		arrX[i] = nextX(arrX[i - 1], arrY[i - 1]);
		arrY[i] = nextY(arrX[i - 1], arrY[i - 1]);
		glVertex2d(arrX[i], arrY[i]);
	}
	glEnd();
	glFinish();

}

void init(void) {
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glLoadIdentity();
}



int main(int argc, char **argv) {

	glutInit(&argc, argv);
	glutInitWindowSize(700, 700);
//	glutInitWindowPosition(100, 100);
//	glutInitDisplayMode(GL_RGB | GL_DEPTH | GL_DOUBLE);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);

	glutCreateWindow("Euler's Method");
	init();

	glutDisplayFunc(display);
	glutMainLoop();
	return 0;
}