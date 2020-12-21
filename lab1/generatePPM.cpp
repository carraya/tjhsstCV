#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <math.h>
using namespace std;
string ls[800][800];
void makeCircle(double centerX, double centerY, double radius)
{
	for (double theta=0; theta<=2*3.141592; theta+=0.001)
	{
		int x = (int) (centerX*800) + (radius*800) * cos(theta);
		int y = (int) (centerY*800) + (radius*800) * sin(theta);

		if (y < 800 && y >= 0 && x < 800 && x >= 0)
		{
            ls[x][y] = "0 0 0 ";
		}
	}
}
void makeCircles(double x1,double x2,double x3,double y1,double y2,double y3)
{
	double a = sqrt(pow(x2-x1, 2)+pow(y2-y1, 2));
	double b = sqrt(pow(x3-x2, 2)+pow(y3-y2, 2));
	double c = sqrt(pow(x1-x3, 2)+pow(y1-y3, 2));

	double s = .5 * (c + a + b);
	double inscribedRadius = sqrt(((s-c) * (s-a) * (s - b))/s);
	double circumscribedRadius = (c * a * b) / (4 * inscribedRadius * s);
	
	double t = 2*(y1 * (x2-x3)+y2 * (x3-x1) + y3 * (x1 - x2));
	double hyp1 = pow(y1, 2) + pow(x1, 2);
	double hyp2 = pow(y2, 2) + pow(x2, 2);
	double hyp3 = pow(y3, 2) + pow(x3, 2);
	double circumscribedCenterX = (1/t)*(hyp1 * (x2-x3)+hyp2*(x3-x1)+(hyp3*(x1-x2)));
	double circumscribedCenterY = (1/t)*(hyp1 * (y3-y2)+hyp2*(y1-y3)+(hyp3*(y2-y1)));
	
	double inscribedCenterX = (b * y1 + c * y2 + a * y3) / (c+a+b);
	double inscribedCenterY = (b * x1 + c * x2 + a * x3) / (c+a+b);

    makeCircle(inscribedCenterX, inscribedCenterY, inscribedRadius);
	makeCircle(circumscribedCenterX, circumscribedCenterY, circumscribedRadius);
}
void bresenham(double x1, double y1, double x2, double y2)
{
	bool check = (abs(y2-y1) > abs(x2-x1));
	if(check)
	{
		double temp = x1;
		x1=y1;
		y1=temp;
		double temp2=x2;
		x2=y2;
		y2=temp2;
	}
	if(x1 > x2)
	{
		double temp = x1;
		x1=x2;
		x2=temp;
		double temp2=y1;
		y1=y2;
		y2=temp2;
	}
	double dx = x2 - x1;
	double dy = abs(y2 - y1);
    int dir;
    if(y1 < y2)
        dir=1;
    else
        dir=-1;
	int j  = (int)(y1*800);
    double e = dy - dx;
	for(int i=(int)(x1*800); i<(int)(x2*800); i++)
	{
		if(check)
			ls[j][i] = "0 0 0 ";
		else
			ls[i][j] = "0 0 0 ";
		if(e >= 0)
		{
			j += dir;
			e -= dx;
		}
		e += dy;
	}
}
void generate()
{
	srand(time(0));	

    double x1=((double) rand() / (RAND_MAX));
    double y1=((double) rand() / (RAND_MAX));
    double x2=((double) rand() / (RAND_MAX));
    double y2=((double) rand() / (RAND_MAX));
    double x3=((double) rand() / (RAND_MAX));
    double y3=((double) rand() / (RAND_MAX));
    ls[(int)x1*800][(int)y1*800] = "0 0 0 ";
    ls[(int)x2*800][(int)y2*800] = "0 0 0 ";
    ls[(int)x3*800][(int)y3*800] = "0 0 0 ";
	bresenham(y1, x1, y2, x2);
	bresenham(y2, x2, y3, x3);
	bresenham(y3, x3, y1, x1);
	makeCircles(x1,x2,x3,y1,y2,y3);
}
int main(){
    ofstream myfile;
    myfile.open("writeHere.ppm");
    myfile << "P3 800 800 1\n";
    for(int r=0;r<800;r++){
        for(int c=0;c<800;c++){
            ls[r][c]="1 1 1 ";
        }
    }
    generate();
    ofstream file;
    file.open("writeHere.ppm");
    file << "P3 800 800 1" << endl;

	for (int row=0; row<800; row++){
        for(int x=0; x<799; x++){
			file << ls[row][x] <<  " ";
		}
		file << ls[row][799] << endl;
	}
	file.close();
} 