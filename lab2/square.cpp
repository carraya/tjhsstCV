#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <math.h>
#include <bits/stdc++.h>
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

void bresenham(double x1,double y1, double x2, double y2){
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

void generate(){
    srand(time(0));

    double x1=abs(((double) rand() / (RAND_MAX)));
    double y1=abs(((double) rand() / (RAND_MAX)));
    double x2=abs(((double) rand() / (RAND_MAX)));
    double y2=abs(((double) rand() / (RAND_MAX)));
    double x3=abs(((double) rand() / (RAND_MAX)));
    double y3=abs(((double) rand() / (RAND_MAX)));
    double x4=abs(((double) rand() / (RAND_MAX)));
    double y4=abs(((double) rand() / (RAND_MAX)));

    makeCircle(((x1+x2)/2), ((y1+y2)/2), sqrt(pow((x2-x1),2) + pow((y2-y1),2))/2);
    makeCircle(((x3+x4)/2), ((y3+y4)/2), sqrt(pow((x3-x4),2) + pow((y3-y4),2))/2);

    double centerx1 = (x2+x1)/2;
	double centery1 = (y2+y1)/2;
	double centerx2 = (x4+x3)/2;
	double centery2 = (y4+y3)/2;
	double midx1 = centerx1-((y2-y1)/2);
	double midy1 = centery1+((x2-x1)/2);
	double midx2 = centerx2-((y4-y3)/2);
	double midy2 = centery2+((x4-x3)/2);

	bresenham(x1,y1,x2,y2);
	bresenham(x3,y3,x4,y4);
	bresenham(centerx1,centery1,midx1,midy1);
	bresenham(centerx2,centery2,midx2,midy2);
	bresenham(midx1,midy1,midx2,midy2);

    
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