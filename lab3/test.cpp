#include<iostream>
#include<cstdio>
#include<vector>
#include<cmath>
#include<algorithm>
#include<fstream>
#include<ctime>
#include<utility>
#include <unordered_map>
#include <bits/stdc++.h> 
using namespace std;
typedef pair<double,double> Point;
typedef vector<Point> Points;
string ls[800][800];
Point solution[1];

void makeCircle(double centerX, double centerY, double radius)
{
	for (double theta=0; theta<=2*3.141592; theta+=0.001)
	{
		int x = (int) (centerX*750) + (radius*750) * cos(theta);
		int y = (int) (centerY*750) + (radius*750) * sin(theta);

		if (y < 800 && y >= 0 && x < 800 && x >= 0)
		{
            ls[x][y] = "0 0 0 ";
		}
	}
}

void showSolution(double centerX, double centerY, double radius)
{
	for (double theta=0; theta<=2*3.141592; theta+=0.001)
	{
		int x = (int) (centerX*800) + (radius*800) * cos(theta);
		int y = (int) (centerY*800) + (radius*800) * sin(theta);

		if (y < 800 && y >= 0 && x < 800 && x >= 0)
		{
            ls[x][y] = "255 0 0 ";
		}
	}
}

double bruteForce(Points a) {
    int n = a.size();
    double minDis = 1e9;
    for(int i=0; i<n; i++) {
        for(int j=i+1; j<n; j++) {
            Point p1 = a[i]; Point p2 = a[j];
            double dis = sqrt(pow(p1.first-p2.first, 2.0) + pow(p1.second-p2.second, 2.0));
            if(dis < minDis) {
                minDis = dis;
            }
        }
    }
    return minDis;
}

double randomizedAlgo(Points myPoints) {
    double s = sqrt(pow(myPoints[1].first-myPoints[0].first, 2.0) + pow(myPoints[1].second-myPoints[0].second, 2.0));
    unordered_map<int, Points> myHash;
    for(int i=0; i<myPoints.size(); i++) {
        int sstx = (int) myPoints[i].first/(s/2);
        int ssty = (int) myPoints[i].second/(s/2);
        int sst = sstx*ssty;
        
        int N = (int) 1/(2*s);
        int numOfSS = (int) pow((float) N, 2.0);
        int up = 0;
        int down = 0;
        int left = 0;
        int right = 0;
        if(numOfSS%sst<=N-2) {
            right = 2;
        }
        else if(numOfSS%sst<=N-1) {
            right = 1;
        }
        else {
            right = 0;
        }
        if(numOfSS%sst>=2) {
            left = 2;
        }
        else if(numOfSS%sst>=1) {
            left = 1;
        }
        else {
            left = 0;
        }
        if((numOfSS%sst)+(2*N)<=numOfSS) {
            down = 2;
        }
        else if((numOfSS%sst)+(2*N)<=numOfSS-N) {
            down = 2;
        }
        else {
            down = 2;
        }
        if((numOfSS%sst)+(2*N)>=numOfSS) {
            up = 2;
        }
        if((numOfSS%sst)+(2*N)>=numOfSS+N) {
            up = 2;
        }
        else {
            up = 2;
        }
        double min = 1e9;
        // int rows = up+down;
        // int cols = right+left;
        vector<vector<Point>> searchThis;
        for(int j=0; j<=up; j++) {
            for(int z=0; z<=left; z++) {
                searchThis.push_back(myHash[sst-(N*j)-z]);
            }
            for(int z=0; z<right; z++) {
                searchThis.push_back(myHash[sst-(N*j)+z]);
            }
        }
        for(int j=0; j<down; j++) {
            for(int z=0; z<=left; z++) {
                searchThis.push_back(myHash[sst+(N*j)-z]);
            }
            for(int z=0; z<right; z++) {
                searchThis.push_back(myHash[sst+(N*j)+z]);
            }
        }
        double min = 1e9;
        for(int j=0; j<searchThis.size(); j++) {
            double temp = bruteForce(searchThis[j]);
            if(temp < s) {
                min = temp;
            }
        }
        if(min < s) {
            myHash.clear();
            for(int j=0; j<i; j++) {
                int sstxP = (int) myPoints[j].first/(min/2);
                int sstyP = (int) myPoints[j].second/(min/2);
                int sstP = sstxP*sstyP;
                myHash[sstP].push_back(myPoints[j]);
            }
        }
        else {
            myHash[sst].push_back(myPoints[i]);
        }

    }
    
}

int main(){
    srand(time(0));	
    ofstream myfile;
    myfile.open("writeHere.ppm");
    myfile << "P3 800 800 255\n";
    clock_t start, end;
    start = clock();
    int numP=0;
    Points a;
    for(int n=0;n<10000;n++) {
        double x=((double) rand() / (RAND_MAX));
        double y=((double) rand() / (RAND_MAX));
        Point b = make_pair(x,y);
        a.push_back(b);
        numP++;
    }
    for(int r=0;r<800;r++){
        for(int c=0;c<800;c++){
            ls[r][c]="255 255 255 ";
        }
    }   
    for(int n=0;n<29;n++) {
        Point p = a[n];
        makeCircle(p.first, p.second, 0.01);
    }
    showSolution(solution[0].first, solution[0].second, 0.01);
    showSolution(solution[1].first, solution[1].second, 0.01);
    for (int row=0; row<800; row++){
        for(int x=0; x<799; x++){
			myfile << ls[row][x] <<  " ";
		}
		myfile << ls[row][799] << endl;
	}
	myfile.close();
    printf("Shortest distance: %f\n",randomizedAlgo(a));
    end = clock();
    double time_taken = double(end-start) / double(CLOCKS_PER_SEC);
    cout << "\nTime taken: " << fixed << time_taken << setprecision(8);
    return 0;
}