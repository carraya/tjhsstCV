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

vector<Points> lloydAlgo(Points pts, int k) {
    vector<Points> finalCluster;
    double minVariance = 1e9;
    for(int x=0; x<50; x++) {    
        Points means;
        vector<Points> cluster;
        srand(time(0));
        for(int i=0; i<k; i++) {
            means.push_back(pts[rand()%pts.size()]);
        }
        
        for(int j=0; j<pts.size(); j++) {
            double meanClosest = 1e9;
            for(int i = 0; i<means.size(); i++) {
                double d = sqrt(pow(pts[j].first-means[i].first, 2.0)+pow(pts[j].second-means[i].second, 2.0));
                if(d < meanClosest) {
                    meanClosest = d;
                }
            }
            cluster[meanClosest].push_back(pts[j]);
        }
        double variance = 0;
        for(int i=0; i<cluster.size(); i++) {
            vector<double> distances;
            for(int j=0; j<cluster[i].size(); j++) {
                double d = sqrt(pow(cluster[i][j].first-means[i].first, 2.0)+pow(cluster[i][j].second-means[i].second, 2.0));
                distances.push_back(d);
            }
            double sumOfD = 0;
            for(int j=0; j<distances.size(); j++) {
                sumOfD += distances[j];
            }
            double mean = sumOfD/distances.size();
            double sumSqDif = 0;
            for(int j=0; j<distances.size(); j++) {
                double sqDif = pow(distances[j] - mean, 2.0);
                sumSqDif += sqDif;
            }
            variance += sumSqDif/(distances.size()-1);
        }
        variance /= cluster.size();
        if(variance < minVariance) {
            minVariance = variance;
            finalCluster.clear();
            finalCluster = cluster;
        }
    }
    return finalCluster;
}

int main(){
    srand(time(0));	
    ofstream myfile;
    myfile.open("writeHere.ppm");
    myfile << "P3 800 800 255\n";
    int numP=0;
    Points a;
    for(int n=0;n<100;n++) {
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
    vector<Points> clusters = lloydAlgo(a, 6);
    for(int j=0; j<clusters[0].size(); j++) {
        int x = (int) (clusters[0][j].first * 800);
        int y = (int) (clusters[0][j].second * 800);
        ls[y][x] = "0 0 0 ";
    }
    for(int j=0; j<clusters[1].size(); j++) {
        int x = (int) (clusters[1][j].first * 800);
        int y = (int) (clusters[1][j].second * 800);
        ls[y][x] = "255 0 0 ";
    }
    for(int j=0; j<clusters[2].size(); j++) {
        int x = (int) (clusters[2][j].first * 800);
        int y = (int) (clusters[2][j].second * 800);
        ls[y][x] = "0 255 0 ";
    }
    for(int j=0; j<clusters[3].size(); j++) {
        int x = (int) (clusters[3][j].first * 800);
        int y = (int) (clusters[3][j].second * 800);
        ls[y][x] = "0 0 255 ";
    }
    for(int j=0; j<clusters[4].size(); j++) {
        int x = (int) (clusters[4][j].first * 800);
        int y = (int) (clusters[4][j].second * 800);
        ls[y][x] = "255 51 51 ";
    }
    for(int j=0; j<clusters[5].size(); j++) {
        int x = (int) (clusters[5][j].first * 800);
        int y = (int) (clusters[5][j].second * 800);
        ls[y][x] = "153 153 255 ";
    }
    for (int row=0; row<800; row++){
        for(int x=0; x<799; x++){
			myfile << ls[row][x] <<  " ";
		}
		myfile << ls[row][799] << endl;
	}
	myfile.close();
    return 0;
}