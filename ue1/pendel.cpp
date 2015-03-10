#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>

using namespace std;

int main(int argn, char* argv[]){

ofstream f;
f.open("data.file");

double h=1/atof(argv[1]);  //1/Schrittweite
double T=atof(argv[2]);  //Zeitintervall (Obergrenze)
double x=M_PI/4; //Anfangswinkel
double y=0;  //Anfangsgeschwindigkeit
double t=0;
int i=0;

cout<<"Hallo"<<h<<T<<x<<y<<h+T+x+y<<endl;

for(t=0; t<=T; t=t+h){

	if((int)t>i){
		f<< t << " " << x << " " << y << " " << -cos(x)+y*y/2<<endl;//übergibt Zeit t, Winkel y und Geschwindigkeit x
		printf("\rBereits %f \%...",t/T);
 		i++;
	}
	x+=h*y;
	y+= -h*sin(x);
}
f.close();
return 0;
}
