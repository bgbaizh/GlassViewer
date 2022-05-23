#include <iostream>
#include <exception>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <vector>
#include <complex>
using namespace std;
int main()
{
    double box[3][3]={1,1,5,1,5,1,5,1,1};
    double Height[3]={0};
    double iCrossj[3][3]={0};
    double iCrossjnorm[3]={0};
    double kdotiCrossj[3]={0};
    int index[3][3]={0,1,2,1,2,0,2,0,1};

    for(auto &m:index)
    {
        int k= m[0];
        int i =m[1];
        int j =m[2];
        for(auto &m2:index)
        {
            int k2= m2[0];
            int i2 =m2[1];
            int j2 =m2[2];
            iCrossj[k][k2]=box[i][i2]*box[j][j2]-box[i][j2]*box[j][i2];
        }
    }
    
    for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
            iCrossjnorm[i]+= iCrossj[i][j]*iCrossj[i][j];
            kdotiCrossj[i]+=box[i][j]*iCrossj[i][j];
            }
            iCrossjnorm[i]=pow(iCrossjnorm[i],0.5);
            Height[i]=abs(kdotiCrossj[i]/iCrossjnorm[i]);
        }


        
}