/*********

A Gillespie Stochastic Simulation of the model described in MJ Taylor et. al.; BioRxiv, 062877.
(core algorithm from Gillespie, Daniel T. (1977). "Exact Stochastic Simulation of Coupled Chemical Reactions". The Journal of Physical Chemistry. 81 (25): 2340–2361.)

-For each value of ligand density L and ligand affinity \tau_u, this code simulates 250 cells for 500s each.
This particular example generates a text file with, for each cell, the largest number of ligands at a cluster, the largest number of ZAP70 positive receptors at a cluster (both over the entire 500s span) and the total number of clusters.

Compiled with gcc using the -O3 flag for compiler optimisation.

********/


# include <stdlib.h>
# include <stdio.h>
# include <math.h>
// # include <time.h.>
# include <string.h>
# include <dirent.h>
# include <sys/stat.h>
// # include "pcg-c-basic-0.9/pcg_basic.h"

int main();
double *bindergillespie(double kb0, double kb_c, double ku, int numberofruns, int maximalcluster);
int denovoBinding(double denovoRate, double totalTime, double *denovoTimes);
int runAndWriteCluster(double StartTime, double totalTime, double kb, double ku, double kZAPon, double kZAPoff, int *lEvents, int *zEvents);
float ran2(long *idum);
double gaussrnd(double mu, double sigmaval);
double exprnd(double krate);
int log10Spacer(double startval, double endval, int numberofpoints, double *arrayVals);

// Random Number Seed
long idum = -1;

int main ( )
// Run The 'Gillespie'
{
	// Parameters
	int NumberOfCells = 250;
	double totalTime = 500.0; //total time to run the simulation, in seconds
	double deNovoRate = 1.0/6.0; // Extracted from 16mer
	double kZAPon = 1.0/51.0; //extracted from my stochastic model
	double kZAPoff = 1.0/15.0; //estimated by Marcus
	int numberofL = 40;
	int numberofku = 30;

	// Master Directory
	char supermaster[300];
	sprintf(supermaster, "ScanSmart/");
	mkdir(supermaster,0777);

	// Generate variables to scan
	int stupid;
	int totalNumberEvents = 1000000;
	// Generate logarithmically spaced kb values.
	double *ArrL = (double *)malloc((numberofL+2)*sizeof(double));
	stupid = log10Spacer(0.01, 600.0,numberofL,ArrL);
	// Generate logarithmically spaced ku values.
	double *ArrTau = (double *)malloc(numberofku*sizeof(double));
	stupid = log10Spacer(1, 1000.0,numberofku,ArrTau);

	ArrL[40] = 795.5;
	ArrL[41] = 1054.8;
	
	numberofL = 42;

	// Run the simulation
	for(int lnum = 0; lnum < numberofL; lnum ++)
	{
		double kb = (1.0/30.0)*ArrL[lnum];
		for(int kunum = 0; kunum < numberofku; kunum++)
		{
			double ku = 1.0/ArrTau[kunum];
			char masterdirname[300];
			//Save Parameters
			sprintf(masterdirname,"%s%s%0.2f%s%0.2f%s%0.5f%s%0.5f%s%0.5f%s%0.5f%s%0.5f%s%0.0f",supermaster,"tau_",ArrTau[kunum],"_L_",ArrL[lnum],"_Data_Rates_kb0_",deNovoRate,"_kb_", kb , "_ku_", ku, "_kZon_", kZAPon, "_kZAPoff_", kZAPoff, "_totalTime_", totalTime);
			// mkdir(masterdirname,0777);

			char ligfile[500];
			char zapfile[500];
			sprintf(ligfile,"%s%s",masterdirname,".txt");
			FILE *lfile;
			lfile = fopen(ligfile,"wt");

			int clusterNumAll[250];
			for(int cellno = 0; cellno < NumberOfCells; cellno++)
			{
				//generate denovo times
				double *denovoTimes = (double *)malloc(sizeof(double)*600);
				stupid = denovoBinding(deNovoRate, totalTime,denovoTimes);
				// Make variable
				int *lNum = (int *)malloc(sizeof(int));
				int *zNum = (int *)malloc(sizeof(int));
				lNum[0] = 0;
				zNum[0] = 0;
				int maxNCell = 0;
				int maxMCell = 0;
				int clusAtEnd = 0;
				//run cluster by cluster, get profiles
				int clustersLeft = 1;
				int clusterID = 0;
				do
				{
					double nextTime = denovoTimes[clusterID];
					if(nextTime > 0)
					{
						clusAtEnd = clusAtEnd + runAndWriteCluster(nextTime, totalTime, kb, ku, kZAPon, kZAPoff, lNum, zNum);
						if(lNum[0] > maxNCell)
						{
							maxNCell = lNum[0];
						}
						if(zNum[0] > maxMCell)
						{
							maxMCell = zNum[0];
						}
						clusterID++;
					}
					else
					{
						clustersLeft = 0;
					}
				}
				while(clustersLeft == 1);

				fprintf(lfile, "%d \t %d \t %d \n", maxNCell, maxMCell,clusAtEnd);

				free(denovoTimes);

				free(lNum);
				free(zNum);

				printf("%s%d%s%d%s%d%s%d%s%d%s%d\n", "Finished Cell ", cellno+1, " of ", NumberOfCells, " and ", kunum + 1, " of ", numberofku, " and ", lnum + 1, " of ", numberofL );
			}
			fclose(lfile);
		}
	}

	return 0;
}



/****************
	Gillespie
	********************/

int denovoBinding(double denovoRate, double totalTime, double *denovoTimes)
{
	int amIdoneYet = 0;
	double t = 0;
	int clusterssofar = 0;
	do
	{
		double nextTime = exprnd(denovoRate);
		t = t + nextTime;
		if(t <= totalTime)
		{
			denovoTimes[clusterssofar] = t;
			clusterssofar++;
		}
		else
		{
			denovoTimes[clusterssofar] = -1.0;
			amIdoneYet = 1;
		}
	}
	while(amIdoneYet == 0);

	return clusterssofar;
}

int runAndWriteCluster(double StartTime, double totalTime, double kb, double ku, double kZAPon, double kZAPoff, int *lEvents, int *zEvents) // int runAndWriteCluster(double StartTime, double totalTime, double kb, double ku, double kZAPon, double kZAPoff, int CellNo, int ClusterNo, char celldirname[], double *lEvents, double *zEvents)
{
	//First, the filename stuff
	int stupid;

	//Now, run until we're done. Internal variables are m = ligand number and mzap = zapated receptors.
	int m = 1;
	int mzap = 0;
	int m2;
	int mzap2;
	int arewedone = 0;
	int lCount = 0;
	int zCount = 0;
	int lHappen = 0;
	int zHappen = 0;
	double t = StartTime;
	int maxM = 10;

	// lEvents[lCount] = StartTime;
	lCount++;
	int maxLig = 0;
	int maxZap = 0;
	do
	{
		//what are the rates?
		double totalrate = kb + ku*((double) m) + kZAPoff*((double) mzap) + kZAPon*((double) (m-mzap));
		double totalwaitingTime = exprnd(totalrate);
		t = t+totalwaitingTime;
		if(t > totalTime)
		{
			arewedone = 1;
			lEvents[0] = maxLig;
			zEvents[0] = maxZap;
			return 1;
			// printf("%d %d \n", m, lEvents[0]);
		}
		else
		{
			//which reaction was it?
			double randoman = ran2(&idum)*totalrate;
			if(randoman <= kb)
			{
				m2 = m+1;
				mzap2 = mzap;
				lHappen = 1;
			}
			else if(randoman <= kb + ku*((double) m))
			{
				m2 = m-1;
				//does a fellow with ZAP unbind?
				double randoman2 = ran2(&idum)*((double) m );
				if(randoman2 <= ((double) mzap))
				{
					mzap2 = mzap - 1;
					lHappen = -1;
					zHappen = -1;
				}
				else
				{
					mzap2 = mzap;
					lHappen = -1;
				}
			}
			else if(randoman <= kb + ku*((double) m) + kZAPoff*((double) mzap))
			{
				mzap2 = mzap-1;
				m2 = m;
				zHappen = -1;
			}
			else
			{
				mzap2 = mzap+1;
				m2 = m;
				zHappen = 1;
			}
			m = m2;
			mzap = mzap2;
			lHappen = 0;
			zHappen = 0;

			if(m > maxLig)
			{
				maxLig = m;
			}
			if(mzap > maxZap)
			{
				maxZap = mzap;
			}
		}
		if(m < 1)
		{
			arewedone = 1;
			lEvents[0] = maxLig;
			zEvents[0] = maxZap;
		}
	}
	while(arewedone == 0);

	return 0;
}



/******************8
	Log Spacer
	****************8***/

int log10Spacer(double startval, double endval, int numberofpoints, double *arrayVals)
{
	double logstartval = log10(startval);
	double logendval = log10(endval);
	double spacer = (logendval - logstartval)/((double )(numberofpoints - 1));

	for(int i = 0; i < numberofpoints; i++)
	{
		arrayVals[i] = pow(10.0, logstartval + i*spacer);
	}

	return 0;
}




/********************
	Random Number Bits - ran2 taken from Numerical Recipes in C
				********************/

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
// Taken from Numerical Recipes in C
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;
	
	if (*idum <= 0) {
	if (-(*idum) < 1) *idum=1;
	else *idum = -(*idum);
	idum2=(*idum);
	for (j=NTAB+7;j>=0;j--) {
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	if (j < NTAB) iv[j] = *idum;
	}
	iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}


double exprnd (double krate)
// Generate an exponentially distributed deviate.
{
	double expnum = (1.0/krate)*log(1.0/ran2(&idum));
	return expnum;
}

int alreadyGeneratedNormal = 0;
double storedGaussianVal;
double gaussrnd (double mu, double sigmaval)
// polar method for generating normally distributed numbers
{
	if(alreadyGeneratedNormal == 1)
	{
		alreadyGeneratedNormal = 0;
		return (mu + sigmaval * storedGaussianVal);
	}
	else
	{
		double U1, U2, W, mult;
		double X1, X2;
		 
		do
		{
			double r1 = ran2(&idum);
			double r2 = ran2(&idum);
			U1 = -1 + ((double) r1 ) * 2;
			U2 = -1 + ((double) r2 ) * 2;
			W = pow (U1, 2) + pow (U2, 2);
	    }
	    while (W >= 1.0 || W == 0.0);
		 
		mult = sqrt ((-2 * log (W)) / W);
		X1 = U1 * mult;
		X2 = U2 * mult;

		storedGaussianVal = ((double) X2); 
		return (mu + sigmaval * (double) X1);
	}
}
