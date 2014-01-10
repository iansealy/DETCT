// steve qin.
// 07/01/08
//#include "stdafx.h"
#include<stdio.h>
#include<cstdlib>
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<string>
#include<math.h>
#include <time.h>
using namespace std;

#define MAX_LENGTH 500
#define PI 3.14159265
#define LIMIT 100
#define MIN(a,b) ((a) < (b) ? (a):(b))
#define MAX(a,b) ((a) > (b) ? (a):(b))

void readData(const char dataFileName[], vector<int> &order,
		vector<double> &data, int *nRow, double *totalcount);

double logIntPoisson(const int k, const double lambda);

double logPoisson(const double x, const double lambda);

double genPoisson(const double y, const double mu, const double alpha);

double logGenPoisson(const double y, const double mu, const double alpha);

double logIntGenPoisson(const int y, const double mu, const double alpha);

double logIntTrunc0GenPoisson(const int y, const double mu, const double alpha);

double logTrunc0GenPoisson(const double y, const double mu, const double alpha);

void pathFinder(const int count, const vector<int> &order,
		const vector<double> &data, const int chromosomeLengthInBins,
		int *&path, double *&proba, double *&logproba, double *&hits,
		const double totalMapReads, const double totalPeakReads,
		const double totalPeakArea, const double medianPeakBinCount,
		const int numPeaks, const int readcoverage, const int threshold);

void parameterEstimate(const double *peaksdata, const int peakscount,
		const int threshold, const double mu, const double alpha,
		double *mufore, double *alphafore);

double likelihood(const double *y, const int chromosomeLengthInBins,
		const int threshold, const double mu, const double alpha);

double unran(int *na, int *nb, int *nc);

double gammaln(double xx);

int main(int argc, char **argv) {
	int j, numCycles = 1;
	int count, chromosomeLengthInBins, binSize;
	int readcoverage = 0;//added 04/13/08
	vector<int> order;
	vector<double> data;
	int *path;
	double *hits;
	double *proba, *logproba;
	char inputname[MAX_LENGTH] = "10.select.chr21.txt";//"jy9.10.select.chr22.txt";//"sle.txt";
	char parasname[MAX_LENGTH] = "jy10.paras.txt";
	char outputname[MAX_LENGTH] = "out.txt";
	ofstream outPutFile;
	double totalMapReads = 0;
	double totalPeakReads = 0;
	double totalPeakArea = 0;
	double medianPeakWidth = 0;
	double medianPeakBinCount = 0;
	int numPeaks = 0;
	istringstream iss;
	string lineString;
	double temval;
	vector<double> datas;
	double totalcount;
	int threshold;

	if (argc != 4) {
		printf(
				"3 options need to be specified:\n\tinput file name,\n\tinformation file name,\n\toutputfile name.\n");
		exit(0);
	}
	for (j = 0; j < MAX_LENGTH; j++) {
		inputname[j] = argv[1][j];
		parasname[j] = argv[2][j];
		outputname[j] = argv[3][j];
	}

	ifstream inFile(parasname);
	if (!inFile) {
		cout << "Error opening input parameter file" << parasname << endl;
		exit(0);
	}
	for (j = 0; j < 9; j++) {
		getline(inFile, lineString);
		iss.clear();
		iss.str(lineString + " ");
		iss >> temval;
		datas.push_back(temval);
	}
	totalMapReads = datas[0];
	totalPeakReads = datas[1];
	totalPeakArea = datas[2];
	medianPeakWidth = datas[3];
	//	if(medianPeakWidth >800)
	//		medianPeakWidth = 800;
	numPeaks = (int) datas[4];
	readcoverage = (int) datas[5];
	threshold = (int) datas[6];
	binSize = (int) datas[7];
	chromosomeLengthInBins = (int) datas[8];

	medianPeakBinCount = medianPeakWidth / (double) binSize;
	path = new int[chromosomeLengthInBins];
	proba = new double[chromosomeLengthInBins];
	logproba = new double[chromosomeLengthInBins];
	srand((unsigned) time(NULL));
	readData(inputname, order, data, &count, &totalcount);
	for (j = 0; j < numCycles; j++) {
		pathFinder(count, order, data, chromosomeLengthInBins, path, proba,
				logproba, hits, totalMapReads, totalPeakReads, totalPeakArea,
				medianPeakBinCount, numPeaks, readcoverage, threshold);
	}
	outPutFile.open(outputname);
	if (!outPutFile) {
		cout << "ERROR: Unable to open file: " << outputname << endl;
		exit(30);
	}//end of if
	for (j = 0; j < chromosomeLengthInBins; j++) {
		if (proba[j] > 0.5)//0.01
		{
			outPutFile << j << "	" << path[j] << "	" << proba[j] << "	"
					<< hits[j] << "	" << logproba[j] << endl;
		}
	}
	delete[] path;
	delete[] proba;
	delete[] hits;
	outPutFile.close();
	return 0;
}//end of main

void pathFinder(const int count, const vector<int> &order,
		const vector<double> &data, const int chromosomeLengthInBins,
		int *&path, double *&proba, double *&logproba, double *&hits,
		const double totalMapReads, const double totalPeakReads,
		const double totalPeakArea, const double medianPeakBinCount,
		const int numPeaks, const int readcoverage, const int threshold) {
	int j;
	double (*logfnh)[2], trp[2][2];
	double p[2], p0, p1;
	double ratio, compa;
	double dif = 0, inside1, inside2, inside3, inside4;
	int na, nb, nc;
	//double lambdaback, lambdafore;
	double sum = 0, sum2 = 0, mean, var, nsize;
	double muback, alphaback;
	double mufore, alphafore;
	//double mufore,mualpha;
	double mu, alpha;
	int bgtotal, number;
	double *peaksdata;
	//-	int threshold;
	//ofstream outParameterFile;

	//outParameterFile.open("poisson.out");
	//if (!outParameterFile) {
	//	cout << "ERROR: Unable to open file poisson.out." << endl;
	//	exit(30);
	//}//end of if
	na = rand() + 1;
	nb = rand() - 1;
	nc = rand();
	// hits are for all 25bp window on the genome.
	// only some of them are non-zero.
	hits = new double[chromosomeLengthInBins + 1]; // to fix the glibc
	sum = 0;
	sum2 = 0;
	bgtotal = 0;
	//-	threshold = 6;
	for (j = 0; j < chromosomeLengthInBins; j++) {
		hits[j] = 0;
	}
	for (j = 0; j < count; j++) {
		if (order[j] > chromosomeLengthInBins) {
			cout << "read bins extend further than the chromosome size in bins " << order[j] << " "<< chromosomeLengthInBins << endl;
			exit(1);
		}
		hits[order[j]] = data[j];
		if (data[j] < threshold) {
			//			sum = sum + (double) floor(data[j]);
			//			sum2 = sum2 + (double) floor(data[j])*floor(data[j]);
			sum = sum + data[j];
			sum2 = sum2 + data[j] * data[j];
			bgtotal++;
		}
	}
	nsize = (double) chromosomeLengthInBins - count + bgtotal;
	mean = sum / nsize;
	var = (sum2 - (double) nsize * mean * mean) / (nsize - 1);
	muback = mean;
	alphaback = (sqrt(var / mean) - 1) / mean;
	cout << "background: mu = " << muback << " alpha = " << alphaback << endl;
	//outParameterFile << muback << "	" << alphaback << endl;
	//	double aa = logGenPoisson(5,5,2);
	peaksdata = new double[count - bgtotal];
	number = 0;
	sum = 0;
	sum2 = 0;
	for (j = 0; j < count; j++) {
		if (data[j] >= threshold) {
			peaksdata[number] = data[j];
			//			sum = sum + (double) floor(data[j]);
			//			sum2 = sum2 + (double) floor(data[j])*floor(data[j]);
			sum = sum + data[j];
			sum2 = sum2 + data[j] * data[j];
			number++;
		}
	}
	mean = sum / number;
	var = (sum2 - (double) number * mean * mean) / (number - 1);
	mu = mean;
	alpha = (sqrt(var / mean) - 1) / mean;
	cout << "foreground (raw): mu = " << mu << " alpha = " << alpha << endl;
	//outParameterFile << mu << "	" << alpha << endl;
	//double ff = likelihood(peaksdata,number,threshold,6,0.3);
	parameterEstimate(peaksdata, number, threshold, mu, alpha, &mufore,
			&alphafore);
	cout << "foreground: mu = " << mufore << " alpha = " << alphafore << endl;
	//outParameterFile << mufore << "	" << alphafore << endl;
	//outParameterFile.close();
	//	exit(0);

	for (j = 0; j < 10; j++)//200
	{
		double aa = logIntGenPoisson(j, muback, alphaback);
		double bb = logIntGenPoisson(j, mufore, alphafore);
		double cc = logIntTrunc0GenPoisson(j, mufore, alphafore);
		cout << "j= " << j << " " << exp(aa) << "	" << exp(bb) << "	"
				<< exp(cc) << endl;
	}//end of j
	//	exit(0);

	//+        lambdaback = readcoverage * totalMapReads /(3100*0.9);
	//        cout <<"lambda foreground = "<<lambdafore<<", lambda background = "<<lambdaback<<endl;
	//***** initial probabilities.
	p[1] = totalPeakArea / 1357;
	p[0] = 1 - p[1];
	//***** transition probabilities.
	trp[0][0] = p[0];
	trp[0][1] = p[1];
	//double a = exp(1.0/medianPeakBinCount * log(0.5));
	trp[1][1] = exp(1.0 / medianPeakBinCount * log(0.5));
	trp[1][0] = 1 - trp[1][1];
	//***** emission probabilities.
	logfnh = new double[chromosomeLengthInBins][2];
	logfnh[0][0] = 0;
	logfnh[0][1] = -LIMIT;
	/*	double aaa = logIntPoisson(1,1);
	 double bbb = logIntPoisson(2,1);
	 double ccc = logPoisson(1.8,1);
	 cout  <<"aaa= "<<aaa<<" bbb= "<<bbb<<" ccc= "<<ccc<<endl;
	 exit(0);
	 */
	for (j = 1; j < chromosomeLengthInBins; j++) {
		if (hits[j] == 0) {
			p0 = 1;
			p1 = 0;
		} else if (hits[j] > LIMIT) {
			p0 = 0;
			p1 = 1;
		} else {
			p0 = logGenPoisson(hits[j], muback, alphaback);
			//+++			p0 = logPoisson(hits[j],muback);
			//+ p0 = exp(logPoisson(hits[j],lambdaback));
			//p0 = logPoisson(8,lambda);
			//++                        p1 = logGenPoisson(hits[j],mufore,alphafore)
			p1 = logTrunc0GenPoisson(hits[j], mufore, alphafore);
			//+ p1 = exp(logPoisson(hits[j],lambdafore));
			//+++			if(hits[j]>2)
			//+++				cout <<"j= "<<j<<" hits= "<<hits[j]<<" p0= "<<p0<<" p1= "<<p1<<endl;
			p0 = 1 / (1 + exp(p1 - p0));
			p1 = 1 - p0;
		}
		//-		if((j>143889)&&(j<143895))
		//-			cout <<"hits = "<<hits[j]<<" "<<logfnh[j-1][0]<<" "<<logfnh[j-1][1]<<endl;
		inside3 = logfnh[j - 1][1] + log(trp[1][0]) - logfnh[j - 1][0] - log(
				trp[0][0]);
		inside4 = logfnh[j - 1][1] + log(trp[1][1]) - logfnh[j - 1][0] - log(
				trp[0][1]);
		if (p0 == 0) {
			logfnh[j][0] = -LIMIT;
			logfnh[j][1] = 0;
		} else if (p1 == 0) {
			logfnh[j][1] = -LIMIT;
			logfnh[j][0] = 0;
		} else {
			if (inside3 < 0) {
				logfnh[j][0] = log(p0) + logfnh[j - 1][0] + log(trp[0][0])
						+ log(1 + exp(inside3));
			} else
				logfnh[j][0] = log(p0) + logfnh[j - 1][1] + log(trp[1][0])
						+ log(1 + exp(-inside3));
			if (inside4 < 0)
				logfnh[j][1] = log(p1) + logfnh[j - 1][0] + log(trp[0][1])
						+ log(1 + exp(inside4));
			else
				logfnh[j][1] = log(p1) + logfnh[j - 1][1] + log(trp[1][1])
						+ log(1 + exp(-inside4));
			//			cout <<logfnh[j][0]<<" "<<logfnh[j][1]<<endl;
		}
	}
	//+++	exit(0);
	ratio = 1.0 / (1.0 + exp(
			logfnh[chromosomeLengthInBins - 1][1]
					- logfnh[chromosomeLengthInBins - 1][0]));
	proba[chromosomeLengthInBins - 1] = 1 - ratio;
	compa = unran(&na, &nb, &nc);
	if (compa <= ratio)
		path[chromosomeLengthInBins - 1] = 0;
	else
		path[chromosomeLengthInBins - 1] = 1;
	logproba[chromosomeLengthInBins - 1] = log(ratio);
	for (j = chromosomeLengthInBins - 2; j >= 0; j--) {
		ratio = 1 / (1 + (trp[1][path[j + 1]] / trp[0][path[j + 1]]) * exp(
				logfnh[j][1] - logfnh[j][0]));
		logproba[j] = log(ratio);
		/*
		 if((j>=143890)&&(j<143895))
		 {
		 cout << "path= "<<path[j+1]<<endl;
		 double aa = trp[1][path[j+1]];
		 double bb = trp[0][path[j+1]];
		 double cc = logfnh[j][1];
		 double dd = logfnh[j][0];
		 cout << "j= "<<j<<" "<<"aa= "<<aa<<" bb= "<<bb<<" cc= "<<cc<<" dd= "<<dd<<" ratio= "<<ratio<<endl;
		 }
		 */
		proba[j] = 1 - ratio;
		compa = unran(&na, &nb, &nc);
		if (compa <= ratio)
			path[j] = 0;
		else
			path[j] = 1;
	}
}//end of pathFinder

double logIntPoisson(const int k, const double lambda) {
	double x;

	x = k * log(lambda) - lambda - gammaln((double) k + 1);
	return (x);
}//end of logIntPoisson

double logPoisson(const double x, const double lambda) {
	double y, ax, apart, tem01, tem02;

	ax = floor(x);
	apart = x - ax;
	tem01 = logIntPoisson((int) ax, lambda);
	tem02 = exp(logIntPoisson((int) (ax + 1), lambda) - tem01);
	y = tem01 + log((1 - apart) + apart * tem02);
	return y;
}//end of logPoisson

void parameterEstimate(const double *peaksdata, const int peakscount,
		const int threshold, const double mu, const double alpha,
		double *mufore, double *alphafore) {// use Metropolis to get MLE of truncated generalized Possion distribution parameters.
	int j;
	int nn;
	double *aa, *bb, *out;
	double stepsize, theta;
	int na, nb, nc;
	double rando, test, ratio;
	double aatry, bbtry;
	double aamax, bbmax, outmax;
	double oldlike, newlike;

	na = rand() + 1;
	nb = rand() - 1;
	nc = rand();
	nn = 200;
	aa = new double[nn];
	bb = new double[nn];
	out = new double[nn];
	stepsize = 0.5;
	aa[0] = mu;
	bb[0] = alpha;
	out[0] = likelihood(peaksdata, peakscount, threshold, mu, alpha);
	oldlike = out[0];
	aamax = mu;
	bbmax = alpha;
	outmax = out[0];
	for (j = 1; j < nn; j++) {
		rando = unran(&na, &nb, &nc);
		theta = rando * 2 * PI;
		aatry = aa[j - 1] + stepsize * cos(theta);
		bbtry = bb[j - 1] + 0.1 * stepsize * sin(theta);
		//double tempa = likelihood(peaksdata,peakscount, threshold, aatry,bbtry);
		//double tempb = likelihood(peaksdata,peakscount, threshold, aa[j-1],bb[j-1]);
		newlike = likelihood(peaksdata, peakscount, threshold, aatry, bbtry);
		ratio = exp(newlike - oldlike);
		if (ratio > 1) {
			aa[j] = aatry;
			bb[j] = bbtry;
			out[j] = newlike;
			outmax = newlike;
			aamax = aatry;
			bbmax = bbtry;
			oldlike = newlike;
		} else {
			test = unran(&na, &nb, &nc);
			if (test < ratio) {
				aa[j] = aatry;
				bb[j] = bbtry;
				out[j] = newlike;
				oldlike = newlike;
			} else {
				aa[j] = aa[j - 1];
				bb[j] = bb[j - 1];
				out[j] = out[j - 1];
			}
		}//end of else
		//cout <<"j= "<<j<<" aa= "<<aa[j]<<" bb= "<<bb[j]<<" out= "<<out[j]<<endl;
	}//end of j
	//exit(0);
	*mufore = aamax;
	*alphafore = bbmax;
}//end of parameterEstimate

double likelihood(const double *y, const int chromosomeLengthInBins,
		const int threshold, const double mu, const double alpha) {
	int j;
	double con, sum;

	con = 0;
	for (j = 0; j < threshold; j++) {
		con = con + exp(logIntGenPoisson(j, mu, alpha));
	}
	sum = 0;
	//	double a = y[1862];
	//	a = logGenPoisson(y[1862],mu,alpha);

	for (j = 0; j < chromosomeLengthInBins; j++) {
		sum = sum + logGenPoisson(y[j], mu, alpha);
		//cout <<"j= "<<j<<" "<<logGenPoisson(y[j],mu,alpha)<<" "<<sum<<endl;
	}
	//exit(0);
	sum = sum - chromosomeLengthInBins * (log(1 - con));
	return (sum);
}//end of likelihood

double logGenPoisson(const double y, const double mu, const double alpha) {
	double result, ay, apart, tem01, tem02;

	ay = floor(y);
	apart = y - ay;
	tem01 = logIntGenPoisson((int) ay, mu, alpha);
	tem02 = exp(logIntGenPoisson((int) (ay + 1), mu, alpha) - tem01);
	result = tem01 + log((1 - apart) + apart * tem02);
	return result;
}//end of logGenPoisson

double logTrunc0GenPoisson(const double y, const double mu, const double alpha) {
	double result, ay, apart, tem01, tem02;

	if (y < 0.0001)
		return -100;
	else if (y < 1) {
		result = log(y) + logIntTrunc0GenPoisson(1, mu, alpha);
		return result;
	} else {
		ay = floor(y);
		apart = y - ay;
		tem01 = logIntTrunc0GenPoisson((int) ay, mu, alpha);
		tem02 = exp(logIntTrunc0GenPoisson((int) (ay + 1), mu, alpha) - tem01);
		result = tem01 + log((1 - apart) + apart * tem02);
		return result;
	}
}//end of logTrunc0GenPoisson

double logIntTrunc0GenPoisson(const int y, const double mu, const double alpha) {
	double result, con;
	if (y == 0) {
		return -100;
	} else {
		con = logIntGenPoisson(0, mu, alpha);
		result = logIntGenPoisson(y, mu, alpha) - log(1 - exp(con));
		return result;
	}
}//end of logIntTrunc0GenPoisson

double genPoisson(const double y, const double mu, const double alpha) {
	double result, ay, apart, tem01, tem02;

	ay = floor(y);
	apart = y - ay;
	tem01 = logGenPoisson((int) ay, mu, alpha);
	tem02 = logGenPoisson((int) (ay + 1), mu, alpha);
	result = (1 - apart) * exp(tem01) + apart * exp(tem02);
	return result;
}//end of genPoisson

double logIntGenPoisson(const int y, const double mu, const double alpha) {
	int j;
	double result, logyfac;

	if (y == 0)
		logyfac = 0;
	else if (y > 0) {
		logyfac = 0;
		for (j = 2; j <= y; j++)
			logyfac = logyfac + log((double) j);
	} else {
		cout << "error, y is negative. " << y << endl;
		exit(0);
	}
	result = log(mu) - log(1 + alpha * mu);
	result = y * result + (y - 1) * log(1 + alpha * y) - logyfac;
	result = result - mu * (1 + alpha * y) / (1 + alpha * mu);
	return result;
}//end of logIntGenPoisson

void readData(const char dataFileName[], vector<int> &order,
		vector<double> &data, int *nRow, double *totalcount) {
	int count = 0;
	int temOrder;
	double temVal;
	istringstream iss;
	string lineString;
	double sum = 0;

	ifstream inFile(dataFileName);
	if (!inFile) {
		cout << "Error opening input file" << dataFileName << endl;
		exit(0);
	}
	count = 0;
	sum = 0;
	while (inFile) {
		if (inFile) {
			getline(inFile, lineString);
			iss.clear();
			iss.str(lineString + " ");
			iss >> temOrder >> temVal;
			if (iss) {
				order.push_back(temOrder);
				data.push_back(temVal);
				//07/04/08				sum = sum + (double) floor(temVal);
				sum = sum + temVal;
			}//end of if
		}//end of if
		count++;
	}//end of while
	*nRow = count - 1;
	*totalcount = sum;
	cout << "There are " << *nRow << " nonzero counts." << endl;
}//end of readData

double unran(int *na, int *nb, int *nc) {
	double random;
	*na = (171 * (*na)) % 30269;
	*nb = (172 * (*nb)) % 30307;
	*nc = (170 * (*nc)) % 30323;
	random = (double) *na / 30269.0 + (double) *nb / 30307.0 + (double) *nc
			/ 30323.0;
	random = random - floor(random);
	return random;
}

double gammaln(double xx) {
	double ser, stp, tmp, x, y, cof[6], gam;
	int j;
	cof[0] = 76.18009172947146;
	cof[1] = -86.50532032941677;
	cof[2] = 24.01409824083091;
	cof[3] = -1.231739572450155;
	cof[4] = 0.1208650973866179 * 0.01;
	cof[5] = -0.5395239384953 * 0.00001;
	stp = 2.5066282746310005;
	x = xx;
	y = x;
	tmp = x + 5.5;
	tmp = (x + 0.5) * log(tmp) - tmp;
	ser = 1.000000000190015;
	for (j = 0; j < 6; j++) {
		y = y + 1.0;
		ser = ser + cof[j] / y;
	}
	gam = tmp + log(stp * ser / x);
	return gam;
}
