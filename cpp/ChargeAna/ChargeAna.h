#pragma once
#include <iostream>
#include <iostream>
#include <exception>
#define _USE_MATH_DEFINES 
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <vector>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <mutex> 
#include <algorithm>
#include<omp.h>
#include<atomic>

namespace py = pybind11;
using namespace std;
class ChargeAnalysis {
public:

	~ChargeAnalysis();
	ChargeAnalysis();

	//所有距离单位采用bohr
	vector<double> boxs;
	vector<int> nxyz;
	vector<vector<vector<double>>>rho;
	vector<double> dx3s;
	vector<vector<vector<int>>> gridtest;
	double dxvol=0;
	int natoms=0;
	bool debug = false;
	int threadnum = 1;
	mutex resmutex;


	void ChargeAnalysisInit(string rhofilenamein, vector<double> boxsin, int natomsin);
	void readrho(string rhofilename);
	void debug_gridtestadd(int i, int j, int k);
	double SumMethod_voro(vector<vector<vector<double>>> v3s, vector<vector<double>> vertex_positions);
	double SumMethod_cut(vector<double> atompos, double cutoff);
	vector<double>  CD_pdf(vector<double> atompos, double cutoff, double histlow = 0,int histbin=1000);;
	double getden(int i, int j, int k);
	double get_abs_distance(vector<double> atompos, vector<double> gridmidpoint);
	vector<vector<vector<vector< double>>>>  calculate_q_sumYdotCharge(vector <int> qs, vector<double> atompos, double cutoff, double histlow = 0, int histbin = 1000);
	void QLM(int l, int m, double theta, double phi, double& realYLM, double& imgYLM);
	void YLM(int l, int m, double theta, double phi, double& realYLM, double& imgYLM);
	double dfactorial(int l, int m);
	double PLM(int l, int m, double x);
	void convert_to_spherical_coordinates(double x, double y, double z, double& r, double& phi, double& theta);
};


