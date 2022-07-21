#include "ChargeAna.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdio.h>
#include "string.h"
#include <chrono>
#include <pybind11/stl.h>
#include <thread>
#include <mutex>

ChargeAnalysis::~ChargeAnalysis() {}
ChargeAnalysis::ChargeAnalysis() {}


void ChargeAnalysis::ChargeAnalysisInit(string rhofilenamein,vector<double> boxsin, int natomsin) {
	boxs = boxsin;
	natoms = natomsin;
	readrho(rhofilenamein);

}


void ChargeAnalysis::readrho(string rhofilename) {
	ifstream file;
	vector<stringstream>buff;
	vector<double> rho0;
	nxyz.resize(3);
	buff.resize(threadnum);
	file.open(rhofilename, ios::in|ios::binary);
	if (file){

		stringstream ss;
		string str;
		string temp;
		int index = 0;
		int datanumperline=0;
		ss.str("");
		if(getline(file,temp)){
			ss.clear();
			ss.str("");
			ss << temp;
			while (ss >> str)
			{
				nxyz[index] = int(atof(str.c_str()));
				index++;
				if (index == 3)
				{
					break;
				}
			}
		}
		if (debug) {
			gridtest = vector<vector<vector<int>>>(nxyz[0], vector<vector<int>>(nxyz[1], vector<int>(nxyz[2], 0)));
		}
		int linesbegin = file.tellg();
		if (getline(file, temp))
		{
			ss.clear();
			ss.str("");
			ss << temp;
			while (ss >> str)
			{
				datanumperline++;
			}
		}
		int linelength = int(file.tellg()) - linesbegin;
		file.seekg(0, file.end);
		unsigned long long  linesend = file.tellg();
		file.seekg(linesbegin, file.beg);
		int linenum = (linesend - linesbegin) / linelength + 1;
		int linenumperthread = linenum / threadnum + 1;
		unsigned long long p = linesbegin;
		char* temp2 = new char[linenumperthread * linelength];
		for (int i = 0; i < threadnum;i++)
		{
			file.read(temp2, linenumperthread * linelength);
			buff[i] << temp2;
			if (i != threadnum - 1)
			{
				file.seekg(p + linenumperthread * linelength, file.beg);
			}
			else{
				file.close();
			}
			p += linenumperthread * linelength;

		}
		delete[]temp2;

	
		rho = vector<vector<vector<double>>>(nxyz[0], vector<vector<double>>(nxyz[1], vector<double>(nxyz[2], 0)));
		omp_set_num_threads(threadnum);
		#pragma omp parallel
		{	
			stringstream ss3;
			string str3;
			string temp3;
			int threadid = omp_get_thread_num();
			//int threadid = 1;
			int iread = -1, jread = -1, kread = -1, all = (threadid * linenumperthread*datanumperline),count=0;
			jread = all / nxyz[0];
			iread = all % nxyz[0];
			kread = jread / nxyz[1];
			jread = jread % nxyz[1];
			while (getline(buff[threadid], temp3 )&& count< linenumperthread * datanumperline) {
				ss3.clear();
				ss3.str("");
				ss3 << temp3;
				while (ss3 >> str3 && count < linenumperthread * datanumperline) {
					rho[iread][jread][kread] = (atof(str3.c_str()));
					count++;
					//gridtest[iread][jread][kread]++;
					iread++;
					if (iread >= nxyz[1]) { iread = 0; jread++; }
					if (jread >= nxyz[0]) { jread = 0; kread++; }
				}
			}
		}
	}



	dx3s.resize(3);
	for (int i=0; i < 3; i++)
	{
		dx3s[i] = boxs[i] / nxyz[i];
	}
	dxvol = dx3s[0] * dx3s[1] * dx3s[2];


	//cout << test << "\n";
}

double ChargeAnalysis::SumMethod_voro(vector<vector<vector<double>>> v3s, vector<vector<double>> vertex_positions) {


	vector<double> minxyz, maxxyz;
	int facesnum;
	vector<vector<double>> X, Y, Z;


	vector<int> minxyzn(3, 0);
	vector<int> maxxyzn(3, 0);
	vector<double> centerofvetex(3,0);
	vector<vector<double>> nvector;
	vector<double> surfaceD,nvectornorm, centerDistance,gridmidpoint, centertopoint_project_nvector;
	double result = 0;
	bool kstartflag = 0;
	bool kendflag = 0;
	int kstart = 0;
	int kend = 0;
	double tempcpn = 0;
	int boollist=0;


	facesnum = int(v3s.size());
	X.resize(0);
	Y.resize(0);
	Z.resize(0);
	for (int n = 0; n < v3s.size(); n++) {
		X.emplace_back(vector<double>());
		Y.emplace_back(vector<double>());
		Z.emplace_back(vector<double>());

		for (int i = 0; i < 3; i++) {
			X[n].emplace_back(v3s[n][i][0]);
			Y[n].emplace_back(v3s[n][i][1]);
			Z[n].emplace_back(v3s[n][i][2]);
		}
	}
	minxyz.resize(3);
	maxxyz.resize(3);
	for (int j = 0; j < 3; j++) {
		maxxyz[j] = vertex_positions[0][j];
		minxyz[j] = vertex_positions[0][j];
	}


	for (int i = 0; i < vertex_positions.size(); i++) {
		for (int j = 0; j < 3; j++) {
			if (maxxyz[j] < vertex_positions[i][j]) {
				maxxyz[j] = vertex_positions[i][j];
			}
			if (minxyz[j] > vertex_positions[i][j]) {
				minxyz[j] = vertex_positions[i][j];
			}
		}
	}


	for (int i = 0; i < 3; i++)
	{
		minxyzn[i] = int(floor(minxyz[i] / dx3s[i]));
		maxxyzn[i] = int(floor(maxxyz[i] / dx3s[i]));
	}
	for (int n = 0; n < v3s.size(); n++)
	{
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				centerofvetex[j] += v3s[n][i][j];
			}
		}
	}
	for (int j = 0; j < 3; j++)
	{
		centerofvetex[j] /= v3s.size()*3;
	}
	for (int n = 0; n < v3s.size(); n++)
	{
		vector<double>temp;
		double distemp=0;
		temp.emplace_back(Y[n][1] * Z[n][2] - Y[n][1] * Z[n][0] - Y[n][0] * Z[n][2] - Y[n][2] * Z[n][1] + Y[n][0] * Z[n][1] + Y[n][2] * Z[n][0]);
		temp.emplace_back(X[n][2] * Z[n][1] - X[n][0] * Z[n][1] - X[n][2] * Z[n][0] - X[n][1] * Z[n][2] + X[n][1] * Z[n][0] + X[n][0] * Z[n][2]);
		temp.emplace_back(X[n][1] * Y[n][2] - X[n][1] * Y[n][0] - X[n][0] * Y[n][2] - X[n][2] * Y[n][1] + X[n][2] * Y[n][0] + X[n][0] * Y[n][1]);
		nvector.emplace_back(temp);
		//surfaceD.emplace_back(inner_product(begin(nvector[n]), end(nvector[n]),begin(v3s[n][0]), 0));
		//nvectornorm.emplace_back(inner_product(begin(nvector[n]), end(nvector[n]), begin(nvector[n]), 0));
		surfaceD.emplace_back(0);
		nvectornorm.emplace_back(0);
		for (int m = 0; m < 3; m++) {
			surfaceD[n] += -nvector[n][m] * v3s[n][0][m];
			nvectornorm[n] += nvector[n][m] * nvector[n][m];
		}
		nvectornorm[n] = sqrt(nvectornorm[n]);
		for (int i = 0; i < 3; i++)
		{
			distemp += nvector[n][i] * centerofvetex[i];
		}
		distemp = (distemp + surfaceD[n]) / nvectornorm[n];

		if (distemp > 0) {
			for (int i = 0; i < 3; i++) {
				nvector[n][i] *= -1;
			}
			surfaceD[n] *= -1;
		}
		else if(distemp <= 0) {
			distemp *= -1;
		}
		centerDistance.emplace_back(distemp);
	}
	gridmidpoint.resize(3);
	for (int i = minxyzn[0]; i < maxxyzn[0] + 1; i++) {
		for (int j = minxyzn[1]; j < maxxyzn[1] + 1; j++) {
			kstartflag = 0;
			kendflag = 0;
			for (int k = minxyzn[2]; k < maxxyzn[2] + 1; k++) {

				gridmidpoint[0] = i * dx3s[0] + dx3s[0] / 2;
				gridmidpoint[1] = j * dx3s[1] + dx3s[1] / 2;
				gridmidpoint[2] = k * dx3s[2] + dx3s[2] / 2;
				centertopoint_project_nvector.resize(0);
				boollist = 0;
				for (int n = 0; n < v3s.size(); n++)
				{
					tempcpn=0;
					for (int m = 0; m < 3; m++)
					{
						tempcpn += (gridmidpoint[m] - centerofvetex[m]) * nvector[n][m];
					}
					tempcpn /= nvectornorm[n];

					centertopoint_project_nvector.emplace_back(tempcpn);
					if (centertopoint_project_nvector[n] <= centerDistance[n]) {
						boollist++;
					}
				}
				if (boollist == v3s.size())
				{
					result += getden(i, j, k) * dxvol;
					debug_gridtestadd(i, j, k);
					if (k + 1 <= maxxyzn[2])
					{
						kstart = int(k + 1);
						kstartflag = true;
						break;
					}
				}



			}
			if (kstartflag == true) {
				for (int k = maxxyzn[2]; k >= kstart; k--)
				{
					gridmidpoint[0] = i * dx3s[0] + dx3s[0] / 2;
					gridmidpoint[1] = j * dx3s[1] + dx3s[1] / 2;
					gridmidpoint[2] = k * dx3s[2] + dx3s[2] / 2;
					centertopoint_project_nvector.resize(0);
					boollist = 0;
					for (int n = 0; n < v3s.size(); n++)
					{
						tempcpn = 0;
						for (int m = 0; m < 3; m++)
						{
							tempcpn += (gridmidpoint[m] - centerofvetex[m]) * nvector[n][m];
						}
						tempcpn /= nvectornorm[n];

						centertopoint_project_nvector.emplace_back(tempcpn);
						if (centertopoint_project_nvector[n] <= centerDistance[n]) {
							boollist++;
						}
					}
					if (boollist == v3s.size())
					{
						result += getden(i, j, k) * dxvol;
						debug_gridtestadd(i, j, k);
						if (k - 1 >= kstart)
						{
							kend = int(k - 1);
							kendflag = true;
							break;
						}
					}
				}
				if (kendflag == true) {
					for (int k = kstart; k < kend + 1;k++) {
						result += getden(i, j, k) * dxvol;
						debug_gridtestadd(i, j, k);
					}
				}
			}
		}
		
	}
	return result;
}

double ChargeAnalysis::SumMethod_cut(vector<double> atompos, double cutoff) {

	vector<double> minxyz, maxxyz;
	vector<int> minxyzn(3, 0);
	vector<int> maxxyzn(3, 0);


	vector<double> gridmidpoint;
	double result = 0;
	bool kstartflag = 0;
	bool kendflag = 0;
	int kstart = 0;
	int kend = 0;

	minxyz.resize(3);
	maxxyz.resize(3);
	for (int j = 0; j < 3; j++) {
		maxxyz[j] = atompos[j] + cutoff;
		minxyz[j] = atompos[j] - cutoff;
	}


	for (int i = 0; i < 3; i++)
	{
		minxyzn[i] = int(floor(minxyz[i] / dx3s[i]));
		maxxyzn[i] = int(floor(maxxyz[i] / dx3s[i]));
	}

	gridmidpoint.resize(3);
	for (int i = minxyzn[0]; i < maxxyzn[0] + 1; i++) {
		for (int j = minxyzn[1]; j < maxxyzn[1] + 1; j++) {
			kstartflag = 0;
			kendflag = 0;
			for (int k = minxyzn[2]; k < maxxyzn[2] + 1; k++) {

				gridmidpoint[0] = i * dx3s[0] + dx3s[0] / 2;
				gridmidpoint[1] = j * dx3s[1] + dx3s[1] / 2;
				gridmidpoint[2] = k * dx3s[2] + dx3s[2] / 2;

				if (get_abs_distance(atompos,gridmidpoint) <= cutoff)
				{
					result += getden(i, j, k) * dxvol;
					debug_gridtestadd(i, j, k);
					if (k + 1 <= maxxyzn[2])
					{
						kstart = int(k + 1);
						kstartflag = true;
						break;
					}
				}



			}
			if (kstartflag == true) {
				for (int k = maxxyzn[2]; k >= kstart; k--)
				{
					gridmidpoint[0] = i * dx3s[0] + dx3s[0] / 2;
					gridmidpoint[1] = j * dx3s[1] + dx3s[1] / 2;
					gridmidpoint[2] = k * dx3s[2] + dx3s[2] / 2;
					
					if (get_abs_distance(atompos, gridmidpoint) <= cutoff)
					{
						result += getden(i, j, k) * dxvol;
						debug_gridtestadd(i, j, k);
						if (k - 1 >= kstart)
						{
							kend = int(k - 1);
							kendflag = true;
							break;
						}
					}
				}
				if (kendflag == true) {
					for (int k = kstart; k < kend + 1; k++) {
						result += getden(i, j, k) * dxvol;
						debug_gridtestadd(i, j, k);
					}
				}
			}
		}

	}
	return result;
}

vector<double>  ChargeAnalysis::CD_pdf(vector<double> atompos, double cutoff, double histlow, int histbin) {

	vector<vector<double>> res;
	vector<double>		 resall	(histbin, 0);
	vector<vector<int>>	gridindex;
	vector<int>			gridtemp(3, 0);
	vector<double>		minxyz	(3, 0);
	vector<double>		maxxyz	(3, 0);
	vector<int>			minxyzn	(3, 0);
	vector<int>			maxxyzn	(3, 0);
	const double		deltacut = (cutoff - histlow) / histbin;
	double				r = 0;


	res.resize(threadnum);
	for (int t=0;t<threadnum;t++){
		res[t].resize(histbin,0);}
	for (int j = 0; j < 3; j++) {
		maxxyz[j] = atompos[j] + cutoff;
		minxyz[j] = atompos[j] - cutoff;}
	for (int i = 0; i < 3; i++)
	{
		minxyzn[i] = int(floor(minxyz[i] / dx3s[i]));
		maxxyzn[i] = int(floor(maxxyz[i] / dx3s[i]));}


	//for (int i = minxyzn[0]; i < maxxyzn[0] + 1; i++) {
	for (int j = minxyzn[1]; j < maxxyzn[1] + 1; j++) {
		for (int k = minxyzn[2]; k < maxxyzn[2] + 1; k++) {
			//gridtemp[0] = i;
			gridtemp[1] = j;
			gridtemp[2] = k;
			gridindex.emplace_back(gridtemp);
		}
	}
	//}	

	omp_set_num_threads(threadnum);
#	pragma omp parallel
	{	
		vector<double> gridmidpoint(3, 0);
		double d = 0;
		int threadid = omp_get_thread_num();
		//int i = -1;
		int	j = -1;
		int	k = -1;
#pragma omp for 
		for (int g = 0; g < gridindex.size(); g++) {
#pragma omp  simd
			for (int i = minxyzn[0]; i < maxxyzn[0] + 1; i++) {
				//i = gridindex[g][0];
				j = gridindex[g][1];
				k = gridindex[g][2];
				gridmidpoint[0] = i * dx3s[0] + dx3s[0] / 2;
				gridmidpoint[1] = j * dx3s[1] + dx3s[1] / 2;
				gridmidpoint[2] = k * dx3s[2] + dx3s[2] / 2;
				d = get_abs_distance(atompos, gridmidpoint);
				if (d <= cutoff && d >= histlow) {
					res[threadid][int(floor((d - histlow) / deltacut))]+= getden(i, j, k);
					debug_gridtestadd(i, j, k);
				}
			}
		}
	}
	
	for (int i = 0; i < histbin; i++) {
		for (int j = 0; j < threadnum; j++) {
			resall[i] += res[j][i];
		}
	}
	for (int i = 0; i < histbin; i++) {
		r = 0.5 * deltacut + histlow + i * deltacut;
		resall[i] = resall[i] * dxvol /( 4 * M_PI * r*r * deltacut);
	}
	return resall;
}

double ChargeAnalysis::get_abs_distance(vector<double> atompos, vector<double> gridmidpoint){
	double abs,diffx,diffy,diffz;


	diffx = gridmidpoint[0] - atompos[0];
	diffy = gridmidpoint[1] - atompos[1];
	diffz = gridmidpoint[2] - atompos[2];

	if (diffx >  boxs[0] / 2.0) { diffx -= boxs[0]; };
	if (diffx < -boxs[0] / 2.0) { diffx += boxs[0]; };
	if (diffy >  boxs[1] / 2.0) { diffy -= boxs[1]; };
	if (diffy < -boxs[1] / 2.0) { diffy += boxs[1]; };
	if (diffz >  boxs[2] / 2.0) { diffz -= boxs[2]; };
	if (diffz < -boxs[2] / 2.0) { diffz += boxs[2]; };
	abs = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);

	return abs;

}

double ChargeAnalysis::getden(int i, int j, int k) {
	vector<int>i3s = { i,j,k };
	for (int m=0; m < 3; m++)
	{
		if ( i3s[m]  < 0)
		{
			i3s[m] = i3s[m] % nxyz[m];
			if(i3s[m]<0)
			{
				i3s[m] += nxyz[m];
			}
		}
		else if (i3s[m] > nxyz[m] - 1)
		{
			i3s[m] = i3s[m] % nxyz[m];
		}
			
	}
	return rho[i3s[0]][i3s[1]][i3s[2]];
}

void ChargeAnalysis::debug_gridtestadd(int i, int j, int k) {

	if (debug && threadnum==1){
		vector<int>i3s = { i,j,k };
		for (int m = 0; m < 3; m++)
		{
			if (i3s[m] < 0)
			{
				i3s[m] = i3s[m] % nxyz[m];
				if (i3s[m] < 0)
				{
					i3s[m] += nxyz[m];
				}
			}
			else if (i3s[m] > nxyz[m] - 1)
			{
				i3s[m] = i3s[m] % nxyz[m];
			}

		}
		gridtest[i3s[0]][i3s[1]][i3s[2]]++;
	}
	return;
}
/*
int main() {
	vector<double> box = { 5,5,5 };
	string filename = "atlas.den";
	ChargeAnalysis(filename,box, 10);
}
*/

vector<vector<vector<vector< double>>>>  ChargeAnalysis::calculate_q_sumYdotCharge(vector <int> qs, vector<double> atompos, double cutoff, double histlow, int histbin) {

	
	vector<vector<vector<vector<vector< double>>>>> res;
	vector<vector<vector<vector< double>>>>			resall;
	vector<vector<int>>								gridindex;
	vector<int>										gridtemp(3,0);
	vector<double>									minxyz(3, 0);
	vector<double>									maxxyz(3, 0);
	vector<int>										minxyzn(3, 0);
	vector<int>										maxxyzn(3, 0);
	const double deltacut = (cutoff - histlow) / histbin;

	res.resize(threadnum);
	for (int t=0;t<threadnum;t++){
		res[t].resize(3);
		for (int i = 0; i < 3; i++){
			res[t][i].resize(qs.size());
			for (int j = 0; j < qs.size(); j++){
				res[t][i][j].resize(1 + 2 * qs[j]);
				for (int k = 0; k < 1 + 2 * qs[j]; k++){
					res[t][i][j][k].resize(histbin);
				}
			}
		}
	}

	resall.resize(3);
	for (int i = 0; i < 3; i++) {
		resall[i].resize(qs.size());
		for (int j = 0; j < qs.size(); j++) {
			resall[i][j].resize(1 + 2 * qs[j]);
			for (int k = 0; k < 1 + 2 * qs[j]; k++) {
				resall[i][j][k].resize(histbin);
			}
		}
	}
	
	for (int j = 0; j < 3; j++) {
		maxxyz[j] = atompos[j] + cutoff;
		minxyz[j] = atompos[j] - cutoff;
	}
	for (int i = 0; i < 3; i++)
	{
		minxyzn[i] = int(floor(minxyz[i] / dx3s[i]));
		maxxyzn[i] = int(floor(maxxyz[i] / dx3s[i]));
	}
	//for (int i = minxyzn[0]; i < maxxyzn[0] + 1; i++) {
		for (int j = minxyzn[1]; j < maxxyzn[1] + 1; j++) {
			for (int k = minxyzn[2]; k < maxxyzn[2] + 1; k++) {
				//gridtemp[0] = i;
				gridtemp[1] = j;
				gridtemp[2] = k;
				gridindex.emplace_back(gridtemp);
			}
		}
	//}	
	omp_set_num_threads(threadnum);
	#pragma omp parallel
	{
		vector<double> gridmidpoint(3, 0);
		vector <double>theta_phi(2, 0);
		double realYLM, imgYLM;
		double den;
		double d = 0;
		int threadid = omp_get_thread_num();
		//int i = -1;
		int	j = -1;
		int	k = -1;
		int index;
		int q;


#pragma omp for 
		for (int g = 0; g < gridindex.size(); g++) {
#pragma omp  simd
			for (int i = minxyzn[0]; i < maxxyzn[0] + 1; i++) {
				//i = gridindex[g][0];
				j = gridindex[g][1];
				k = gridindex[g][2];
				gridmidpoint[0] = i * dx3s[0] + dx3s[0] / 2;
				gridmidpoint[1] = j * dx3s[1] + dx3s[1] / 2;
				gridmidpoint[2] = k * dx3s[2] + dx3s[2] / 2;
				d = get_abs_distance(atompos, gridmidpoint);
				if (d <= cutoff && d >= histlow) {
					theta_phi[0] = acos((gridmidpoint[2] - atompos[2]) / d);//acos z/r
					theta_phi[1] = atan2((gridmidpoint[1] - atompos[1]), (gridmidpoint[0] - atompos[0]));//atan2(y,x)
					den = getden(i, j, k);
					for (int tq = 0; tq < qs.size(); tq++) {
						q = qs[tq];
						for (int mi = -q; mi < q + 1; mi++) {
							QLM(q, mi, theta_phi[0], theta_phi[1], realYLM, imgYLM);
							index = int(floor((d - histlow) / deltacut));
							res[threadid][0][tq][mi + q][index] += (den * realYLM);
							res[threadid][1][tq][mi + q][index] += (den * imgYLM);
							res[threadid][2][tq][mi + q][index] += (den);
							debug_gridtestadd(i, j, k);
						}
					}
				}
			}
		}
	}
	
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < qs.size(); j++) {
			for (int k = 0; k < 2 * qs[j] + 1; k++) {
				for (int m = 0; m < histbin; m++) {
					for (int t = 0; t < threadnum; t++)
					{
						resall[i][j][k][m] += res[t][i][j][k][m] * dxvol;
					}
				}
			}
		}
	}
	return resall;
}


void ChargeAnalysis::QLM(int l, int m, double theta, double phi, double& realYLM, double& imgYLM) {

	realYLM = 0.0;
	imgYLM = 0.0;
	if (m < 0) {
		YLM(l, abs(m), theta, phi, realYLM, imgYLM);
		realYLM = pow(-1.0, m) * realYLM;
		imgYLM = pow(-1.0, m + 1) * imgYLM;
	}
	else {
		YLM(l, m, theta, phi, realYLM, imgYLM);
	}
}

void ChargeAnalysis::YLM(int l, int m, double theta, double phi, double& realYLM, double& imgYLM) {

	double factor;
	double m_PLM;
	m_PLM = PLM(l, m, cos(theta));
	factor = ((2.0 * double(l) + 1.0) / (4.0 * M_PI)) * dfactorial(l, m);
	realYLM = sqrt(factor) * m_PLM * cos(double(m) * phi);
	imgYLM = sqrt(factor) * m_PLM * sin(double(m) * phi);
}

double ChargeAnalysis::PLM(int l, int m, double x) {

	double fact, pll, pmm, pmmp1, somx2;
	int i, ll;
	pll = 0.0;
	if (m < 0 || m > l || fabs(x) > 1.0)
		cerr << "impossible combination of l and m" << "\n";
	pmm = 1.0;
	if (m > 0) {
		somx2 = sqrt((1.0 - x) * (1.0 + x));
		fact = 1.0;
		for (i = 1; i <= m; i++) {
			pmm *= -fact * somx2;
			fact += 2.0;
		}
	}

	if (l == m)
		return pmm;
	else {
		pmmp1 = x * (2 * m + 1) * pmm;
		if (l == (m + 1))
			return pmmp1;
		else {
			for (ll = m + 2; ll <= l; ll++) {
				pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
				pmm = pmmp1;
				pmmp1 = pll;
			}
			return pll;
		}
	}
}

double ChargeAnalysis::dfactorial(int l, int m) {

	double fac = 1.00;
	for (int i = 0; i < 2 * m; i++) {
		fac *= double(l + m - i);
	}
	return (1.00 / fac);
}
