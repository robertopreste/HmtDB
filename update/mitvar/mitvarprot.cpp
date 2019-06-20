#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <vector>
#include "./includes/aminoacids.h"
#include "./includes/mtrev_100.h"
#include "./includes/mtrev_120.h"
#include "./includes/mtrev_160.h"
#include "./includes/mtrev_200.h"

using namespace std;

#define CUTOFF				5.7
#define GAP					'-'
#define GAP2					'*'
#define DISTANCE_DELTA		0.0001
#define RATE_FILE_EXTENSION	".rates.tdt"

typedef struct Sequence {
	string name;
	string data;
	Sequence(string& inName, string& inData) {
		name.assign(inName);
		data.assign(inData);
	}
} Sequence;

typedef struct Distance {
	Sequence* a;
	Sequence* b;
	double identity;
	double distance;
	Distance(Sequence* inA, Sequence* inB, double inIdent, double inDist) 
		: a(inA), b(inB) {
			identity = inIdent;
			distance = inDist;
	}
} Distance;

/******************************************************************/
/***************** Prototypes ************************************/
/******************************************************************/

void ParseInputFile(string& inFileName, vector<Sequence>& ioSequences);
void ProcessData(vector<Sequence>& inSequences, const double inCutoff,
				 vector<double>& ioRates, int& outValues, double& outMaxRate,
				 double& outMeanRate);
void CalcSitesRate(vector<Sequence>& inSequences, vector<Distance>& inDistances,
				   double inCutoff, vector<double>& ioRates, double& outMaxRate,
				   double& outMeanRate);
double SubstitutionWeight(char inA, char inB, double inDist);
void PrintRateResults(string& inFileName, vector<double>& inRates, int inValueCount,
					  double inMaxRate, double inMeanRate);
void PrintComp(vector<Sequence>& inSequences);

/******************************************************************/
/***************** Functions **************************************/
/******************************************************************/

int main(int argc, char* argv[]) {

	string fileName;
	string outputName;
	fileName.assign(argv[1] + 4, strlen(argv[1]) - 4);
	outputName.assign(argv[2] + 5, strlen(argv[2]) - 5);
	double cutoff = CUTOFF;
	
	vector<Sequence> sequences;
	ParseInputFile(fileName, sequences);
	if (sequences.empty()) {
		return 1;
	}
	
	//cout << sequences[9].name << endl << sequences[9].data << endl << sequences.size();
		
	
	vector<double> rates(sequences[0].data.size());
	//vector<double> ratescoding;
	double maxRate;
	double meanRate;

	int valueCount;
	ProcessData(sequences, cutoff, rates, valueCount, maxRate, meanRate);
	PrintRateResults(outputName, rates, valueCount, maxRate, meanRate);

	PrintComp(sequences);

	return 0;

}

void ParseInputFile(string& inFileName, vector<Sequence>& ioSequences) {
	ifstream is(inFileName.c_str());
	if (!is) {
		cout << "Problems reading file " << inFileName;
		return;
	}

	cout << "Loading Alignement sequences...";
	cout.flush();
	char line[BUFSIZ];
	string name;
	//////////////
	string temp;
	unsigned int t = 0;
	//////////////
	string sequence;

	sequence.reserve(BUFSIZ);
	
	while (t == 0) {
	
		is.getline(line, BUFSIZ);
	
		if (is.eof()) {

			ioSequences.push_back(Sequence(name, sequence));
			//cout << "read : " << name << "\n";
			t++;

		}
	
		if (line[0] == '>') {
			temp.assign(line+1,strlen(line)-1);
			//cout << temp << "\n";
			if (name != temp){
			
				if (!name.empty() and !sequence.empty()) {

					ioSequences.push_back(Sequence(name, sequence));
					//cout << "read : " << name << "\n";

				}

				name.assign(temp);
				sequence = "";

			}

		} else if (line != ""){

			sequence.append(line, strlen(line));
	
		}

	}
	cout << "done" << endl;

}

void ProcessData(vector<Sequence>& inSequences, const double inCutoff, 
				 vector<double>& ioRates, int& outValues,
				 double& outMaxRate,double& outMeanRate){
	vector<Distance> distances;
	const int seqCount = inSequences.size();
	cout << "Number of Sequences: " << seqCount << " of Length: " << inSequences[0].data.size() << endl;
	cout << "Calculating Distances...";
	cout.flush();

	for (int i = 0; i < seqCount - 1; i++) {
		for (int j = i + 1; j < seqCount; j++) {
			double p = 0;
			double P = 0;
			int seqLen = inSequences[i].data.size();
			double pLen = (double)seqLen;
			for (unsigned int k = 0; k < seqLen; k++) {
				const char x = inSequences[i].data[k];
				const char y = inSequences[j].data[k];
				if (x != y) {
					p++;
				} else if (x == y == '-' | x == y == '*') {
					pLen--;
				}
			}
			
			P = (pLen - p) / pLen;

			if (P > 0.85) {
				P = 0.85;
			 	const double dist = (pLen - p) / pLen;
				const double T = -log(1-P-((P*P)/5));
				distances.push_back(Distance(&inSequences[i], &inSequences[j], T, dist));
				//cout << inSequences[i].name << "\t" << inSequences[j].name << "\tp = \t" << p << "\tP = \t" << P << "\tdist = \t" << T << endl;
			}
		}
	}
	
	outValues = distances.size();
	cout << "done\n";
	
	CalcSitesRate(inSequences, distances, inCutoff, ioRates, outMaxRate, outMeanRate);

}

void CalcSitesRate(vector<Sequence>& inSequences, vector<Distance>& inDistances,
				   double inCutoff, vector<double>& ioRates, double& outMaxRate,
				   double& outMeanRate) {
	int seqLen = inSequences[0].data.size();
	double totRate = 0;
	double distmed = 0;
	outMaxRate = 0;
	char gaps[inSequences[0].data.size()];
	cout << "Calculating Rates...";
	cout.flush();
	// Calcolo Kmedia
	
	for (unsigned int f = 0; f < inDistances.size(); f++){
	
		distmed += inDistances[f].identity + DISTANCE_DELTA;
	
	}
	
	distmed = distmed / inDistances.size();

	// Fine Calcolo

	for (int i = 0; i < seqLen; i++) {
		for (unsigned int j = 0; j < inDistances.size(); j++) {
			const double w = SubstitutionWeight(inDistances[j].a->data[i], inDistances[j].b->data[i], inDistances[j].distance);
			const double d = inDistances[j].identity + DISTANCE_DELTA;
			
			// Controlla i gap
			if ((inDistances[j].a->data[i] == GAP or inDistances[j].b->data[i] == GAP) or (inDistances[j].a->data[i] == GAP2 or inDistances[j].b->data[i] == GAP2)){
				gaps[i] = 'Y';
			}
			/////////////////////
			const double adRate = w / d;
			
			////////////////////////
			if (adRate != 0){
				//cout << adRate << endl;
			}
			
			//cout << inDistances[j].distance << endl;
			////////////////////////
			
			ioRates[i] += adRate;
		}
		// vi = vi + 2 / Kmedia
		if (gaps[i] == 'Y'){
		
			//cout << "Gap on site: " << i+1 << "\n";
			ioRates[i] += 2 / distmed;

		}
		
		//Calcolo absolute rate
		ioRates[i] = ioRates[i] / inDistances.size();
		
		//////////////////////////
		outMaxRate = max(ioRates[i], outMaxRate);
		totRate += ioRates[i];
		
		//cout << (ioRates[i]/outMaxRate) << endl;

		
	}
	outMeanRate = totRate / ioRates.size();
	///
	cout << "done\n";
	cout.setf(ios::fixed, ios::floatfield);
	cout << "MaxRate Coding: " << outMaxRate << endl;
	///

}

double SubstitutionWeight(char inA, char inB, double inDist){
	
	unsigned int idx1 = 0;
	unsigned int idx2 = 0;
	
	if (inA == GAP or inB == GAP or inA == GAP2 or inB == GAP2){

		return 0;
	}
	
	switch (inA){
		case 'A':
		case 'a':
			idx1 = AA_A;
			break;
		case 'C':
		case 'c':
			idx1 = AA_C;
			break;
		case 'D':
		case 'd':
			idx1 = AA_D;
			break;
		case 'E':
		case 'e':
			idx1 = AA_E;
			break;
		case 'F':
		case 'f':
			idx1 = AA_F;
			break;
		case 'G':
		case 'g':
			idx1 = AA_G;
			break;
		case 'H':
		case 'h':
			idx1 = AA_H;
			break;
		case 'I':
		case 'i':
			idx1 = AA_I;
			break;
		case 'K':
		case 'k':
			idx1 = AA_K;
			break;
		case 'L':
		case 'l':
			idx1 = AA_L;
			break;
		case 'M':
		case 'm':
			idx1 = AA_M;
			break;
		case 'N':
		case 'n':
			idx1 = AA_N;
			break;
		case 'P':
		case 'p':
			idx1 = AA_P;
			break;
		case 'Q':
		case 'q':
			idx1 = AA_Q;
			break;
		case 'R':
		case 'r':
			idx1 = AA_R;
			break;
		case 'S':
		case 's':
			idx1 = AA_S;
			break;
		case 'T':
		case 't':
			idx1 = AA_T;
			break;
		case 'V':
		case 'v':
			idx1 = AA_V;
			break;
		case 'W':
		case 'w':
			idx1 = AA_W;
			break;
		case 'Y':
		case 'y':
			idx1 = AA_Y;
			break;
	}
	
	switch (inB){
		case 'A':
		case 'a':
			idx2 = AA_A;
			break;
		case 'C':
		case 'c':
			idx2 = AA_C;
			break;
		case 'D':
		case 'd':
			idx2 = AA_D;
			break;
		case 'E':
		case 'e':
			idx2 = AA_E;
			break;
		case 'F':
		case 'f':
			idx2 = AA_F;
			break;
		case 'G':
		case 'g':
			idx2 = AA_G;
			break;
		case 'H':
		case 'h':
			idx2 = AA_H;
			break;
		case 'I':
		case 'i':
			idx2 = AA_I;
			break;
		case 'K':
		case 'k':
			idx2 = AA_K;
			break;
		case 'L':
		case 'l':
			idx2 = AA_L;
			break;
		case 'M':
		case 'm':
			idx2 = AA_M;
			break;
		case 'N':
		case 'n':
			idx2 = AA_N;
			break;
		case 'P':
		case 'p':
			idx2 = AA_P;
			break;
		case 'Q':
		case 'q':
			idx2 = AA_Q;
			break;
		case 'R':
		case 'r':
			idx2 = AA_R;
			break;
		case 'S':
		case 's':
			idx2 = AA_S;
			break;
		case 'T':
		case 't':
			idx2 = AA_T;
			break;
		case 'V':
		case 'v':
			idx2 = AA_V;
			break;
		case 'W':
		case 'w':
			idx2 = AA_W;
			break;
		case 'Y':
		case 'y':
			idx2 = AA_Y;
			break;
	}

	if (inDist > 0.90){

		return mtrev_100[idx1][idx2];		//mtrev100 ==> Blosum90
	
	} else if (0.90 >= inDist >0.80) {
		
		return mtrev_120[idx1][idx2];		//mtrev120 ==> Blosum80
	
	} else if (0.80 >= inDist > 0.60) {
		
		return mtrev_160[idx1][idx2];		//mtrev160 ==> Blosum60
	
	} else{
	
		return mtrev_200[idx1][idx2];		//mtrev200 ==> Blosum52
	
	}

}

void PrintRateResults(string& inFileName, vector<double>& inRates, int inValueCount,
					  double inMaxRate, double inMeanRate) {
    ofstream fos((inFileName + RATE_FILE_EXTENSION).c_str());
    if (!fos) {
    	cout << "Problems writing to file " << inFileName << ". Print to standard output";
	}
	
	const int rateCount = inRates.size();
	ostream& os = fos != NULL ? fos : cout;
	
	os << "site\t rel_rate\t normalized_rate\n";
	for (int i = 0; i < rateCount; i++) {
		os << (i + 1) << '\t' << (inRates[i] / inMaxRate) << '\t';
		os << (inRates[i] / inValueCount / inMeanRate) << '\n';
	}
}


void PrintComp(vector<Sequence>& inSequences){
	
	int seqLen = inSequences[0].data.size();	
	
	ofstream fos("debug.log");
    	if (!fos) {
    		cout << "Problems writing to file debug.log. Print to standard output";
	}
	
	ostream& os = fos != NULL ? fos : cout;

	os << "site;A;C;D;E;F;G;H;I;K;L;M;N;P;Q;R;S;T;V;W;Y;-;\n";

	for (int i = 0; i < seqLen; i++) {
		
		vector<int> comp(21);
		
		for (unsigned int j = 0; j < inSequences.size(); j++) {
			
		switch (inSequences[j].data[i]){
		case 'A':
		case 'a':
			comp[0]++;
			break;
		case 'C':
		case 'c':
			comp[1]++;
			break;
		case 'D':
		case 'd':
			comp[2]++;
			break;
		case 'E':
		case 'e':
			comp[3]++;
			break;
		case 'F':
		case 'f':
			comp[4]++;
			break;
		case 'G':
		case 'g':
			comp[5]++;
			break;
		case 'H':
		case 'h':
			comp[6]++;
			break;
		case 'I':
		case 'i':
			comp[7]++;
			break;
		case 'K':
		case 'k':
			comp[8]++;
			break;
		case 'L':
		case 'l':
			comp[9]++;
			break;
		case 'M':
		case 'm':
			comp[10]++;
			break;
		case 'N':
		case 'n':
			comp[11]++;
			break;
		case 'P':
		case 'p':
			comp[12]++;
			break;
		case 'Q':
		case 'q':
			comp[13]++;
			break;
		case 'R':
		case 'r':
			comp[14]++;
			break;
		case 'S':
		case 's':
			comp[15]++;
			break;
		case 'T':
		case 't':
			comp[16]++;
			break;
		case 'V':
		case 'v':
			comp[17]++;
			break;
		case 'W':
		case 'w':
			comp[18]++;
			break;
		case 'Y':
		case 'y':
			comp[19]++;
			break;
		case '-':
			comp[20]++;
			break;
		}
			
		}
		
		os << i << ";";
		
		for (int q = 0; q < 21; q++){
			os << comp[q] << ";";
		}
		
		os << endl;
		
	}
}
