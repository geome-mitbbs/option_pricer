#ifndef GNUPLOT_H_
#define GNUPLOT_H_
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
using namespace std;
class Gnuplot {
public:
	Gnuplot();
	~Gnuplot();
	void operator ()(const string & command);
	void operator ()(const vector<double> x, const vector<double> y);
	void operator ()(const vector<double> x_list, const vector<vector<double>> y_list, double null_number = -1);
	void operator ()(const vector<double> x, const vector<vector<double>> y_list, const vector<vector<double>> z_list, double null_number = -1);
	// send any command to gnuplot
protected:
	FILE *gnuplotpipe;
};
Gnuplot::Gnuplot() {
	// with -persist option you will see the windows as your program ends
	//gnuplotpipe=_popen("gnuplot -persist","w");
	//without that option you will not see the window
	// because I choose the terminal to output files so I don't want to see the window
	gnuplotpipe = _popen("gnuplot", "w");
	if (!gnuplotpipe) {
		cerr << ("Gnuplot not found !");
	}
}
Gnuplot::~Gnuplot() {
	fprintf(gnuplotpipe, "exit\n");
	_pclose(gnuplotpipe);
}
void Gnuplot::operator()(const string & command) {
	fprintf(gnuplotpipe, "%s\n", command.c_str());
	fflush(gnuplotpipe);
	// flush is necessary, nothing gets plotted else
};
// plot x-y with lines
void Gnuplot::operator()(const vector<double> x, const vector<double> y){
	ofstream myfile;
	myfile.open("tmp_gnu_plot.dat");
	myfile << "X Y\n";
	for (int i = 0; i < x.size(); ++i)
		myfile << x[i] << " " << y[i] << endl;
	myfile.close();
	this->operator()("plot  \"tmp_gnu_plot.dat\" using 1:2 with lines");
};
// plot x-yn with lines
void Gnuplot::operator()(const vector<double> x, const vector<vector<double>> y_list, double null_number){
	ofstream myfile;
	myfile.open("tmp_gnu_plot.dat");
	int num_pairs = y_list.size();
	myfile << "X" << " ";
	for (int i = 0; i < num_pairs; ++i)
		myfile<< "Y" << i << " ";
	myfile << endl;
	
	for (int i = 0; i < x.size(); ++i){
		myfile << x[i] << " ";
		for (int j = 0; j < num_pairs; ++j){
			if (i < y_list[j].size()){
				if (y_list[j][i]!=null_number)
					myfile << y_list[j][i] << " ";
				else
					myfile << "." << " ";
			}
			else{
				myfile << "." << " ";
			}
		}
		myfile << endl;
	}
	myfile.close();

	stringstream ss;
	ss << "plot ";
	for (int i = 0; i < num_pairs; ++i){
		ss << "\"tmp_gnu_plot.dat\" using " << "1:" << i+2 << " with lines";
		if (i < num_pairs - 1)
			ss << ",";
	}
	this->operator()(ss.str());
};
// plot x-yn with lines, x-zn with dots.
void Gnuplot::operator()(const vector<double> x, const vector<vector<double>> y_list, const vector<vector<double>> z_list, double null_number){
	ofstream myfile;
	myfile.open("tmp_gnu_plot.dat");
	myfile << "X" << " ";
	for (int i = 0; i < y_list.size(); ++i)
		myfile << "Y" << i << " ";
	for (int i = 0; i < z_list.size(); ++i)
		myfile << "Z" << i << " ";
	myfile << endl;

	for (int i = 0; i < x.size(); ++i){
		myfile << x[i] << " ";
		for (int j = 0; j < y_list.size(); ++j){
			if (i < y_list[j].size()){
				if (y_list[j][i] != null_number)
					myfile << y_list[j][i] << " ";
				else
					myfile << "." << " ";
			}
			else{
				myfile << "." << " ";
			}
		}
		
		for (int j = 0; j < z_list.size(); ++j){
			if (i < z_list[j].size()){
				if (z_list[j][i] != null_number)
					myfile << z_list[j][i] << " ";
				else
					myfile << "." << " ";
			}
			else{
				myfile << "." << " ";
			}
		}

		myfile << endl;
	}
	myfile.close();

	stringstream ss;
	ss << "plot ";
	for (int i = 0; i < y_list.size(); ++i){
		ss << "\"tmp_gnu_plot.dat\" using " << "1:" << i + 2 << " with lines";
		ss << ",";
	}
	for (int i = 0; i < z_list.size(); ++i){
		ss << "\"tmp_gnu_plot.dat\" using " << "1:" << i + 2 + y_list.size();
		if (i < z_list.size() - 1)
			ss << ",";
	}
	this->operator()(ss.str());
};

#endif

