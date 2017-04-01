#include<vector>
#include<algorithm>
using namespace std;
#ifndef trigonal_matrix_h
#define trigonal_matrix_h
class trigonal_matrix
{
public:
	trigonal_matrix();
	trigonal_matrix(vector<double> diag, vector<double> upper, vector<double> lower);
	double& operator()(int i, int j);
	double size();
	void set(vector<double> diag, vector<double> upper, vector<double> lower);
private:
	vector<double> diag, upper, lower;
};
vector<double> solve_equation_trigonal(trigonal_matrix A, vector<double> b);
vector<double> operator*(vector<double>& left, double scaler);
#endif

