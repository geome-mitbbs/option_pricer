#include "stdafx.h"
#include<vector>
#include<stdexcept>
#include<sstream>
#include<iostream>
#include "Lin_Algebra.h"
using namespace std;

trigonal_matrix::trigonal_matrix(){}
trigonal_matrix::trigonal_matrix(vector<double> diag, vector<double> upper, vector<double> lower) :diag(diag), upper(upper), lower(lower){}
void trigonal_matrix::set(vector<double> diag, vector<double> upper, vector<double> lower)
{
	this->diag = diag;
	this->upper = upper;
	this->lower = lower;
}

double& trigonal_matrix::operator()(int i, int j)
{
	if((j - i)>1 || (j-i) <-1 )
	{
		stringstream msg;
		msg << i << " and " << j << " out of triagonal matrix range" << endl;
		throw invalid_argument(msg.str());
	};

	if (j == i)
	{
		return diag[i];
	}
	else if (j == i - 1)
	{
		return lower[j];
	}
	else
	{
		return upper[i];
	}
}

double trigonal_matrix::size()
{
	return diag.size();
}

vector<double> solve_equation_trigonal(trigonal_matrix A, vector<double> b)
{
	double ratio;
	for (int base_row = 0; base_row < A.size() - 1; ++base_row)
	{
		ratio = A(base_row+1,base_row) / A(base_row,base_row);
		A(base_row+1,base_row) = 0;
		A(base_row + 1,base_row + 1) -= A(base_row,base_row + 1) * ratio;
		b[base_row + 1] -= b[base_row] * ratio;
	}

	for (int base_row = A.size() - 1; base_row>0; --base_row)
	{
		ratio = A(base_row - 1,base_row) / A(base_row,base_row);
		b[base_row - 1] -= b[base_row] * ratio;
	}

	for (int base_row = 0; base_row < A.size(); ++base_row)
	{
		b[base_row] /= A(base_row,base_row);
	}
	
	return b;
}

vector<double> operator*(vector<double>& left, double scaler)
{
	vector<double> ret(left.size());
	for (int i = 0; i < left.size(); ++i)
		ret[i] = left[i] * scaler;
	return ret;
}
