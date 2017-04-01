#include "stdafx.h"
#include<vector>
#include<iostream>
#include "Lin_Algebra.h"
#include "Option_Pricer.h"
#include "gnuplot.h"

using namespace std;

void test_trig_mat()
{
	double m[] = { 1, 1 }, u[] = { 0 }, d[] = { 0 }, b[] = { 1, 2 };
	vector<double> mm(m, m + 2), uu(u, u + 1), dd(d, d + 1), bb(b, b + 2);
	trigonal_matrix A(mm, uu, dd);
	vector<double> c = solve_equation_trigonal(A, bb);
	cout << c[0] << c[1] << endl;
}

double payoff(double(spot))
{
	if (spot < 3000)
	{
		return(3000 - spot);
	}
	else
	{
		return 0;
	}
}

double forward(double(time))
{
	return(3100 * exp(0.01*time));
}

double vol(double(spot), double(time))
{
	return 0.2;
}

double discount(double(time))
{
	return(exp(-0.01*time));
}

double lb(double(time))
{
	return(1500);
}

double lr(double(time))
{
	return(0);
}

double hb(double(time))
{
	return(3800);
}

double hr(double(time))
{
	return(0);
}

double ep(double(spot), double(time))
{
	if (spot < 3000)
	{
		return(3000 - spot);
	}
	else
	{
		return 0;
	}
}

double pp(vector<double> spots, vector<double> times)
{
	double spot = spots[spots.size() - 1];
	
	for (int i = 0; i < spots.size(); ++i)
	{
		if (spots[i]>3800)
			return 0;
		else if (spots[i] < 1500)
			return 0;
	}

	if (spot < 3000)
	{
		return(3000 - spot);
	}
	else
	{
		return 0;
	}
}

double autocall_func(double spot, double time)
{
	if (spot > 3600){
		if (time > 0.98 || (time > 0.48&&time < 0.5))
			return 3600;
		else
			return 0;
	}
	return 0;
}

bool ko_func(double spot)
{
	if (spot < 2000)
		return true;
	return false;
}

double autocall_sec_payoff_func(double spot)
{
	if (spot>3600)
		return 3600;
	else
		return spot;
}

double ko_payoff_func(double spot)
{
	if (spot > 3600)
		return 0;
	else if (spot < 2000)
		return 0;
	else
		return 3600 - spot;
}

void test_plot()
{
	Gnuplot plot;
	vector<double> x = { 1, 2, 3 }, y = { 1, 2, 3 };
	plot(x, y);
	system("pause");
}

void test()
{
	double price,mc_error;
	pde_pricer_auto_call pde_pm;
	pde_pm.forward_func = forward;
	pde_pm.vol_func = vol;
	pde_pm.payoff_func = autocall_sec_payoff_func;
	pde_pm.discount_func = discount;
	pde_pm.expire_time = 1;
	pde_pm.autocall_func = autocall_func;
	pde_pm.ko_func = ko_func;
	pde_pm.is_continuous_barrier = true;
	pde_pm.ko_payoff_func = ko_payoff_func;
	cout << pde_pm.get_price() << endl;
	
	vector<double> y1 = pde_pm.get_grid_slice();
	vector<double> y2 = pde_pm.ko_pricer_ptr->get_grid_slice();
	vector<double> x = pde_pm.get_spot_slice();
	
	vector<double> y(y1.size());
	for (int i = 0; i < y.size(); ++i)
		y[i] = y1[i] + y2[i];


	//pde_pricer pde_pm;
	//pde_pm.forward_func = forward;
	//pde_pm.vol_func = vol;
	//pde_pm.payoff_func = payoff;
	//pde_pm.discount_func = discount;
	//pde_pm.expire_time = 1;
	//price = pde_pm.get_price();
	//cout << price << endl;
	//
	//pde_pricer_with_barrier pde_pm1;
	//pde_pm1.forward_func = forward;
	//pde_pm1.vol_func = vol;
	//pde_pm1.payoff_func = payoff;
	//pde_pm1.discount_func = discount;
	//pde_pm1.expire_time = 1;
	//pde_pm1.lower_b_func = lb;
	//pde_pm1.upper_b_func = hb;
	//pde_pm1.lower_r_func = lr;
	//pde_pm1.upper_r_func = hr;
	//price = pde_pm1.get_price();
	//cout << price << endl;
	//
	//pde_pricer_american pde_pm2;
	//pde_pm2.forward_func = forward;
	//pde_pm2.vol_func = vol;
	//pde_pm2.payoff_func = payoff;
	//pde_pm2.discount_func = discount;
	//pde_pm2.expire_time = 1;
	//pde_pm2.early_payoff_func = ep;
	//price = pde_pm2.get_price();
	//cout << price << endl;
	//
	//monte_carlo_pricer mc_pm(40000);
	//mc_pm.forward_func = forward;
	//mc_pm.vol_func = vol;
	//mc_pm.path_payoff_func = pp;
	//mc_pm.discount_func = discount;
	//mc_pm.expire_time = 1;
	//mc_pm.simulate();
	//price = mc_pm.get_price();
	//mc_error = mc_pm.get_error();
	//cout << price << endl << mc_error << endl;
	
	Gnuplot plot;
	//vector<vector<double>> y(3);
	//vector<vector<double>> z(1);
	//y[0] = pde_pm.get_grid_slice();
	//y[1] = pde_pm1.get_grid_slice();
	//y[2] = pde_pm2.get_grid_slice();
	//vector<double> spots = pde_pm.get_spot_slice();
	//z[0].resize(spots.size());
	//
	//for (int i = 0; i < spots.size(); ++i)
	//{
	//	if (i % 10 == 0)
	//		z[0][i] = mc_pm.get_shift_price(spots[i]);
	//	else
	//		z[0][i] = -1;
	//}

	//plot(spots,y,z);
	plot(x, y);
	system("pause");
	
}
