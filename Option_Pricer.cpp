// Option_Pricer.cpp : Defines the entry point for the console application.
#include "stdafx.h"
#include<iostream>
#include<vector>
#include "Lin_Algebra.h"
#include "Unit_Test.h"
#include "Option_Pricer.h"
#include<random>
using namespace std;

//////////////////////////////////////////Base Pricer/////////////////////////////////////////

double pricer::get_price()
{
	return 0;
}

/////////////////////////////////////////pde pricer///////////////////////////////////////////

pde_pricer::pde_pricer()
{
	time_step = 1 / 365.2425;
	spot_range_min = -5; // 5 vol
	spot_range_max = 5; // 5 vol
	spot_step = 0.01;
};
	
double pde_pricer::get_price()
{
	set_payoff_slice();
	while (slice_time > time_step)
			step_back_one_dt();
	//step the last non-standard step.
	time_step = slice_time;
	step_back_one_dt();

	int lower_index = floor((-lower_x) / spot_step);
	double px_lower = grid_slice[lower_index], px_higher = grid_slice[lower_index + 1];
	double spot_lower = spot_slice[lower_index], spot_higher = spot_slice[lower_index + 1];
	double pct_higher = -spot_lower / (spot_higher - spot_lower);
	double px = px_lower * (1-pct_higher) + px_higher * pct_higher;
	px *= discount_func(expire_time);
	return px;
}

void pde_pricer::mush_grid()
{}

void pde_pricer::set_payoff_slice()
{
	double fwd = forward_func(expire_time);
	double ref_vol = vol_func(fwd, expire_time);
	double grid_vol;
	if (ref_vol*sqrt(expire_time) > spot_step)
		grid_vol = ref_vol*sqrt(expire_time);
	else
		grid_vol = spot_step;

	upper_x = spot_range_max*grid_vol;
	lower_x = spot_range_min*grid_vol;
	int N = int((upper_x - lower_x) / spot_step);
	grid_slice.resize(N);
	spot_slice.resize(N);
	double log_spot;
	for (int i = 0; i < N; ++i)
	{
		log_spot = lower_x + i * spot_step;
		spot_slice[i] = log_spot;
		grid_slice[i] = payoff_func(fwd*exp(log_spot));
	}
	slice_time = expire_time;
	upper_x = log_spot; // re-line up
	mush_grid();
}

void pde_pricer::set_matrix()
{
	double new_slice_time = slice_time - time_step;
	double vol, c1, c2;
	vector<double> diag(grid_slice.size()), upper(grid_slice.size() - 1), lower(grid_slice.size() - 1);

	// diag is -g_ii - dt/dx^2 sigma_ii^2
	// upper is -dt/(4dx)sigma(i-1,i)^2+dt/(2dx^2)sigma(i-1,i)^2
	// lower is dt/(4dx)sigma(i-1,i)^2+dt/(2dx^2)sigma(i-1,i)^2
	for (int i = 1; i < diag.size() - 1; ++i)
	{
		vol = vol_func(forward_func(new_slice_time)*exp(spot_slice[i]), new_slice_time);
		diag[i] = -1 - time_step / (spot_step*spot_step) * vol * vol;
		c1 = time_step / (4 * spot_step)*vol*vol;
		c2 = time_step / (2 * spot_step*spot_step)*vol*vol;

		upper[i] = -c1 + c2;
		lower[i - 1] = c1 + c2;
	}

	//lower boundary and upper boundary
	// g(x,t) = E[P(X,T)] = E[Payoff(exp(X),T)] where we assume Payoff(y,T) is approximately linear in y
	// so we approx the boundary condition by
	// g(x,t)~ Payoff( E(exp(X)), T ) ~ Payoff( exp(x), T ) as exp(x) is martingale.
	double fwd = forward_func(new_slice_time);
	double lower_vol, upper_vol, lower_term, upper_term;
	lower_vol = vol_func(forward_func(new_slice_time)*exp(lower_x), new_slice_time);
	upper_vol = vol_func(forward_func(new_slice_time)*exp(upper_x), new_slice_time);
	lower_term = time_step / (4 * spot_step)*lower_vol*lower_vol + time_step / (2 * spot_step*spot_step)*lower_vol*lower_vol;
	upper_term = -time_step / (4 * spot_step)*upper_vol*upper_vol + time_step / (2 * spot_step*spot_step)*upper_vol*upper_vol;
	diag[0] = -1 - time_step / (spot_step*spot_step) * lower_vol * lower_vol;
	diag[diag.size() - 1] = -1 - time_step / (spot_step*spot_step) * upper_vol * upper_vol;
	upper[0] = -time_step / (4 * spot_step)*lower_vol*lower_vol + time_step / (2 * spot_step*spot_step)*lower_vol*lower_vol;
	lower[diag.size() - 2] = time_step / (4 * spot_step)*upper_vol*upper_vol + time_step / (2 * spot_step*spot_step)*upper_vol*upper_vol;

	trig_mat.set(diag, upper, lower);
}

double pde_pricer::lower_price_func()
{
	double fwd = forward_func(slice_time);
	double lower_b = payoff_func(fwd*exp(lower_x));
	return lower_b;
}

double pde_pricer::upper_price_func()
{
	double fwd = forward_func(slice_time);
	double upper_b = payoff_func(fwd*exp(upper_x));
	return upper_b;
}

void pde_pricer::set_b()
{
	double fwd = forward_func(slice_time);
	double lower_vol, upper_vol, lower_term, upper_term;
	double lower_b = lower_price_func(), upper_b = upper_price_func();

	lower_vol = vol_func(forward_func(slice_time)*exp(lower_x), slice_time);
	upper_vol = vol_func(forward_func(slice_time)*exp(upper_x), slice_time);
	lower_term = time_step / (4 * spot_step)*lower_vol*lower_vol + time_step / (2 * spot_step*spot_step)*lower_vol*lower_vol;
	upper_term = -time_step / (4 * spot_step)*upper_vol*upper_vol + time_step / (2 * spot_step*spot_step)*upper_vol*upper_vol;

	grid_slice[0] += lower_term *lower_b;
	grid_slice[grid_slice.size() - 1] += upper_term *upper_b;

	for (int i = 0; i < grid_slice.size(); ++i)
		grid_slice[i] *= -1;
}

void pde_pricer::step_back_one_dt()
{
	set_matrix();
	slice_time -= time_step;
	set_b();
	grid_slice = solve_equation_trigonal(this->trig_mat, grid_slice);
	mush_grid();
}

///////////////////////////////////////////////////////////////barrier pricer///////////////////////////////////////////////

void pde_pricer_with_barrier::mush_grid()
{
	double lower_b = lower_b_func(slice_time), upper_b = upper_b_func(slice_time);
	double fwd = forward_func(slice_time), spot;
	for (int i = 0; i < spot_slice.size(); ++i)
	{
		spot = fwd*exp(spot_slice[i]);
		if ( spot < lower_b)
		{
			grid_slice[i] = lower_r_func(slice_time);
		}
		else if (spot > upper_b)
		{
			grid_slice[i] = upper_r_func(slice_time);
		}
	}
	double testtest = 1;
}

double pde_pricer_with_barrier::lower_price_func()
{
	double fwd = forward_func(slice_time);
	double lower_spot = fwd*exp(lower_x);
	double lower_b = payoff_func(lower_spot);

	if (lower_spot < lower_b_func(slice_time))
		lower_b = lower_r_func(slice_time);

	return lower_b;
}

double pde_pricer_with_barrier::upper_price_func()
{
	double fwd = forward_func(slice_time);
	double upper_spot = fwd*exp(upper_x);
	double upper_b = payoff_func(upper_spot);

	if (upper_spot > upper_b_func(slice_time))
		upper_b = upper_r_func(slice_time);
	return upper_b;
}

/////////////////////////////////////////////////////////////////////////American Pricer///////////////////////////////////////////

void pde_pricer_american::mush_grid()
{
	double fwd = forward_func(slice_time);
	double this_spot,this_payoff,expected_payoff;
	for (int i = 0; i < spot_slice.size(); ++i)
	{
		this_spot = fwd*exp(spot_slice[i]);
		this_payoff = early_payoff_func(this_spot, slice_time);
		expected_payoff = this_payoff * (discount_func(slice_time) / discount_func(expire_time));
		if (grid_slice[i] < expected_payoff)
			grid_slice[i] = expected_payoff;
	}
}

///////////////////////////////////////////////////////////////American Barrier Pricer///////////////////////////////////////////////
void pde_pricer_american_barrier::mush_grid()
{
	pde_pricer_with_barrier::mush_grid();
	pde_pricer_american::mush_grid();
}

///////////////////////////////////////////////////////////////Monte Carlo Pricer///////////////////////////////////////////////////

monte_carlo_pricer::monte_carlo_pricer(int num_paths): number_paths(num_paths)
{}

void monte_carlo_pricer::simulate()
{
	set_times();
	gen_random_numbers();
	gen_paths();
}

double monte_carlo_pricer::get_price()
{
	double price = 0;
	for (int i = 0; i < paths.size(); ++i)
		price = (price * i + path_payoff_func(paths[i], times)) / (i + 1);
	return price *discount_func(expire_time);
}

double monte_carlo_pricer::get_shift_price(double log_shift)
{
	double price = 0;
	for (int i = 0; i < paths.size(); ++i)
		price = (price * i + path_payoff_func(paths[i]*exp(log_shift), times)) / (i + 1);
	return price *discount_func(expire_time);
}

double monte_carlo_pricer::get_error()
{
	double price = get_price();
	double std_dev = 0;
	for (int i = 0; i < paths.size(); ++i)
		std_dev = (std_dev * i + (path_payoff_func(paths[i],times) - price) * (path_payoff_func(paths[i],times) - price)) / (i + 1);
	return sqrt(std_dev)/sqrt(number_paths); // the Gaussian convergence theorem.
}

void monte_carlo_pricer::set_times()
{
	double min_time_step = 1 / 365.2425; // less than 3m
	double max_time_step = 7 / 365.2425; // larger than 3m

	times.clear();
	times.push_back(0);
	double time_step = min_time_step;
	double current_time = 0;
	int event_index = 0;
	while (current_time < expire_time)
	{
		if (current_time > 0.25)
			time_step = max_time_step;
		if (event_index >= event_times.size())
		{
			if (expire_time - current_time < time_step)
				current_time = expire_time;
			else
				current_time += time_step;
		}
		else if (event_times[event_index] - current_time < time_step)
		{
			current_time = event_times[event_index];
			++event_index;
		}
		else
		{
			current_time += time_step;
		}
		times.push_back(current_time);
	}
}

void monte_carlo_pricer::gen_random_numbers()
{
	default_random_engine generator;
	normal_distribution<double> distribution(0.0, 1.0);

	random_numbers.resize(number_paths);
	for (int i = 0; i < number_paths; ++i)
	{
		random_numbers[i].resize(times.size());
		random_numbers[i][0] = 0;
		for (int j = 1; j < times.size(); ++j)
			random_numbers[i][j] = distribution(generator) * sqrt(times[j] - times[j - 1]);
	}
}

void monte_carlo_pricer::gen_paths()
{
	paths.resize(number_paths);
	double x,vol;
	for (int i = 0; i < number_paths; ++i)
	{
		paths[i].resize(times.size());
		paths[i][0] = forward_func(times[0]);
		x = 0;
		for (int j = 1; j < times.size(); ++j)
		{
			vol = vol_func(paths[i][j - 1], times[j - 1]);
			x += -0.5 * vol * vol * (times[j] - times[j - 1]) + vol * random_numbers[i][j];
			paths[i][j] = forward_func(times[j]) * exp(x);
		}
	}
}

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\Auto Pricer in Grid\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

pde_pricer_auto_call::pde_pricer_auto_call(bool is_main_pricer) :is_main_pricer(is_main_pricer)
{
	if (is_main_pricer)
		ko_pricer_ptr = new pde_pricer_auto_call(false);
	else
		ko_pricer_ptr = nullptr;
}

pde_pricer_auto_call::~pde_pricer_auto_call()
{
	if (is_main_pricer)
		delete ko_pricer_ptr;
}

double pde_pricer_auto_call::lower_price_func()
{
	return 0;
}

double pde_pricer_auto_call::upper_price_func()
{
	if (!is_main_pricer)
		return 0;
	else
		return payoff_func(forward_func(expire_time)*exp(upper_x));
}

double pde_pricer_auto_call::get_price()
{
	if (is_main_pricer)
	{
		ko_pricer_ptr->payoff_func = ko_payoff_func;
		ko_pricer_ptr->forward_func = forward_func;
		ko_pricer_ptr->discount_func = discount_func;
		ko_pricer_ptr->vol_func = vol_func;
		ko_pricer_ptr->autocall_func = autocall_func;
		ko_pricer_ptr->expire_time = expire_time;
		ko_pricer_ptr->is_continuous_barrier = is_continuous_barrier;
		ko_pricer_ptr->ko_func = ko_func;
	}

	double main_price = pde_pricer::get_price();
	double ko_price = ko_pricer_ptr->pde_pricer::get_price();

	return main_price + ko_price;
}

void pde_pricer_auto_call::mush_grid()
{
	double fwd = forward_func(slice_time);

	if (!is_main_pricer)
	{
		for (int i = 0; i < spot_slice.size(); ++i)
		{
			if (autocall_func(fwd*exp(spot_slice[i]), slice_time)!=0.0)
				grid_slice[i] = 0;
		
			if (is_continuous_barrier && ko_func(fwd*exp(spot_slice[i])))
				grid_slice[i] = 0;
		}
	}
	else
	{
		for (int i = 0; i < spot_slice.size(); ++i)
		{
			if (autocall_func(fwd*exp(spot_slice[i]), slice_time) != 0.0)
				grid_slice[i] = autocall_func(fwd*exp(spot_slice[i]), slice_time)*discount_func(slice_time)/discount_func(expire_time);
		}
	}
}
