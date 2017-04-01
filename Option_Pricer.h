#include "Lin_Algebra.h"
class pricer
{
public:
	virtual double get_price();
	double(*forward_func)(double);
	double(*vol_func)(double, double);
	double(*payoff_func)(double);
	double(*discount_func)(double);
	double expire_time;
};

class pde_pricer : public pricer
{
	/*
	S(t) = F(t)exp(x(t))
	dx = -0.5\sigma(x,t)^2 dt + sigma(x,t) dw
	g(x,t) = E[P(X,T)|(x,t)]
	dg/dt + (-0.5\sigma^2)dg/dx + 0.5 \sigma^2 d^2g/dx^2 = 0
	price(x,t) = discount(t,T)g(x,t)
	*/

public:
	pde_pricer();
	virtual double get_price();
	virtual void mush_grid();
	virtual vector<double> get_grid_slice(){ return grid_slice; };
	virtual vector<double> get_spot_slice(){ return spot_slice; };
	virtual double lower_price_func();
	virtual double upper_price_func();
protected:
	trigonal_matrix trig_mat;
	double slice_time, time_step, spot_step, spot_range_min, spot_range_max, lower_x, upper_x;
	vector<double> grid_slice;
	vector<double> spot_slice;

	virtual void set_payoff_slice();
	virtual void step_back_one_dt();
	void set_matrix();
	virtual void set_b();
};

class pde_pricer_with_barrier : public pde_pricer
{
public:
	virtual void mush_grid();
	double(*lower_b_func)(double);
	double(*upper_b_func)(double);
	double(*lower_r_func)(double);
	double(*upper_r_func)(double);
	virtual double lower_price_func();
	virtual double upper_price_func();
};

class pde_pricer_american : public pde_pricer
{
public:
	virtual void mush_grid();
	double(*early_payoff_func)(double, double);
};

class pde_pricer_american_barrier : virtual public pde_pricer_american, virtual public pde_pricer_with_barrier
{
public:
	virtual void mush_grid();
};

class pde_pricer_auto_call : public pde_pricer
{
public:
	pde_pricer_auto_call(bool is_main_pricer=true);
	~pde_pricer_auto_call();
	virtual double get_price();
	virtual void mush_grid();

	double(*autocall_func)(double, double);
	bool(*ko_func)(double);
	double(*ko_payoff_func)(double);
	bool is_continuous_barrier;
	virtual double lower_price_func();
	virtual double upper_price_func();
	pde_pricer_auto_call* ko_pricer_ptr;
	bool is_main_pricer;
};

class monte_carlo_pricer : public pricer
{
	/*
	S(t) = F(t)exp(x(t))
	dx = -0.5\sigma(x,t)^2 dt + sigma(x,t) dw
	g(x,t) = E[P(X,T)|(x,t)]
	dg/dt + (-0.5\sigma^2)dg/dx + 0.5 \sigma^2 d^2g/dx^2 = 0
	price(x,t) = discount(t,T)g(x,t)
	*/
public:
	monte_carlo_pricer( int num_paths = 40000 );
	virtual double get_price();
	virtual double get_error();
	double get_shift_price(double log_shift);// quick and dirty - assume no vol skew - just scale all the spots path up.
	void simulate();
	double(*path_payoff_func)(vector<double>,vector<double>);
protected:
	vector<vector<double>> random_numbers; // on w.
	vector<vector<double>> paths; // on S.
	vector<double> event_times;
	vector<double> times;
	void set_times();
private:
	int number_paths;
	void gen_random_numbers();
	void gen_paths();
};

