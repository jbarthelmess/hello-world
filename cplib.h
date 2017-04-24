#ifndef __cplib_h_
#define __cplib_h_ 
/*
	Computational Physics Library Class
	
	Will include:
	Numeric Integration (Trapezoid, Simpsons, Gaussian 5-point Quadrature)
	Numeric Differentiation (5-point Central Difference)
	
	Root-Finding(Bisection, Newton-Rhaphson)
	Extrema-Finding(Golden Ratio)
	
	Factorial 
	
	Error Approximation
*/
#define DBL_MAX 1000000000000
#include <vector>
#include <cmath>
#include <fstream>
#include <cassert>
#include <iostream>


class cplib
{
	public:
		// Numeric Integration
		double trap_int(double (*func)(double), int steps, double min, double max);
		double simp_int(double (*func)(double), int steps, double min, double max);
		double gauss_int3(double (*func)(double), int steps, double min, double max);
		double gauss_int5(double (*func)(double), int steps, double min, double max);
		
		// Numeric Derivatives
		double forward_deriv(double (*func)(double), double point, double h);
		double cent_deriv_3(double (*func)(double), double point, double h);
		double cent_deriv_5(double (*func)(double), double point, double h);
		
		// Root Finding 
		double bisect(double (*func)(double), double max, double min);
		double secant(double (*func)(double), double x_0, double x_i);
		double nr(double (*func)(double), double x_0);
		
		// Solving differential equations
		void simple_diff(double (*func)(double, double), double y_0, double t_0, double max, std::ofstream& out);
		void n_order_diff(int n, double* coef, double (*func)(double, double), double* init, double t_0, double step, double length, std::ofstream& out);
		void n_coupled_diff(int n, double* x, double t, double h, double length, void (*derivs)(double, double*, double*), std::ofstream& out);
		
		// Extrema-Finding
		//double gold(double (*func)(double), double max, double min);
};


double cplib::trap_int(double (*func)(double), int steps, double min, double max)
{
	double total = 0;
	double width = (max-min)/(steps);
	for(int i=0; i<steps; i++)
		total += (((*func)(i*(width) + min)) + ((*func)(((i+1)*width)+min)))*width/2;
	return total;
}


double cplib::simp_int(double (*func)(double), int steps, double min, double max)
{
	double total = 0;
	double width = (max-min)/(steps-1);
	std::vector<double> weight;
	weight.push_back(1.0);
	for(int i=1; i<steps-1; i++)
	{
		if(i%2 != 0)
			weight.push_back(4.0);
		else
			weight.push_back(2.0);
	}
	weight.push_back(1.0);
	assert(weight.size() == steps);
	for(int i=0; i<steps; i++)
	{
		total= total + ((((*func)((i*width)+min))*weight[i]/3)*width);
	}
	return total;
}

double cplib::gauss_int3(double (*func)(double), int steps, double min, double max)
{
	double c1, c2, c3, x1, x2, x3;
	c1 = 5.0/9;
	c2 = 8.0/9;
	c3 = 5.0/9;
	
	x1 = -sqrt(0.6);
	x2 = 0;
	x3 = sqrt(0.6);
	
	double total = 0.0;
	for(int i=0; i<steps; i++)
	{
		total += c1*((max-min)/(2*steps))*(*func)(((max-min)/(2*steps))*x1 + ((max-min)*(2*i +1)/(2*steps)) + min);
		total += c2*((max-min)/(2*steps))*(*func)(((max-min)/(2*steps))*x2 + ((max-min)*(2*i +1)/(2*steps)) + min);
		total += c3*((max-min)/(2*steps))*(*func)(((max-min)/(2*steps))*x3 + ((max-min)*(2*i +1)/(2*steps)) + min);
	}
	return total;
}

double cplib::gauss_int5(double (*func)(double), int steps, double min, double max)
{
	double c1, c2, c3, c4, c5, x1, x2, x3, x4, x5;
	c1 = 0.236926885;
	c2 = 0.478628670;
	c3 = 0.568888889;
	c4 = 0.478628670;
	c5 = 0.236926885;
	
	x1 = -0.906179846;
	x2 = -0.538469310;
	x3 = 0.0;
	x4 = 0.538469310;
	x5 = 0.906179846;
	
	double total = 0.0;
	
	for(int i=0; i<steps; i++)
	{
		total += c1*((max-min)/(2*steps))*(*func)(((max-min)/(2*steps))*x1 + ((max-min)*(2*i +1)/(2*steps)) + min);
		total += c2*((max-min)/(2*steps))*(*func)(((max-min)/(2*steps))*x2 + ((max-min)*(2*i +1)/(2*steps)) + min);
		total += c3*((max-min)/(2*steps))*(*func)(((max-min)/(2*steps))*x3 + ((max-min)*(2*i +1)/(2*steps)) + min);
		total += c4*((max-min)/(2*steps))*(*func)(((max-min)/(2*steps))*x4 + ((max-min)*(2*i +1)/(2*steps)) + min);
		total += c5*((max-min)/(2*steps))*(*func)(((max-min)/(2*steps))*x5 + ((max-min)*(2*i +1)/(2*steps)) + min);
	}
	
	return total;
}

double cplib::secant(double (*func)(double), double x_0, double x_i)
{
	double x_start = x_0;
	double y_i = (*func)(x_i);
	double y_0 = (*func)(x_0);
	double slope = (y_i - y_0)/(x_i - x_0);
	while(fabs(y_i) > 0.00001 && fabs(x_0 - x_start) < 100000)
	{
		double x_t = x_i - ((y_i)/slope);
		double y_t = (*func)(x_t);
		slope = (y_t - y_i)/(x_t - x_i);
		x_0 = x_i;
		x_i = x_t;
		y_0 = y_i;
		y_i = y_t;
	}
	if(x_0 - x_start >= 100000)
		return DBL_MAX;
	return x_i;
}

double cplib::nr(double (*func)(double), double x_0)
{
	double x_start = x_0;
	double y_i = (*func)(x_0);
	double slope = cent_deriv_5(func, x_0, 0.0000001);
	while(fabs(y_i) > 0.0000001 && fabs(x_0 - x_start) < 100000)
	{
		x_0 = -(y_i/slope) + x_0;
		slope = cent_deriv_5(func, x_0, 0.0000001);
		y_i = (*func)(x_0);
	}
	if(x_0 - x_start >= 100000)
		return DBL_MAX;
	return x_0;
}

double cplib::forward_deriv(double (*func)(double), double point, double h)
{
	return ((*func)(point+h) - (*func)(point))/h;
}

double cplib::cent_deriv_3(double (*func)(double), double point, double h)
{
	return ((*func)(point+(h/2)) - (*func)(point-(h/2)))/h;
}

double cplib::cent_deriv_5(double (*func)(double), double point, double h)
{
	return ((*func)(point - h) + 8*(*func)(point + (h/2)) - 8*(*func)(point-(h/2)) - (*func)(point + h))/(6*h);
}

void cplib::simple_diff(double (*func)(double, double), double y_0, double t_0, double max, std::ofstream& out)
{
	double k[4];
	out << t_0 << "\t" << y_0 << std::endl;
	double step = 0.01;
	
	for(double i = t_0; i< max; i+=step)
	{
		k[0] = (*func)(t_0, y_0);
		k[1] = (*func)(t_0 + (0.5)*(step), y_0 + step*(0.5)*(k[0]));
		k[2] = (*func)(t_0 + (0.5)*(step), y_0 + step*(0.5)*(k[1]));
		k[3] = (*func)(t_0 + step, y_0 + step*k[2]);
		
		y_0 = y_0 + (step/6)*(k[0] + 2*k[1] + 2*k[2] + k[3]);
		
		out << i+step << "\t" << y_0 << std::endl;
	}
}

void cplib::n_order_diff(int n, double* coef, double (*func)(double, double), double* init, double t_0, double step, double length, std::ofstream& out)
{
	/*
	This code solves any n-th order linear differential equation with constant coefficients
	
	n - order of the differential equation
	coef - pointer to an array of coefficients, coef[0] is the coefficient in front of the 0th derivative and so on, should be n+1 items
	func - function pointer to the differential equation, assuming the form, (linear combo of derivatives) = func(x, t)
	init - pointer to array of initial values for derivatives, assuming all values are at time t_0,
	t_0 - inital time
	length - amount of time beyond inital time to take derivative
	step - desired step size
	out - file to send data points to
	
	*/
	double **k;
	k = new double*[n];
	for(int i = 0; i< n; i++)
		k[i] = new double[4];

	out << t_0 << "\t" << init[0] << std::endl;
	
	for(double h = t_0; h < t_0+length; h+= step)
	{
		// k_1 calculations: first, we calculate the highest order variable, 
		k[n-1][0] = (*func)(init[0], h);
		for(int i = 0; i< n; i++)
		{
			k[n-1][0] = k[n-1][0] - coef[i]*init[i];
		}
		k[n-1][0] = k[n-1][0]/coef[n];
		// we find it by moving all the other terms to the right hand side of the equation 
		
		// the other terms are found by equating them to the current value of the previous derivative
		for(int i = 0; i < n-1; i++)
			k[i][0] = init[i+1];
		
		// k_2 calculations, done similar to k_1
		k[n-1][1] = (*func)(init[0] + step*(k[0][0]/2), h + (step/2));
		for(int i = 0; i< n; i++)
		{
			k[n-1][1] = k[n-1][1] - coef[i]*(init[i] + step*(k[i][0]/2));
		}
		k[n-1][1] = k[n-1][1]/coef[n];
		
		for(int i = 0; i < n-1; i++)
			k[i][1] = init[i+1] + step*(k[i+1][0]/2);
		
		
		// k_3 calculations, done similar to k_2
		k[n-1][2] = (*func)(init[0] + step*(k[0][1]/2), h + (step/2));
		for(int i = 0; i< n; i++)
		{
			k[n-1][2] = k[n-1][2] - coef[i]*(init[i] + step*(k[i][1]/2));
		}
		k[n-1][2] = k[n-1][2]/coef[n];
		
		for(int i = 0; i < n-1; i++)
			k[i][2] = init[i+1] + step*(k[i+1][1]/2);
		
		
		// k_4 calculations
		k[n-1][3] = (*func)(init[0] + step*(k[0][2]), h + step);
		for(int i = 0; i< n; i++)
		{
			k[n-1][3] = k[n-1][3] - coef[i]*(init[i] + step*k[i][2]);
		}
		k[n-1][3] = k[n-1][3]/coef[n];
		
		for(int i = 0; i < n-1; i++)
			k[i][3] = init[i+1] + step*(k[i+1][2]);
		
		
		for(int i =0; i< n; i++)
			init[i] = init[i] + (step/6)*(k[i][0] + 2*k[i][1] + 2*k[i][2] + k[i][3]);
		
		out << h + step << "\t" << init[0] << std::endl;
	}		
	for(int i = 0; i < n; i++)
		delete [] k[i];
	delete [] k;
}

// Using numerical recipes version of rk4
// More Adaptable for coupled differential equations. 

// Hope to generalize to handle higher order coupled differential equations 

/*
n - number of coupled equations to solve
x - initial values for the function
t - starting t value
h - desired step size
length - defines domain over which to integrate function [t, t+length]
derivs - computes derivatives based on current t, current x values, and stores value in derivative array
out - output stream to send data
*/
void cplib::n_coupled_diff(int n, double* x, double t, double h, double length, void (*derivs)(double, double*, double*), std::ofstream& out)
{
	double *k1, *k2, *k3; // There are 4 k's but only need to store the first three
	k1 = new double[n];	  // The fourth will be stored in dxdt after calling derivs the fourth time
	k2 = new double[n];
	k3 = new double[n];
	double* dxdt = new double[n];
	//for(int i =0; i<n; i++)
	//	dxdt[i] = 0;
	
	
	double* x_0 = new double[n]; // Copy of initial values, so that we can change them when passing to derivs 
	
	out << t;
	for(int i=0; i<n; i++)
		out << "\t" << x[i];
	out << std::endl;
	
	for(double j = t; j < t + length; j+= h)
	{
		for(int i = 0; i< n; i++) // copy x values
			x_0[i] = x[i];
		
		(*derivs)(j, x_0, dxdt);
		for (int i=0; i<n; i++) // First k1 step. Take derivative calculated earlier and set k1 equal to that
		{	
			k1[i] = dxdt[i];
			x_0[i] = x[i] + (h/2.0)*k1[i];
		}
		
		(*derivs)(j + (h/2.0), x_0, dxdt); // Second k2 step. Reevaluate derivative based on k1 calculated earlier. 
		for (int i=0;i<n;i++)
		{
			k2[i] = dxdt[i];
			x_0[i] = x[i] + (h/2.0)*k2[i];
		}
		
		(*derivs)(j+(h/2.0), x_0, dxdt); // Third (k3) step
		for (int i=0;i<n;i++) 
		{
			k3[i] = dxdt[i];
			x_0[i] = x[i] + h*k3[i];
		}
		
		(*derivs)(j+h, x_0, dxdt);		// Fourth (k4) step
		for (int i=0; i<n; i++)
			x[i]=x[i]+ h*(dxdt[i]+k1[i]+2.0*(k2[i] + k3[i]))/6.0;
		
		
		out << j;
		for(int i = 0; i< n; i++)
			out << "\t" << x[i];
		out << std::endl;
		/*
		if(int(j)%100 == 0)
		{
			std::cout << "energy of Earth: " << 0.5*0.00000300*(x[1]*x[1] + x[3]*x[3]) - 0.00000300/pow(x[0]*x[0] + x[2]*x[2], 0.5) << std::endl;;
		}
		*/
	}
	delete [] k1;
	delete [] k2;
	delete [] k3;
	delete [] x_0;
	delete [] dxdt;
}

double cplib::bisect(double (*func)(double), double min, double max)
{
	if((*func)(max) == 0)
		return max;
	if((*func)(min) == 0)
		return min;
	
	int loop_count = 0;
	double test = (*func)((max + min)/2);
	if((*func)(max) > 0 && (*func)(min) < 0)
	{
		while(fabs(test) > 0.00001 && loop_count < 1000)
		{
			if(test > 0)
				max = (max+min)/2;
			if(test < 0)
				min = (max+min)/2;
			test = (*func)((max+min)/2);
			loop_count++;
		}
		if(loop_count >= 1000)
			return max+1;
		return (max+min)/2;
	}
	else if ((*func)(max) > 0 && (*func)(min) > 0)
	{
		while(fabs(test) > 0.00001 && loop_count < 1000)
		{
			if(test < 0)
				return bisect(func, (max+min)/2, max);
			if(cent_deriv_5(func, (max+min)/2, 0.0000001) > 0)
				max = (max+min)/2;
			else if(cent_deriv_5(func, (max+min)/2, 0.0000001) < 0)
				min = (max+min)/2;
			test = (*func)((max+min)/2);
			loop_count++;
		}
		if(loop_count >= 1000)
			return max+1;
		return (max+min)/2;
	}
	else if ((*func)(max) <0 && (*func)(min) < 0)
	{
		while(fabs(test) > 0.00001 && loop_count < 1000)
		{
			if(test > 0)
				return bisect(func, (max+min)/2, max);
			if(cent_deriv_5(func, (max+min)/2, 0.0000001) > 0)
				min = (max+min)/2;
			else if(cent_deriv_5(func, (max+min)/2, 0.0000001) < 0)
				max = (max+min)/2;
			test = (*func)((max+min)/2);
			loop_count++;
		}
		if(loop_count >= 1000)
			return max+1;
		return (max+min)/2;
	}
	else
	{
		while(fabs(test) > 0.00001 && loop_count < 1000)
		{
			if(test > 0)
				min = (max+min)/2;
			if(test < 0)
				max = (max+min)/2;
			test = (*func)((max+min)/2);
			loop_count++;
		}
		if(loop_count >= 1000)
			return max+1;
		return (max+min/2);
	}
}

#endif

