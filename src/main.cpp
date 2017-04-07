// System includes
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <array>
#include <libalglib/ap.h>
#include <libalglib/specialfunctions.h>

// User includes
#include "P3D_thomas.hpp"

using std::cout;
using std::endl;

const int rhs_coeff = 1;  // коэффициент при правой части

const int start_point = 0;  // левый край отрезка по x
const int end_point = 1;  // правый край отрезка по x

const int start_time = 0;  // начальное время t
const int end_time = 1;  // конечное время t

int n = 0;  // номер слоя по времени

const double alpha = 0.5;  // задайте альфа на полуинтервале (0;1]
const double dx = 0.01;  // шаг по х
const double dt = 0.01; // * pow(dx, 3.0 / alpha) / rhs_coeff;  // шаг по времени
const double s = 0.01;  // самая первая ступенька (тесная зависимость с alpha, по какой формуле?)

const int stepAmountX = int((end_point - start_point) / dx + 1);  // количество точек по пространству
const int stepAmountT = int((end_time - start_time) / dt + 1);  // полные временные слои

const int time_steps_to_count = 1000 < stepAmountT ? 1000 : stepAmountT;  // число временных слоев необходимых для вычисления

const int threadNum = 1;  // number of threads for openmp

const int sums_hyper = 20;  // кол-во слагаемых для гипергеометрической функции

const double sin_pi_alpha_pi = sin(M_PI * alpha) / M_PI;

// массивы для прогонки ------------------------------------------------------------------------------------------------
std::vector<double> A((size_t)stepAmountX);  // Под-под-диагональ
std::vector<double> B((size_t)stepAmountX);  // Под-диагональ
std::vector<double> C((size_t)stepAmountX);  // Диагональ
std::vector<double> D((size_t)stepAmountX);  // Над-диагональ
std::vector<double> E((size_t)stepAmountX);  // Над-над-диагональ
std::vector<double> F((size_t)stepAmountX);  // Здесь правая часть системы
std::vector<double> X((size_t)stepAmountX);  // Решение
// ---------------------------------------------------------------------------------------------------------------------

std::array<double, stepAmountX> func;      // Массив значений функции начального условия
std::array<double, stepAmountX> funcDer1;  // 1 производная функции начального условия
std::array<double, stepAmountX> funcDer2;  // 2 производная функции начального условия
std::array<double, stepAmountX> funcDer3;  // 3 производная функции начального условия

std::array<std::array<double, stepAmountX>, time_steps_to_count> V;  // Массив для значений сеточной функции

std::array<double, stepAmountX> x_steps;          // Массив шагов по пространству
std::array<double, time_steps_to_count> t_steps;  // Массив шагов по времени

std::array<double, stepAmountX> Us;  // Массив значений функции сдвига

std::array<double, time_steps_to_count> grunvald_coeffs;  // массив коэффициентов Грюнвальда-Летникова
std::vector<double> hyper_func;  // Массив значений гипергеометрической функции 2F1(alpha, alpha, alpha + 1, s / t)

//записывает результаты в файл, с нормировкой или без неё
void toFile
	()
{
	{
		std::fstream str;
		str.open("out.txt", std::ios::out);
		//for(int n=0;n<t;n+=49)	//для вывода малого количества графиков (чтобы не захламлять графиками картинку)
		for(size_t n = 0; n < t_steps.size(); ++n)
		{
			str << "#" << n << std::endl;  // для gnuplot
			for(size_t i = 0; i < x_steps.size(); ++i)
			{
				str << t_steps[n] << "\t" << x_steps[i] << "\t" << V[n][i] << endl;
			}
			str << endl;
		}
		str.close();
	}
	
	{
		std::fstream str;
		str.open("grunvald.txt", std::ios::out);
		for(size_t n = 0; n < grunvald_coeffs.size(); ++n)
		{
			str << n << "\t" << grunvald_coeffs[n] << "\t" << fabs(grunvald_coeffs[n]) << endl;
		}
		str.close();
	}
	
	{
		std::fstream str;
		str.open("hyper.txt", std::ios::out);
		for(size_t n = 0; n < hyper_func.size(); ++n)
		{
			str << n << "\t" << hyper_func[n] << endl;
		}
		str.close();
	}
}

double hyper_geometric
	(const double a,
	 const double b,
	 const double c,
	 const double z)
{
	if (z == 0.0)
	{
		return 1.0;
	}
	
	double result = 1.0;
	
	for (int k = 1; k <= sums_hyper; ++k)
	{
		double prod = 1.0;
		for (int l = 0; l <= k - 1; ++l)
		{
			prod *= (a + l) * (b + l) / (1.0 + l) / (c + l);
		}
		result += prod * pow(z, k);
	}
	
	return result;
}

double gamma
	(const double input)
{
	return alglib::gammafunction(input);
}
const double g_alpha = gamma(alpha);
const double g_1_alpha_rev = 1.0 / gamma(1.0 - alpha);

//функция - начальное условие
double initial_condition
	(const double x)
{
	return x * x - x * x * x;
}

double initial_condition_first_derivative
	(const double x)
{
	return 2.0 * x - 3.0 * x * x;
}

double initial_condition_second_derivative
	(const double x)
{
	return 2.0 - 6.0 * x;
}

double initial_condition_third_derivative
	(const double x)
{
	(void)x;
	
	return - 6.0;
}

double Q1
	(const int it,
	 const int ix)
{
	return (1.0 / s - 1 / t_steps[it]) * sin_pi_alpha_pi * pow(s / (t_steps[it] - s), alpha) * func[ix];
}

double Q2
	(const int it,
	 const int ix)
{
	return funcDer1[ix] * sin_pi_alpha_pi / alpha * pow(s / t_steps[it], alpha) * hyper_func[it];
}

double Q3
	(const int it,
	 const double ix)
{
	return funcDer1[ix] * sin_pi_alpha_pi * pow(s/ (t_steps[it] - s), alpha - 1.0) / (1.0 - alpha);
}

double I1
	(const double t,
	 const int    j,
	 const int    j1)
{
	return   pow(t - t_steps[j1], 1.0 - alpha) / (alpha - 1.0)
	       - pow(t - t_steps[j] , 1.0 - alpha) / (alpha - 1.0);
}

double I2
	(const double t,
	 const int    j,
	 const int    j1)
{
	return   pow(t - t_steps[j1], 1.0 - alpha) * ((alpha - 1.0) * t_steps[j1] - (alpha - 2.0) * t_steps[j] - t) / (alpha - 1.0) / (alpha - 2)
	       - pow(t - t_steps[j] , 1.0 - alpha) * ((alpha - 1.0) * t_steps[j]  - (alpha - 2.0) * t_steps[j] - t) / (alpha - 1.0) / (alpha - 2);
}

double fracIn
	(const int    n,
	 const int    i)
{
	double result = 0.0;
	
	for(int j = 0; j <= n - 1; ++j)
	{
		result += V[j][i] * I1(t_steps[n], j, j+1) + (V[j+1][i] - V[j][i]) / dt * I2(t_steps[n], j, j+1);
	}
	
	return g_1_alpha_rev * result;
}

double fracIn1
	(const int n,
	 const int i)
{
	double result = 0.0;
	
	for(int j = 1; j <= n - 1; ++j)
	{
		result += V[j][i] * I1(t_steps[n+1], j, j+1) + (V[j+1][i] - V[j][i]) / dt * I2(t_steps[n+1], j, j+1);
	}
	
	return g_1_alpha_rev * result;
}

double fracD
	(const int n,
	 const int i)
{
	double result = 0.0;
	
	for(int j = 0; j <= n; ++j)
	{
		result += grunvald_coeffs[j] * V[n - j + 1][i];
	}
	
	return pow(dt, - alpha) * result;
}

int main()
{
	omp_set_num_threads(threadNum);
	
	// переменные для подсчета потраченного времени
	const double start = omp_get_wtime();
	
	cout.precision(15);
	cout << "Number of threads: " << threadNum << endl;
	cout << "dt: " << dt << endl;
	cout << "dx: " << dx << endl;
	cout << "Step amount over x: " << stepAmountX << endl;
	cout << "Step amount over t: " << stepAmountT << endl;
	cout << "How much to count: " << time_steps_to_count << endl;
	
	// Заполняю некоторые массивы ----------------------------------------------------------------------------------------
	for(size_t i = 0; i < x_steps.size(); ++i)
	{
		x_steps[i] = i * dx;
	}
	
	t_steps[0] = 0;
	t_steps[1] = s;
	for(size_t i = 2; i < t_steps.size(); ++i)
	{
		t_steps[i] = s + (i - 1) * dt;
	}
	
	for (size_t i = 0; i < func.size(); ++i)
	{
		func[i]     = initial_condition(x_steps[i]);
		funcDer1[i] = initial_condition_first_derivative(x_steps[i]);
		funcDer2[i] = initial_condition_second_derivative(x_steps[i]);
		funcDer3[i] = initial_condition_third_derivative(x_steps[i]);
	}
	
	for (size_t i = 0; i < Us.size(); ++i)
	{
		Us[i] = func[i] * pow(s, alpha - 1.0) / g_alpha;
	}
	
	grunvald_coeffs[0] = pow(-1.0, 0.0) * 1.0;
	grunvald_coeffs[1] = pow(-1.0, 1.0) * alpha;
	for (size_t i = 2; i < grunvald_coeffs.size(); i++)
	{
		grunvald_coeffs[i] = pow(-1.0, i) * grunvald_coeffs[i - 1] * (alpha - i + 1.0) / i;
	}
	
	hyper_func.push_back(0);
	for(size_t i = 1; i < time_steps_to_count; ++i)
	{
		hyper_func.push_back(hyper_geometric(alpha, alpha, alpha + 1.0, t_steps[1] / t_steps[i]));
	}
	// -------------------------------------------------------------------------------------------------------------------
	
	
	// Выполняю обратную подстановку -------------------------------------------------------------------------------------
	for (size_t l = 0; l < V.size(); ++l)
	{
		for (size_t p = 0; p < V[l].size(); ++p)
		{
			V[l][p] = V[l][p] + Us[p];
		}
	}
	// -------------------------------------------------------------------------------------------------------------------
	
	toFile();
	
	const double end = omp_get_wtime();
	const double total_time = end - start;
	cout << "time = " << total_time << " seconds\n";
	
	system("gnuplot ./gnuplot/testOMP.plt");
	system("gnuplot ./gnuplot/plot.plt");
	
	return 0;
}
