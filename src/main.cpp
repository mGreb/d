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

// Дебажный макрос
#define SHOW_N(a) std::cout << #a << ": " << (a) << std::endl
#define SHOW(a) std::cout << #a << ": " << (a) << " "

using std::cout;
using std::endl;

const double rhs_coeff = 1.0;  // коэффициент при правой части

const int start_point = 0;  // левый край отрезка по x
const int end_point = 1;  // правый край отрезка по x

const int start_time = 0;  // начальное время t
const int end_time = 1;  // конечное время t

const double alpha = 0.9;  // задайте альфа на полуинтервале (0;1]
const double dx = 0.1;  // шаг по х
const double dt = 0.1 * pow(dx, 3.0 / alpha) / rhs_coeff;  // шаг по времени
const double s = 0.001;  // самая первая ступенька (тесная зависимость с alpha, по какой формуле?)

const int stepAmountX = int((end_point - start_point) / dx + 1);  // количество точек по пространству
const int stepAmountT = int((end_time - start_time) / dt + 1);  // полные временные слои

const int time_steps_to_count = 1000;  // число временных слоев необходимых для вычисления
//const int time_steps_to_count = stepAmountT;  // число временных слоев необходимых для вычисления

const int threadNum = 8;  // number of threads for openmp

const int sums_hyper = 20;  // кол-во слагаемых для гипергеометрической функции

const int output_step = 1;

const int number_of_explicit_steps = time_steps_to_count - 1;

// массивы для прогонки ------------------------------------------------------------------------------------------------
std::vector<double> A((size_t)stepAmountX);  // Под-под-диагональ
std::vector<double> B((size_t)stepAmountX);  // Под-диагональ
std::vector<double> C((size_t)stepAmountX);  // Диагональ
std::vector<double> D((size_t)stepAmountX);  // Над-диагональ
std::vector<double> E((size_t)stepAmountX);  // Над-над-диагональ
std::vector<double> F((size_t)stepAmountX);  // Здесь правая часть системы
std::vector<double> X((size_t)stepAmountX);  // Решение
// ---------------------------------------------------------------------------------------------------------------------

std::vector<double> func((size_t)stepAmountX);      // Массив значений функции начального условия
std::vector<double> funcDer1((size_t)stepAmountX);  // 1 производная функции начального условия
std::vector<double> funcDer2((size_t)stepAmountX);  // 2 производная функции начального условия
std::vector<double> funcDer3((size_t)stepAmountX);  // 3 производная функции начального условия

std::array<std::array<double, stepAmountX>, time_steps_to_count> V;  // Массив для значений сеточной функции

std::vector<double> x_steps((size_t)stepAmountX);          // Массив шагов по пространству
std::vector<double> t_steps((size_t)time_steps_to_count);  // Массив шагов по времени

std::array<double, stepAmountX> Uxs;  // Массив значений функции сдвига
std::array<double, stepAmountX> UxsDer3;  // Массив третьих производных значений функции сдвига

std::array<double, time_steps_to_count> grunvald_coeffs;  // массив коэффициентов Грюнвальда-Летникова
std::vector<double> hyper_func;  // Массив значений гипергеометрической функции 2F1(alpha, alpha, alpha + 1, s / t)

//записывает результаты в файл, с нормировкой или без неё
void toFile
	()
{
	{
		std::fstream str;
		str.open("out.txt", std::ios::out);
		for(size_t n = 0; n < t_steps.size(); n += output_step)
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

namespace useful_const
{
	const double g_alpha = gamma(alpha);
	const double g_1_alpha_rev = 1.0 / gamma(1.0 - alpha);
	const double sin_pi_alpha_pi = sin(M_PI * alpha) / M_PI;
	const double dx3 = dx * dx * dx;
	const double dt_m_alpha = pow(dt, - alpha);
}

double Q1
	(const int    it,
	 const size_t ix)
{
	return (1.0 / s - 1.0 / t_steps[it]) * useful_const::sin_pi_alpha_pi * pow(s / (t_steps[it] - s), alpha) * func[ix];
}

double Q2
	(const int    it,
	 const size_t ix)
{
	return funcDer1[ix] * useful_const::sin_pi_alpha_pi / alpha * pow(s / t_steps[it], alpha) * hyper_geometric(alpha, alpha, alpha + 1.0, s / t_steps[it]);
}

double Q3
	(const int    it,
	 const size_t ix)
{
	return funcDer1[ix] * useful_const::sin_pi_alpha_pi * pow(s / (t_steps[it] - s), alpha - 1.0) / (1.0 - alpha);
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
	return   pow(t - t_steps[j1], 1.0 - alpha) * ((alpha - 1.0) * t_steps[j1] - (alpha - 2.0) * t_steps[j] - t) / (alpha - 1.0) / (alpha - 2.0)
	       - pow(t - t_steps[j] , 1.0 - alpha) * ((alpha - 1.0) * t_steps[j]  - (alpha - 2.0) * t_steps[j] - t) / (alpha - 1.0) / (alpha - 2.0);
}

double fracIn
	(const int    n,
	 const size_t i)
{
	double result = 0.0;
	
	for(int j = 1; j <= n - 1; ++j)
	{
		result += V[j][i] * I1(t_steps[n], j, j+1) + (V[j+1][i] - V[j][i]) / dt * I2(t_steps[n], j, j+1);
	}
	
	return useful_const::g_1_alpha_rev * result;
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
	
	return useful_const::g_1_alpha_rev * result;
}

double fracD
	(const int    n,
	 const size_t i)
{
	double result = 0.0;
	for(int j = 1; j < n; ++j)
	{
		result += grunvald_coeffs[j] * V[n - j + 1][i];
	}
	
	return useful_const::dt_m_alpha * result;
}

int main()
{
	omp_set_num_threads(threadNum);
	
	// переменные для подсчета потраченного времени
	const double start = omp_get_wtime();
	
	cout.precision(5);
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
	
	for(size_t i = 0; i < t_steps.size(); ++i)
	{
		t_steps[i] = s + i * dt;
	}
	
	for (size_t i = 0; i < func.size(); ++i)
	{
		func[i]     = initial_condition(x_steps[i]);
		funcDer1[i] = initial_condition_first_derivative(x_steps[i]);
		funcDer2[i] = initial_condition_second_derivative(x_steps[i]);
		funcDer3[i] = initial_condition_third_derivative(x_steps[i]);
	}
	
	for (size_t i = 0; i < Uxs.size(); ++i)
	{
		Uxs[i]     = func[i]     * pow(s, alpha - 1.0) / useful_const::g_alpha;
		UxsDer3[i] = funcDer3[i] * pow(s, alpha - 1.0) / useful_const::g_alpha;
	}
	
	grunvald_coeffs[0] = 1.0;
	for (size_t i = 1; i < grunvald_coeffs.size(); ++i)
	{
		grunvald_coeffs[i] = grunvald_coeffs[i - 1] * (alpha - i + 1.0) / i;
	}
	
	for(size_t i = 0; i < grunvald_coeffs.size(); ++i)
	{
		grunvald_coeffs[i] *= pow(- 1.0, i);
	}
	
	hyper_func.push_back(0);
	for(size_t i = 1; i < time_steps_to_count; ++i)
	{
		hyper_func.push_back(hyper_geometric(alpha, alpha, alpha + 1.0, t_steps[1] / t_steps[i]));
	}
	// -------------------------------------------------------------------------------------------------------------------
	
	// Шагаю явной схемой ------------------------------------------------------------------------------------------------
	// it - номер слоя который вычисляем
	for(int it = 1; it < number_of_explicit_steps; ++it)
	{
		// отступаю с краев по две ячейки, так как производная третьего порядка
#pragma omp parallel for
		for(size_t ix = 1; ix < x_steps.size() - 2; ++ix)
		{
			const double Vder3 = (V[it][ix+2] - 3.0 * V[it][ix+1] + 3.0 * V[it][ix] - V[it][ix-1]) / useful_const::dx3;
			
			V[it+1][ix] = - Q1(it, ix) - fracD(it, ix)
			              - (Uxs[ix] + V[it][ix]) * (Q2(it, ix) + Q3(it, ix) + 0.5 / dx * (fracIn(it, ix + 1) - fracIn(it, ix - 1)))
			              + rhs_coeff * (UxsDer3[ix] + Vder3);
			V[it+1][ix] *= 1.0 / useful_const::dt_m_alpha;
			
//			SHOW(it); SHOW(ix); SHOW(V[it+1][ix]); SHOW(Vder3); SHOW(Q1(it, ix)); SHOW(fracD(it, ix)); SHOW(Q2(it, ix)); SHOW(Q3(it, ix)); SHOW(fracIn(it, ix + 1)); SHOW(fracIn(it, ix - 1));
//			getchar();
		}
		
		// условия нулевого градиента слева и справа
		V[it+1][0] = V[it+1][1];
		V[it+1][stepAmountX - 1] = V[it+1][stepAmountX - 3];
		V[it+1][stepAmountX - 2] = V[it+1][stepAmountX - 3];
	}
	// -------------------------------------------------------------------------------------------------------------------
	
	// Вычисление неявной схемой -----------------------------------------------------------------------------------------
//	for(int i = number_of_explicit_steps; i < time_steps_to_count; ++i)
//	{
//		
//		
//		P3D::P3D_thomas::thomas(A, B, C, D, E, F, X);
//	}
	// -------------------------------------------------------------------------------------------------------------------
	
	// Выполняю обратную подстановку -------------------------------------------------------------------------------------
	for (size_t l = 0; l < V.size(); ++l)
	{
		for (size_t p = 0; p < V[l].size(); ++p)
		{
//			V[l][p] = V[l][p] + Uxs[p];
		}
	}
	// -------------------------------------------------------------------------------------------------------------------
	
	toFile();
	
	const double end = omp_get_wtime();
	const double total_time = end - start;
	cout << "time = " << total_time << " seconds\n";
	
	return 0;
}
