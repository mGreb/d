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
const double dx = 0.1;  // шаг по х
const double dt = 0.1 * pow(dx, 3.0 / alpha) / rhs_coeff;  // шаг по времени
const double s = 0.01;  // самая первая ступенька (тесная зависимость с alpha, по какой формуле?)

const int stepAmountX = int((end_point - start_point) / dx + 1);  // количество точек по пространству
const int stepAmountT = int((end_time - start_time) / dt);  // полные временные слои

const int time_steps_to_count = 1000;  // число временных слоев необходимых для вычисления

const int threadNum = 1;  // number of threads for openmp

const int sums_hyper = 20;  // кол-во слагаемых для гипергеометрической функции

// массивы для прогонки ------------------------------------------------------------------------------------------------
std::vector<double> A((size_t)stepAmountX);  // под-под-диагональ
std::vector<double> B((size_t)stepAmountX);  // под-диагональ
std::vector<double> C((size_t)stepAmountX);  // диагональ
std::vector<double> D((size_t)stepAmountX);  // над-диагональ
std::vector<double> E((size_t)stepAmountX);  // над-над-диагональ
std::vector<double> F((size_t)stepAmountX);  // здесь правая часть системы
std::vector<double> X((size_t)stepAmountX);  // решение
// ---------------------------------------------------------------------------------------------------------------------

std::array<double, stepAmountX> func;      // массив значений функции начального условия
std::array<double, stepAmountX> funcDer1;  // 1 производная функции начального условия
std::array<double, stepAmountX> funcDer2;  // 2 производная функции начального условия
std::array<double, stepAmountX> funcDer3;  // 3 производная функции начального условия

std::array<std::array<double, stepAmountX>, time_steps_to_count> V;  // Массив для значений сеточной функции

std::array<double, stepAmountX> x_steps;          // Массив шагов по пространству
std::array<double, time_steps_to_count> t_steps;  // Массив шагов по времени

std::array<double, stepAmountX> Us;  // Массив значений функции сдвига

std::array<double, time_steps_to_count> grunvald_coeffs;  // массив коэффициентов грюнвальда-летникова

//записывает результаты в файл, с нормировкой или без неё
void toFile
	(const int    t,
	 const double d)
{
	std::fstream str;
	str.open("out.txt", std::ios::out);
	
	//for(int n=0;n<t;n+=49)	//для вывода малого количества графиков (чтобы не захламлять графиками картинку)
	for (int n = 0; n < t; n++)
	{
		str << "#" << n << std::endl;  // для gnuplot
		double shag = 0;
		for (int i = 0; i < stepAmountX; i++)
		{
			str << n * dt << "\t" << shag << "\t" << V[n][i] << endl;
			shag = shag + d;
		}
		str << endl;
	}
	str.close();
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

//функция - граничное условие
double g
	(const double t)
{
	(void)t;
	
	return 0.;
	//return gamma(alpha);
	//return 1.;
}

// dt/()
//вычисляет интеграл
//podstanovka - Tn
//up - T(j+1)
//down - Tj
double newIntegral1
	(const double podstanovka,
	 const int    check)
{
	double result = 0.0;
	
	if (podstanovka > check)
	{
		result = (1.0 / (alpha - 1.0))
		         * (  pow(s + podstanovka * dt -  check        * dt, 1.0 - alpha)
		            - pow(s + podstanovka * dt - (check - 1.0) * dt, 1.0 - alpha));
	}
	
	return result;
}

// t*dt/()
// podstanovka - Tn
// up   - T(j+1)
// down - Tj
double newIntegral2
	(const double podstanovka,
	 const int    check)
{
	double result = 0.0;
	
	if (podstanovka > check)
	{
		result = (1.0 / ((alpha - 2.0) * (alpha - 1.0)))
		         * (  pow(s + podstanovka * dt - check        * dt, 1.0 - alpha) * (alpha *  check        * dt - s - podstanovka * dt -  check       * dt)
		            - pow(s + podstanovka * dt -(check - 1.0) * dt, 1.0 - alpha) * (alpha * (check - 1.0) * dt - s - podstanovka * dt - (check - 1.0)* dt));
	}
	
	return result;
}

//подсчет дробного интеграла на n слое
//down - номер временного слоя
//order -  порядок дробного интеграла
//step - номер шага по Ox
double fracInt
	(const double down,
	 const double order,
	 const int    step)
{
	//случай выхода за границу сетки
	if (step >= stepAmountX)
	{
		return 0.;
	}
	
	double result = 0.0;
	
	for (int j = 1; j < down; j++)
	{
		result += newIntegral1(n - 1.0, j + 1)
		          * (V[j][step] + ((s + (j - 1) * dt) / dt) * (- V[j+1][step] + V[j][step]))
		          + 1.0 / dt * newIntegral2(n - 1.0, j + 1) * (  V[j+1][step] - V[j][step]);
	}
	
	return result / gamma(order);
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
	for (size_t i = 0; i < func.size(); ++i)
	{
		const double point = i * dx;
		func[i]     = initial_condition(point);
		funcDer1[i] = initial_condition_first_derivative(point);
		funcDer2[i] = initial_condition_second_derivative(point);
		funcDer3[i] = initial_condition_third_derivative(point);
	}
	
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
	
	for (int i = 0; i < Us.size(); ++i)
	{
		Us[i] = func[i] * pow(s, alpha - 1.0) / g_alpha;
	}
	
	grunvald_coeffs[0] = pow(-1.0, 0.0) * 1.0;
	grunvald_coeffs[1] = pow(-1.0, 1.0) * alpha;
	for (int i = 2; i < grunvald_coeffs.size(); i++)
	{
		grunvald_coeffs[i] = pow(-1.0, i) * grunvald_coeffs[i - 1] * (alpha - i + 1.0) / i;
	}
	// -------------------------------------------------------------------------------------------------------------------
	
	double sum1 = 0.0;
	double sum2 = 0.0;
	double sum3 = 0.0;
	double sum4 = 0.0;
	double sum5 = 0.0;
	double sum6 = 0.0;
	double sumGroup = 0.0;
	double sumGroup2 = 0.0;
	
	// Явная схема
	int i = 0;
	
	const double sin_m_pi = sin(M_PI * alpha);
	
	for (n = 1; n < time_steps_to_count; n++)
	{
		for (i = 2; i < stepAmountX - 1; i++)
		{
			sum2 =   sin_m_pi * pow(s, alpha)       * func[i] / (M_PI * pow(n * dt, alpha) * (s + n * dt))
			       - sin_m_pi * pow(s, alpha - 1.0) * func[i] /  M_PI * pow(n * dt, alpha);
			sum3 = V[n][i] + func[i] * pow(s,alpha - 1.0) / g_alpha;
			sum4 = sin_m_pi / (M_PI * alpha) * funcDer1[i] * pow(s, alpha) * hyper_geometric(alpha, alpha, alpha + 1.0, s / (s + (n - 1) * dt))
			       / pow(s + (n - 1) * dt, alpha) + funcDer1[i] * pow(s, alpha - 1.0) / g_alpha;
			sum5 = 1.0 / (2.0 * dx) * (fracInt(n, 1.0 - alpha, i + 1) - fracInt(n, 1.0 - alpha, i - 1));
			sum6 =   rhs_coeff * (V[n][i - 2] - 3.0 * V[n][i - 1] + 3.0 * V[n][i] - V[n][i + 1]) / (dx * dx * dx)
			       + rhs_coeff * pow(s,alpha - 1.0) * funcDer3[i] / g_alpha;
			
			V[n+1][i] = - sum1 + pow(dt, alpha) * (sum2 - sum3 * (sum4 + sum5) + sum6);
		}
		
		V[n][0] = V[n][1] + Us[0] + Us[1];
		
		V[n][stepAmountX - 2] = V[n][stepAmountX - 3] + Us[stepAmountX - 3] - Us[stepAmountX - 2];
		V[n][stepAmountX - 1] = V[n][stepAmountX - 2] + Us[stepAmountX - 2] - Us[stepAmountX - 1];
	}
	
	int m = time_steps_to_count-1;
	
	V[m][0] = V[m][1] + Us[0] + Us[1];
	
	V[m][stepAmountX - 2] = V[m][stepAmountX - 3] + Us[stepAmountX - 3] - Us[stepAmountX - 2];
	V[m][stepAmountX - 1] = V[m][stepAmountX - 2] + Us[stepAmountX - 2] - Us[stepAmountX - 1];
	
	// выполняю обратную подстановку
	for (int l = 0; l < time_steps_to_count; ++l)
	{
		for (size_t p = 0; p < stepAmountX; ++p)
		{
			V[l][p] = V[l][p] + Us[p];
		}
	}
	
	toFile(time_steps_to_count,dx);
	
	const double end = omp_get_wtime();
	const double total_time = end - start;
	cout << "time = " << total_time << " seconds\n";
	
	system("gnuplot testOMP.plt");
	system("gnuplot plot.plt");
	
	return 0;
}
