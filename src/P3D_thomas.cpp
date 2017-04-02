//
// Created by mgreb on 11/3/16.
//

#include "P3D_thomas.hpp"

// Системные файлы
#include <cstdio>
#include <cmath>
#include <cstdlib>

void P3D::P3D_thomas::thomas
	(const std::vector<double> &top,
	 const std::vector<double> &middle,
	 const std::vector<double> &bottom,
	 const std::vector<double> &right,
	       std::vector<double> &result)
{
	const size_t size = result.size();
	
#ifndef NDEBUG
	for(size_t i = 0; i < size; ++i)
	{
		if((fabs(top[i]) + fabs(bottom[i])) > fabs(middle[i]))
		{
			printf("System matrix is not diagonally dominant!\n");
			exit(1);
		}
	}
#endif
	
	std::vector<double> c(size);
	std::vector<double> d(size);
	
	c[0] = top[0] / middle[0];
	d[0] = right[0] / middle[0];
	for (size_t i = 1; i < size; ++i)
	{
		double denominator = (middle[i] - bottom[i] * c[i-1]);
		c[i] = top[i] / denominator;
		d[i] = (right[i] - bottom[i] * d[i-1]) / denominator;
	}
	
	result.back() = d.back();
	for (int i = (int)size - 2; i >= 0; --i)
	{
		result[i] = d[i] - c[i] * result[i + 1];
	}
}

void P3D::P3D_thomas::test_thomas_3()
{
	std::vector<double> top(3);  // top row, starts from 0 to Nx-1
	std::vector<double> middle(3);  // middle row, 0 to Nx
	std::vector<double> bottom(3);  // bottom row, 1 to Nx
	std::vector<double> f(3);  // equation`s right part
	std::vector<double> result(3);
	
	top[0]    = 1.0; top[1]    = 2.0; top[2]    = 3.0;
	bottom[0] = 3.0; bottom[1] = 2.0; bottom[2] = 0.0;
	middle[0] = 4.0; middle[1] = 5.0; middle[2] = 6.0;
	f[0]      = 1.0; f[1]      = 5.0; f[2]      = 4.0;
	
	thomas(top, middle, bottom, f, result);
	
	double eps0 = middle[0] * result[0] + top[0] * result[1] - f[0];  // i = 0
	double eps1 = bottom[1] * result[0] + middle[1] * result[1] + top[1] * result[2] - f[1];  // i = 1
	double eps2 = bottom[2] * result[1] + middle[2] * result[2] - f[2];  // i = 2
	printf("Eps: %lf\n", eps0);
	printf("Eps: %lf\n", eps1);
	printf("Eps: %lf\n", eps2);
}

void P3D::P3D_thomas::thomas
	(const std::vector<double> &A,
	 const std::vector<double> &B,
	 const std::vector<double> &C,
	 const std::vector<double> &D,
	 const std::vector<double> &E,
	 const std::vector<double> &F,
	       std::vector<double> &result)
{
	const size_t size = result.size();
	
#ifndef NDEBUG
	for(size_t i = 0; i < size; ++i)
	{
		if((fabs(A[i]) + fabs(B[i]) + fabs(D[i]) + fabs(E[i])) > fabs(C[i]))
		{
			printf("System matrix is not diagonally dominant!\n");
			exit(1);
		}
	}
#endif
	
	std::vector<double> alpha(size);
	std::vector<double> beta(size);
	std::vector<double> gamma(size);
	
	alpha[0] = - D[0] / C[0];
	beta [0] = - E[0] / C[0];
	gamma[0] =   F[0] / C[0];
	
	alpha[1] = - (B[1] * beta[0] + D[1])  / (B[1] * alpha[0] + C[1]);
	beta [1] = -  E[1]                    / (B[1] * alpha[0] + C[1]);
	gamma[1] = - (B[1] * gamma[0] - F[1]) / (B[1] * alpha[0] + C[1]);
	
	for(size_t i = 2; i < size; ++i)
	{
		double denominator = A[i] * (alpha[i-2] * alpha[i-1] + beta[i-2]) + B[i] * alpha[i-1] + C[i];
		
		alpha[i] = - (A[i] * alpha[i-2] * beta[i-1] + B[i] * beta[i-1] + D[i])                       / denominator;
		beta[i]  = -  E[i]                                                                           / denominator;
		gamma[i] = - (A[i] * alpha[i-2] * gamma[i-1] + A[i] * gamma[i-2] + B[i] * gamma[i-1] - F[i]) / denominator;
	}
	
	result.back() = gamma.back();
	result[size - 2] = alpha[size - 2] * result[size - 1] + gamma[size - 2];
	for (int i = (int)size - 3; i >= 0; --i)
	{
		result[i] = alpha[i] * result[i+1] + beta[i] * result[i+2] + gamma[i];
	}
}

void P3D::P3D_thomas::test_thomas_5()
{
	const size_t size = 5;
	
	std::vector<double> top_top(size);
	std::vector<double> top(size);
	std::vector<double> middle(size);
	std::vector<double> bot(size);
	std::vector<double> bot_bot(size);
	std::vector<double> f(size);
	std::vector<double> result(size);
	
	for(size_t i = 0; i < size; ++i)
	{
		top_top[i] = 1.0;
		top[i] = 2.0;
		middle[i] = 6.0;
		bot[i] = 1.0;
		bot_bot[i] = 1.0;
		f[i] = 5.0;
	}
	
	thomas(bot_bot, bot, middle, top, top_top, f, result);
	
	double eps0 = middle[0] * result[0] + top[0] * result[1] + top_top[0] * result[2] - f[0];  // i = 0
	double eps1 = bot[1] * result[0] + middle[1] * result[1] + top[1] * result[2] + top_top[1] * result[3] - f[1];  // i = 1
	double eps2 = bot_bot[2] * result[0] + bot[2] * result[1] + middle[2] * result[2] + top[2] * result[3] + top_top[2] * result[4] - f[2];  // i = 2
	double eps3 = bot_bot[3] * result[1] + bot[3] * result[2] + middle[3] * result[3] + top[3] * result[4] - f[3];  // i = 3
	double eps4 = bot_bot[4] * result[2] + bot[4] * result[3] + middle[4] * result[4] - f[4];  // i = 4
	printf("Eps: %lf\n", eps0);
	printf("Eps: %lf\n", eps1);
	printf("Eps: %lf\n", eps2);
	printf("Eps: %lf\n", eps3);
	printf("Eps: %lf\n", eps4);
	
	for(int j = 0; j < size; ++j)
	{
		printf("Res: %lf\n", result[j]);
	}
}

P3D::P3D_thomas::P3D_thomas(){};
P3D::P3D_thomas::~P3D_thomas(){};
