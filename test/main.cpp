#include <iostream>
#include <math.h>
#include <array>
#include <vector>

int size = 100000;
double alpha = 0.5;
double dx = 1.0 / size;
std::vector<double> grunvald_coeffs(size);  // массив коэффициентов Грюнвальда-Летникова
std::vector<double> V(size);

double fracD()
{
	double result = 0.0;
	
	for(int j = 0; j < size; ++j)
	{
		result += grunvald_coeffs[j] * 1.0;
	}
	
	return pow(dx, - alpha) * result;
}

int main(){
	std::cout << "Hello, World!" << std::endl;
	
	grunvald_coeffs[0] = 1.0;
	for (size_t i = 1; i < grunvald_coeffs.size(); ++i)
	{
		grunvald_coeffs[i] = grunvald_coeffs[i - 1] * (alpha - i + 1.0) / i;
	}
	
	for(size_t i = 0; i < grunvald_coeffs.size(); ++i)
	{
		grunvald_coeffs[i] *= pow(- 1.0, i);
	}
	
	double result = fracD();
	
	printf("%lf\n", result);
	
	return 0;
}
