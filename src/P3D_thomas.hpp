//
// Created by mgreb on 11/3/16.
//

namespace P3D
{
	class P3D_thomas;
}

#pragma once

// Системные файлы
#include <vector>

class P3D::P3D_thomas
{
public:
	static void thomas
		(const std::vector<double> &top, 
		 const std::vector<double> &middle, 
		 const std::vector<double> &bottom,
		 const std::vector<double> &right,
		       std::vector<double> &result);
	
public:
	/**
	* 
	* @param A - самая нижняя диагональ матрицы 
	* @param B - ниже центральной
	* @param C - центральная диагональ
	* @param D - выше центральной
	* @param E - самая верхняя диагональ
	* @param F - правая часть
	* @param result 
	*/
	static void thomas
		(const std::vector<double> &A,
		 const std::vector<double> &B,
		 const std::vector<double> &C,
		 const std::vector<double> &D,
		 const std::vector<double> &E,
		 const std::vector<double> &F,
	         std::vector<double> &result);
	
	static void test_thomas_3();
	static void test_thomas_5();
private:
	P3D_thomas();
	~P3D_thomas();
};
