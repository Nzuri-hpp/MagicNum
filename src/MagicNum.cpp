// MagicNum.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include "exprtk.hpp"
#include "methods_collection.h"
#include "Polynomial.h"
#include "integration.h"

#include <chrono>

int main()
{
	integration::Constructor<integration::CONSTRUCTOR::HERMITE> master;
	auto qh = master.construct(10000);
	double res = qh([](double x)->double {return std::sqrt(1-x*x); });
	std::cout << std::setprecision(50) << res << std::endl << std::endl;
	
	
	Polynomial pl({5, 4, 3, 2, 1, 1});
	std::cout << pl.serialize() << endl;
	double a;
	auto begin = std::chrono::high_resolution_clock::now();
	for (double i = 0; i < 1e+6; i+=1.)
	{
		a = pl(10*i);
	}
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << "ns" << std::endl;

	begin = std::chrono::high_resolution_clock::now();
	for (double i = 0; i < 1e+6; i += 1.)
	{
		a = pl.horners_method(10*i);
	}
	end = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << "ns" << std::endl;
	cout << a;
	system("pause");
    return 0;
}

