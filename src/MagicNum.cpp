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
/*	integration::Constructor<integration::CONSTRUCTOR::HERMITE> master;
	auto qh = master.construct(10000);
	double res = qh([](double x)->double {return std::sqrt(1-x*x); });
	std::cout << std::setprecision(50) << res << std::endl << std::endl;
	
	
	Polynomial pl({5, 4, 3, 2, 1, 1});
	std::cout << pl.serialize() << endl;
	double a;
	auto begin = std::chrono::high_resolution_clock::now();
	for (double i = 0; i < 1e+12; i+=1.)
	{
		a = pl(10*i);
	}
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << "ns" << std::endl;

	begin = std::chrono::high_resolution_clock::now();
	for (double i = 0; i < 1e+12; i += 1.)
	{
		a = pl.horners_method(10*i);
	}
	end = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << "ns" << std::endl;
	
	Polynomial pp({-42, 0, -12, 1});
	Polynomial pr({ -3, 1 });

	Polynomial rr = pp / pr;

	cout << rr.serialize() << endl << endl;

	cout << Polynomial({-1, 0 ,1}).sturms_method(-1, 1) << endl << endl;
	cout << Polynomial({ -1, 0 ,1 }).pricisie_root(-1.4, -0.5) << endl << endl;

	auto rv = Polynomial({ 24, -50 , 35, -10, 1 }).pricisie_all_roots(-300000,300000);

	for (auto rr : rv) cout << rr << ' ';
	cout << endl << endl;

	integration::Constructor<integration::CONSTRUCTOR::__ALG> alg;

	for(double i = 1.; i <= 20.; i+=1.)
		cout << alg.factorial(i) << " ";
	cout << endl << endl;

	cout << alg.sum_log(10) << " ";
	cout << endl << endl;

	alg.clear_buffers();
	cout << endl << endl;
*/
	integration::Constructor<integration::CONSTRUCTOR::LAGGER>  master_L;
	integration::Constructor<integration::CONSTRUCTOR::HERMITE> master_H;

	auto qr = master_H.construct(1000);

	auto f = [](double x)->double {return sqrt(1 - x * x)*exp(-x * x); };

	cout << qr(f) << endl;

	cout << integration::Constructor<integration::CONSTRUCTOR::LAGGER>::precision_test(1000) << endl;

	system("pause");
    return 0;
}

