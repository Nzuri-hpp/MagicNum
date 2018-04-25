#include "stdafx.h"
#include "Polynomial.h"


Polynomial::Polynomial() : n(1), coeff({ 0 })
{
}

Polynomial::Polynomial(size_t _n) : n(_n)
{
	coeff.resize(n);
}

Polynomial::Polynomial(const vector<double> & cf) : n(cf.size()), coeff(cf)
{
}

Polynomial::Polynomial(QVector & cf) : n(cf.size())
{
	coeff.resize(cf.size());

	for (size_t i = 0; i < cf.size(); ++i)
	{
		coeff.at(i) = cf.at(i);
	}

}

Polynomial::Polynomial(const Polynomial & right) : n(right.n), coeff(right.coeff)
{
}

Polynomial::Polynomial(const Polynomial && right) : n(right.n), coeff(right.coeff)
{
}


Polynomial::~Polynomial()
{
}

const vector<double>& Polynomial::getCoeff()
{
	return coeff;
}

void Polynomial::setCoeff(vector<double>& newCoeff)
{
	if (newCoeff.size() != n) throw error(errcode::DIFFERENCE_SIZE);

	coeff = newCoeff;
}

Polynomial Polynomial::operator*(Polynomial & right)
{
	size_t dim = max(n, right.n);
	Polynomial res(n+right.n-1);

	for(int p = res.n-1; p >=0; --p)
	{
		for (int i = 0; i < dim; ++i)
		{
			res.coeff.at(p) += safe_at(i)*right.safe_at(p - i);
		}
	}

	return res;
}

Polynomial& Polynomial::operator*=(Polynomial & right)
{
	return (*this) = operator*(right);
}

Polynomial& Polynomial::operator=(Polynomial & right)
{
	n = right.n;
	coeff = right.coeff;

	return *this;
}

Polynomial Polynomial::operator*(double right)
{
	Polynomial res(n);

	res.setCoeff(coeff);

	for (auto& cf : res.coeff)
	{
		cf *= right;
	}

	return res;
}


Polynomial& Polynomial::operator+=(Polynomial & right)
{
	return (*this) = operator+(right);;
}

Polynomial Polynomial::operator+(double right)
{
	Polynomial res(n);

	res.setCoeff(coeff);

	res.coeff.at(0) += right;

	return res;
}

Polynomial Polynomial::operator-(double right)
{
	Polynomial res(n);

	res.setCoeff(coeff);

	res.coeff.at(0) -= right;

	return res;
}

double& Polynomial::at(size_t i)
{
	return coeff.at(i);
}

double Polynomial::safe_at(int i)
{
	if (i >= n || i < 0) return 0.;

	return coeff.at(i);
}

Polynomial Polynomial::der()
{
	if (n == 0) return Polynomial(0);
	if (n == 1) return Polynomial({ 0 });

	Polynomial derivate(n-1);

	for (size_t i = 0; i < n - 1; ++i)
	{
		derivate.coeff.at(i) = (i + 1)*coeff.at(i + 1);
	}

	return derivate;
}

Polynomial Polynomial::der_by_coeff(size_t i)
{
	Polynomial derc(i + 1);

	derc.coeff.at(i) = 1.;

	return derc;
}

Polynomial Polynomial::operator-()
{
	Polynomial res(n);

	res.setCoeff(coeff);

	for_dim(i)
	{
		res.coeff.at(i) *= -1.;
	}

	return res;
}

Polynomial Polynomial::operator-(Polynomial & right)
{
	size_t dim = max(right.n, n);
	Polynomial res(dim);

	if (right.n <= n)
	{
		res.setCoeff(coeff);

		for(size_t i = 0; i < right.n; ++i)
		{
			res.coeff.at(i) -= right.coeff.at(i);
		}
	}
	else
	{
		res.setCoeff((-right).coeff);

		for_dim(i)
		{
			res.coeff.at(i) += coeff.at(i);
		}
	}
	return res;
}

Polynomial Polynomial::operator+(Polynomial & right)
{
	size_t dim = max(right.n, n);
	Polynomial res(dim);

	if (right.n <= n)
	{
		res.setCoeff(coeff);

		for (size_t i = 0; i < right.n; ++i)
		{
			res.coeff.at(i) += right.coeff.at(i);
		}
	}
	else
	{
		res.setCoeff(right.coeff);

		for_dim(i)
		{
			res.coeff.at(i) += coeff.at(i);
		}
	}
	return res;
}

string Polynomial::serialize()
{
	stringstream ss;
	string serial;

	ss.precision(3);
	ss.setf(ios::showpos);
	for_dim(i)
	{
		if(n-i-1 != 0) ss << coeff.at(n-i-1) << "x^" << n-i-1;
		else           ss << coeff.at(n-i-1);
	}
	ss >> serial;

	size_t ind;
	if (serial.at(0) == '+') ind = 1; else ind = 0;

	return serial.substr(ind);
}

double Polynomial::operator()(double x)
{
	double acc = 0;

	for_dim(i)
	{
		acc += coeff.at(i)*pow(x, i);
	}

	return acc;
}
