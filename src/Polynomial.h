#pragma once
#include "QVectorT.h"

using namespace std;

class Polynomial
{
#define for_dim(a) for(size_t (a)=0;a<n;++a)

	vector<double> coeff;
	size_t n;
public:
	struct error
	{
		const size_t code;
		error(size_t code_) : code(code_) {}
	};
	struct errcode
	{
		enum codes { DIFFERENCE_SIZE, BAD_INDEX };
	};

	Polynomial();
	Polynomial(size_t _n);
	Polynomial(const vector<double>& cf);
	Polynomial(QVector& cf);
	Polynomial(const Polynomial& right);
	Polynomial(const Polynomial&& right);
	~Polynomial();

	const vector<double>& getCoeff();
	void  setCoeff(vector<double>& newCoeff);

	double operator()(double x);
	string serialize();

	Polynomial operator*(Polynomial& right);
	Polynomial& operator*=(Polynomial& right);
	Polynomial& operator=(Polynomial& right);
	//Polynomial operator=(Polynomial& right);
	Polynomial operator*(double right);
	Polynomial operator+(Polynomial& right);
	Polynomial& operator+=(Polynomial& right);
	Polynomial operator+(double right);
	Polynomial operator-(Polynomial& right);
	Polynomial operator-(double right);
	Polynomial operator-();
	double& at(size_t i);
	double  safe_at(int i);
	Polynomial der();
	Polynomial der_by_coeff(size_t i);
	friend Polynomial operator*(double left, Polynomial& right)
	{
		return right.operator*(left);
	}
	friend Polynomial operator+(double left, Polynomial& right)
	{
		return right.operator+(left);
	}
	friend Polynomial operator-(double left, Polynomial& right)
	{
		return (-right).operator+(left);
	}

	double horners_method(double x);
};

