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
		enum codes { DIFFERENCE_SIZE, BAD_INDEX, EQUALS_SIGNS, EMPTY_SIZE };
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
	double at_end();
	double at_begin();
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

	Polynomial operator/(Polynomial & right);
	Polynomial mod(Polynomial & right);
	pair<Polynomial, Polynomial> sub(Polynomial & right);

	size_t degree();
	double restrict(double epsilon = DBL_EPSILON);
	static Polynomial x_pow(size_t pow);

	size_t sturms_method(double a, double b);

	double pricisie_root(double a, double b);
	vector<double> pricisie_all_roots(double a, double b, double base_step = 10.);

	static Polynomial zero();
	static Polynomial constant(double c);

	void take_normal();
	double drop_back();

};

