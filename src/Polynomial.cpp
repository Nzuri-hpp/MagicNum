#include "stdafx.h"
#include "Polynomial.h"


Polynomial::Polynomial() : n(1), coeff(1,0)
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
	return coeff[i];
}

double Polynomial::safe_at(int i)
{
	if (i >= n || i < 0) return 0.;

	return coeff.at(i);
}

double Polynomial::at_end()
{
	return coeff.back();
}

double Polynomial::at_begin()
{
	return coeff.front();
}

Polynomial Polynomial::der()
{
	if (n == 0) return Polynomial(0);
	if (n == 1) return Polynomial::zero();

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
	if (n == 0) return std::numeric_limits<double>::quiet_NaN();

	double acc = 0;
	const auto c_cf = coeff.data();
	for_dim(i)
	{
		acc += c_cf[i]*pow(x, i);
	}

	return acc;
}

double Polynomial::horners_method(double x)
{
	if (coeff.size() == 0)
	{
		throw error(errcode::EMPTY_SIZE);
	}

	if (x == std::numeric_limits<double>::infinity() || x == (-1)*std::numeric_limits<double>::infinity())
	{
		double back_coeff = coeff.back();
		size_t k = coeff.size() - 1;
		if (back_coeff == 0) return 0.;
		if (k == 0)          return back_coeff;

		if (x == std::numeric_limits<double>::infinity())
		{
			if (back_coeff > 0.)
			{
				return std::numeric_limits<double>::infinity();
			}
			if (back_coeff < 0.)
			{
				return (-1.)*std::numeric_limits<double>::infinity();
			}
		}
		else
		{
			if (back_coeff > 0. && k % 2 == 0)
			{
				return std::numeric_limits<double>::infinity();
			}
			if (back_coeff < 0. && k % 2 == 0)
			{
				return (-1)*std::numeric_limits<double>::infinity();
			}
			if (back_coeff > 0. && k % 2 == 1)
			{
				return (-1)*std::numeric_limits<double>::infinity();
			}
			if (back_coeff < 0. && k % 2 == 1)
			{
				return std::numeric_limits<double>::infinity();
			}
		}
	}

	double b = coeff.back();
	const auto c_cf = coeff.data();
	for (int i = coeff.size() - 2; i >= 0; --i)
	{
		b = c_cf[i] + b * x;
	}

	return b;
}

Polynomial Polynomial::operator/(Polynomial & right)
{
	return sub(right).first;
}

Polynomial Polynomial::mod(Polynomial & right)
{
	return sub(right).second;
}

pair<Polynomial, Polynomial> Polynomial::sub(Polynomial & right)
{
	if (degree() < right.degree()) return std::make_pair(Polynomial::zero(), (*this));
	if (right.degree() == 0) return std::make_pair((1/right.coeff.front()*(*this)),Polynomial::zero());
	
	Polynomial remainder = *this;
	Polynomial result = Polynomial::zero();

	while (remainder.degree() >= right.degree())
	{
		double C = remainder.at_end()/right.at_end();
		Polynomial xp = x_pow(remainder.degree() - right.degree());;
		result = result + C * xp;

		remainder = remainder - C * xp * right;
		remainder.restrict(DBL_MIN);
	}

	return make_pair(result, remainder);
}

size_t Polynomial::degree()
{
	if (n == 0) throw error(errcode::EMPTY_SIZE);
	return n-1;
}

double Polynomial::restrict(double epsilon)
{
	if (n == 0) return 0;
	if (n == 1) return 0;

	double eps_r = 0;
	size_t k = n - 1;

	while (abs(this->at(k)) <= epsilon)
	{
		eps_r = abs(this->at(k));
		--k;
	}
	n = k + 1;

	if (n == 0)
	{
		n = 1;
		coeff.at(0) = 0.;
	}

	coeff.resize(n);

	return eps_r;
}

Polynomial Polynomial::x_pow(size_t pow)
{
	vector<double> cf(pow + 1, 0);
	cf.at(pow) = 1.;
	return Polynomial(move(cf));
}

size_t Polynomial::sturms_method(double a, double b)
{
	if (degree() < 1) return 0;

	vector<double> sturms_variate_val_a(degree() + 1);
	vector<double> sturms_variate_val_b(degree() + 1);

	Polynomial p_der = der();

	sturms_variate_val_a[0] = this->horners_method(a);
	sturms_variate_val_a[1] = p_der.horners_method(a);

	sturms_variate_val_b[0] = this->horners_method(b);
	sturms_variate_val_b[1] = p_der.horners_method(b);

	Polynomial p_o_old = *this;
	Polynomial p_old(move(p_der));

	for (size_t i = 2; i <= degree(); ++i)
	{
		Polynomial sturms_variate = (-1)*p_o_old.mod(p_old);
		
		sturms_variate_val_a[i] = sturms_variate.horners_method(a);
		sturms_variate_val_b[i] = sturms_variate.horners_method(b);

		p_o_old = move(p_old);
		p_old   = move(sturms_variate);
	}

	size_t sturms_val_a = 0;
	size_t sturms_val_b = 0;

	auto sign = [](double x) {if (x == 0) return 0; return x > 0 ? +1 : -1; };

	int sign_a = sign(sturms_variate_val_a[0]);
	int sign_b = sign(sturms_variate_val_b[0]);
	for (size_t i = 1; i <= degree(); ++i)
	{
		int _sai = sign(sturms_variate_val_a[i]);
		int _sbi = sign(sturms_variate_val_b[i]);

		if (sign_a == 0) sign_a = _sai;
		if (sign_b == 0) sign_b = _sbi;

		if (sign_a != _sai && _sai != 0)
		{
			++sturms_val_a;
			sign_a = _sai;
		}
		if (sign_b != _sbi && _sbi != 0)
		{
			++sturms_val_b;
			sign_b = _sbi;
		}
	}

	return sturms_val_a - sturms_val_b;
}


double Polynomial::pricisie_root(volatile double a, volatile double b)
{
	auto sign = [](double x) {if (x == 0) return 0; return x > 0 ? +1 : -1; };

	double pa = horners_method(a);
	double pb = horners_method(b);

	//if (sign(pa) == sign(pb)) throw error(errcode::EQUALS_SIGNS);

	if (pa == 0.) return a;
	if (pb == 0.) return b;

	size_t i = 0;
	double ab_old = DBL_MAX;
	while ((a + b) / 2 != ab_old)
	{
		++i;
		double c  = (a + b) / 2;
		double pc = horners_method(c);

		if (pc == 0.) return c;

		if (sign(pc) != sign(pa))
		{
			b  = c;
			pb = pc;
		}
		else if (sign(pc) != sign(pb))
		{
			a  = c;
			pa = pc;
		}

		ab_old = c;
	}

	return a;
}

vector<double> Polynomial::pricisie_all_roots(double a, double b, double base_step)
{
	if (a > b) swap(a, b);
	if (a == b)
	{
		if (horners_method(a) == 0.)
			return vector<double>(1, a);
		else
			return vector<double>(0);
	}

	size_t all_r = sturms_method(a - DBL_EPSILON, b + DBL_EPSILON);
	if (all_r == 0) return vector<double>(0);
	if (all_r == 1) return vector<double>(1, pricisie_root(a, b));

	vector<double> Is(all_r + 1);
	Is.front() = a;
	Is.back()  = b;

	double L = (b - a) / static_cast<double>(all_r);
	for (size_t i = 1; i < Is.size() - 1; ++i)
	{
		Is[i] = a + i*L;
	}

	vector<size_t> r_amount(all_r);

	for (size_t i = 0; i < all_r; ++i)
	{
		r_amount[i] = sturms_method(Is[i], Is[i + 1]);
	}


	for (size_t i = 1; i < Is.size() - 1; ++i)
	{
		double step = base_step;
		size_t amount_old;
		if (Is[i - 1] <= Is[i])
		{
			Is[i] = Is[i - 1] + DBL_EPSILON;
			r_amount[i - 1] = sturms_method(Is[i - 1], Is[i]);
		}

		while (r_amount[i - 1] != 1)
		{
			double direct = +1.;
			if (r_amount[i - 1] > 0) direct = -1.;

			double _Is_new = Is[i] + direct*step;

			if (_Is_new <= Is[i - 1])
			{
				step /= 2.;
				continue;
			}
			Is[i] = _Is_new;

			amount_old = r_amount[i - 1];
			r_amount[i - 1] = sturms_method(Is[i - 1], Is[i]);

			if (amount_old < r_amount[i - 1]) step /= 2.;

			if (step <= DBL_MIN) step = base_step;
		}
	}

	vector<double> roots(all_r);

	for (size_t i = 0; i < all_r; ++i)
	{
		roots[i] = pricisie_root(Is[i], Is[i+1]);
	}

	return move(roots);
}

Polynomial Polynomial::zero()
{
	auto p = Polynomial(1);
	p.coeff[0] = 0;
	return move(p);
}

Polynomial Polynomial::constant(double c)
{
	auto p = Polynomial(1);
	p.coeff[0] = c;
	return move(p);
}

void Polynomial::take_normal()
{
	double max = *std::max_element(coeff.begin(), coeff.end());

	for (auto& c : coeff)
	{
		c /= max;
	}
}

double Polynomial::drop_back()
{
	if (coeff.size() == 0) throw error(errcode::EMPTY_SIZE);
	if (coeff.size() == 1) return coeff.front();

	double b = coeff.back();
	coeff.resize(coeff.size() - 1);
	n = coeff.size();

	return b;
}
