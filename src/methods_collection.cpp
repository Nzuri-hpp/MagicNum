#include "stdafx.h"
#include "methods_collection.h"

QVector JacobiIterations(QMatrix& A, QVector& b, QVector x0, bool std_init, double epsilon, size_t MAX_STEPS)
{
	size_t n = A.size();

	if (A.size() != b.size() || b.size() != x0.size())
	{
		throw string("Error! Difference size of matrix and vector!");
	}

	QMatrix L(n);
	QMatrix invD(n);
	QMatrix R(n);

	L.zero();
	R.zero();
	invD.zero();

	for (size_t i = 0; i < n; i++)
	{
		if (abs(A.at(i, i)) < epsilon) throw string("Error! Zeros on diag!");

		invD.at(i, i) = 1 / A.at(i, i);
		for (size_t j = i + 1; j < n; j++)
		{
			L.at(i, j) = A.at(i, j);
			R.at(j, i) = A.at(j, i);
		}
	}

	QMatrix B(n);
	QVector c(n);

	B = (-1)*invD*(L + R);
	c = invD * b;

	if (std_init)
	{
		x0 = c;
	}

	double norm_B = B.cubic_norm();
	double norm_c = c.cubic_norm();
	double norm_x0 = x0.cubic_norm();

	cout << "B: " << B << endl;
	if (A.diagonal_dominance())
	{
		size_t k = 1;
		while (pow(norm_B, k)*(norm_x0 + norm_c / (1 - norm_B)) > epsilon)
		{
			++k;
		}
		cout << "Решение сходится для любого x0, априорное количество шагов: " << k << "  для epsilon = " << epsilon << endl;
	}
	else
	{
		cout << "Сходимость зависит от начальных данных!" << endl;
		throw string("Error! Convergence depends on the initial data!");
	}

	QVector x(n);
	QVector x_next(n);
	x = x0;

	QVector delta(n);

	size_t k = 0;

	do {
		x_next = B * x + c;
		delta = x_next - x;

		swap(x, x_next);
		k++;
	} while (delta.cubic_norm() >= (1 - norm_B) / norm_B * epsilon && k <= MAX_STEPS);

	cout << "Понадобилось шагов: " << k << endl;

	return x;
}

QVector ZeidelIterations(QMatrix& A, QVector& b, QVector x0, bool std_init, double epsilon, size_t MAX_STEPS)
{
	size_t n = A.size();
	if (A.size() != b.size() || b.size() != x0.size())
	{
		throw string("Error! Difference size of matrix and vector!");
	}
	if (std_init)
	{
		QMatrix invD(n);
		invD.zero();

		for (size_t i = 0; i < n; i++)
		{
			if (abs(A.at(i, i)) < epsilon) throw string("Error! Zeros on diag!");
			invD.at(i, i) = 1 / A.at(i, i);
		}

		x0 = invD * b;
	}
	if (A.diagonal_dominance())
	{
	}
	else
	{
		cout << "Сходимость зависит от начальных данных!" << endl;
		throw string("Error! Convergence depends on the initial data!");
	}

	QVector x = x0;
	QVector p(n);
	size_t k = 0;

	auto converge = [epsilon, n](QVector& xk, QVector& xkp)
	{
		double norm = 0;
		for (int i = 0; i < n; i++)
		{
			norm += (xk.at(i) - xkp.at(i))*(xk.at(i) - xkp.at(i));
		}
		if (sqrt(norm) >= epsilon)
			return false;
		return true;
	};

	do
	{
		p = x;

		for (size_t i = 0; i < n; i++)
		{
			double var = 0;
			for (size_t j = 0; j < i; j++)
				var += (A.at(i, j) * x.at(j));
			for (size_t j = i + 1; j < n; j++)
				var += (A.at(i, j) * p.at(j));
			x.at(i) = (b.at(i) - var) / A.at(i, i);
		}
		k++;
	} while (!converge(x, p) && k <= MAX_STEPS);

	cout << "Понадобилось шагов: " << k << "  для epsilon = " << epsilon << endl;

	return x;
}

QVector SimpleIterations(QMatrix& B, QVector& c, QVector x0, bool std_init, double epsilon, size_t MAX_STEPS)
{
	size_t n = B.size();
	if (B.size() != c.size() || c.size() != x0.size())
	{
		throw string("Error! Difference size of matrix and vector!");
	}
	if (std_init)
	{
		x0 = c;
	}
	double norm_B = B.cubic_norm();
	double norm_x0 = x0.cubic_norm();
	double norm_c = c.cubic_norm();
	if (norm_B < 1)
	{
		size_t k = 1;
		while (pow(norm_B, k)*(norm_x0 + norm_c / (1 - norm_B)) > epsilon)
		{
			++k;
		}
		cout << "Решение сходится для любого x0, априорное количество шагов: " << k << "  для epsilon = " << epsilon << endl;
	}
	else
	{
		cout << "Сходимость, возможно, зависит от начальных данных!" << endl;
	}

	QVector x(n);
	QVector x_next(n);
	x = x0;

	QVector delta(n);

	size_t k = 0;

	do {
		x_next = B * x + c;
		delta = x_next - x;

		swap(x, x_next);
		k++;
	} while (delta.cubic_norm() >= (1 - norm_B) / norm_B * epsilon && k <= MAX_STEPS);

	cout << "Понадобилось шагов: " << k << endl;

	return x;
}

double eval(string expression_string, double x)
{
	typedef double T;
	typedef exprtk::symbol_table<T> symbol_table_t;
	typedef exprtk::expression<T>     expression_t;
	typedef exprtk::parser<T>             parser_t;

	symbol_table_t symbol_table;
	symbol_table.add_variable("x", x);
	symbol_table.add_constants();

	expression_t expression;
	expression.register_symbol_table(symbol_table);

	parser_t parser;

	if (!parser.compile(expression_string, expression))
	{
		throw string("Error in expression!");
	}

	return expression.value();
}

exprtk::expression<double> fast_eval(string expression_string, double* x)
{
	typedef double T;
	typedef exprtk::symbol_table<T> symbol_table_t;
	typedef exprtk::expression<T>     expression_t;
	typedef exprtk::parser<T>             parser_t;

	symbol_table_t symbol_table;
	symbol_table.add_variable("x", *x);
	symbol_table.add_constants();

	expression_t expression;
	expression.register_symbol_table(symbol_table);

	parser_t parser;

	if (!parser.compile(expression_string, expression))
	{
		throw string("Error in expression!");
	}

	return expression;
}

double bisection(func_Rt f, double a, double b, double epsilon, size_t MAX_STEPS)
{
	if (a == b)
	{
		throw string("Error! Degenerating interval (a,a)!");
	}
	if (a > b) swap(a, b);
	if (f(a)*f(b) >= 0)
	{
		throw string("Error! For x in (a, b): f(a)f(b) >= 0");
	}

	size_t k = 0;
	do
	{
		k++;
		double c = (a + b) / 2.;

		if (f(c) == 0.) return c;

		if (f(a)*f(c) < 0) b = c;
		if (f(b)*f(c) < 0) a = c;
	} while (abs(a - b) >= epsilon && k <= MAX_STEPS);

	cout << "Понадобилось итераций: " << k << " при epsilon = " << epsilon << endl;

	return (a + b) / 2.;
}

double regula_falsi(func_Rt f, double a, double b, double epsilon, size_t MAX_STEPS)
{
	if (abs(a - b) <= 1000 * epsilon)
	{
		throw string("Error! Too small interval for this epsilon!");
	}
	if (a > b) swap(a, b);
	if (f(a)*f(b) >= 0)
	{
		throw string("Error! For x in (a, b): f(a)f(b) >= 0");
	}

	double h = abs(b - a) / 100;
	if (h > 1000 * epsilon) h = 1000 * epsilon;
	auto dder_f = [&f, h](double x) {return (f(x + h) - 2 * f(x) + f(x - h)) / h; };

	double fix;
	double start;
	if (f(b) > 0 && dder_f(b) > 0 || f(b) < 0 && dder_f(b) < 0)
	{
		fix = b;
		start = a;

	}
	else
	{
		fix = a;
		start = b;
	}

	double f_fix = f(fix);
	double x = start;
	double x_old;

	size_t k = 0;
	do
	{
		swap(x, x_old);
		x = x_old - (f(x_old)*(x_old - fix)) / (f(x_old) - f_fix);
		k++;
	} while (abs(x - x_old) >= epsilon && k <= MAX_STEPS);

	cout << "Понадобилось итераций: " << k << " при epsilon = " << epsilon << endl;

	return x;
}

double newtone(func_Rt f, func_Rt der_f, double a, double b, double multiplicity, double epsilon, size_t MAX_STEPS)
{
	if (abs(a - b) <= 1000 * epsilon)
	{
		throw string("Error! Too small interval for this epsilon!");
	}
	if (a > b) swap(a, b);
	if (f(a)*f(b) >= 0)
	{
		throw string("Error! For x in (a, b): f(a)f(b) >= 0");
	}

	double h = abs(b - a) / 100;
	if (h > 1000 * epsilon) h = 1000 * epsilon;
	auto dder_f = [&der_f, h](double x) {return (der_f(x + h) - der_f(x)) / h; };

	double start;
	if (f(b) > 0 && dder_f(b) > 0 || f(b) < 0 && dder_f(b) < 0)
	{
		start = b;
	}
	else
	{
		start = a;
	}

	double x = start;
	double x_old;

	size_t k = 0;
	do
	{
		swap(x, x_old);

		x = x_old - multiplicity * f(x_old) / der_f(x_old);
		k++;

	} while (abs(x - x_old) >= epsilon && k <= MAX_STEPS);

	cout << "Понадобилось итераций: " << k << " при epsilon = " << epsilon << endl;

	return x;
}

double combnewtone(func_Rt f, func_Rt der_f, double a, double b, double multiplicity, double epsilon, size_t MAX_STEPS)
{
	if (abs(a - b) <= 100 * epsilon)
	{
		throw string("Error! Too small interval for this epsilon!");
	}
	if (a > b) swap(a, b);
	if (f(a)*f(b) >= 0)
	{
		throw string("Error! For x in (a, b): f(a)f(b) >= 0");
	}

	double h = abs(b - a) / 10;
	if (h > 100 * epsilon) h = 100 * epsilon;
	auto dder_f = [&f, h](double x) {return (f(x + h) - 2 * f(x) + f(x - h)) / h; };

	double start_newtone;
	double start_line;
	double fix;
	if (f(b) > 0 && dder_f(b) > 0 || f(b) < 0 && dder_f(b) < 0)
	{
		start_newtone = b;
		start_line = a;
		fix = b;
	}
	else
	{
		start_newtone = a;
		start_line = b;
		fix = a;
	}

	double x1 = start_newtone;
	double x1_old;

	double x2 = start_line;
	double x2_old;

	double f_fix = f(fix);

	size_t k = 0;
	do
	{
		k++;
		swap(x1, x1_old);
		swap(x2, x2_old);

		x2 = x2_old - (f(x2_old)*(x2_old - fix)) / (f(x2_old) - f_fix);
		x1 = x1_old - multiplicity * f(x1_old) / der_f(x1_old);

	} while (abs(x1 - x2) >= epsilon && k <= MAX_STEPS);

	cout << "Понадобилось итераций: " << k << " при epsilon = " << epsilon << endl;

	return (x1 + x2) / 2;
}

double simpleit(func_Rt f, double a, double b, double epsilon, size_t MAX_STEPS)
{
	double x = (a + b) / 2;
	double x_old;

	size_t k = 0;
	do
	{
		swap(x, x_old);
		x = f(x_old);
		k++;

	} while (abs(x - x_old) >= epsilon && k <= MAX_STEPS);

	cout << "Понадобилось итераций: " << k << " при epsilon = " << epsilon << endl;

	return x;
}

QVector PLU_slove(QMatrix& L, QMatrix& U, QVector& b)
{
	size_t n = L.size();
	if (L.size() != U.size() || U.size() != b.size())
	{
		throw string("Error! Difference size of matrix and vector!");
	}
	L.vector_renumerate(b);
	QVector y(b);

	for (size_t i = 0; i < n; i++)
	{
		double sum = 0;
		for (size_t k = 0; k < i; k++) sum += L.at(i, k)*y.at(k);
		y.at(i) -= sum;
	}

	QVector x(y);

	for (int i = n - 1; i >= 0; i--)
	{
		double sum = 0;
		for (size_t k = i + 1; k < n; k++) sum += U.at(i, k)*x.at(k);
		x.at(i) -= sum;
		x.at(i) /= U.at(i, i);
	}

	return x;
}

double eval_RN(string expression_string, QVector x)
{
	size_t n = x.size();

	typedef double T;
	typedef exprtk::symbol_table<T> symbol_table_t;
	typedef exprtk::expression<T>     expression_t;
	typedef exprtk::parser<T>             parser_t;

	symbol_table_t symbol_table;
	for (size_t i = 0; i < n; i++)
	{
		symbol_table.add_variable("x" + to_string(i + 1), x.at(i));
	}
	symbol_table.add_constants();

	expression_t expression;
	expression.register_symbol_table(symbol_table);

	parser_t parser;

	if (!parser.compile(expression_string, expression))
	{
		throw string("Error in expression!");
	}

	return expression.value();
}

exprtk::expression<double> fast_eval_RN(string expression_string, QVector *x)
{
	size_t n = x->size();

	typedef double T;
	typedef exprtk::symbol_table<T> symbol_table_t;
	typedef exprtk::expression<T>     expression_t;
	typedef exprtk::parser<T>             parser_t;

	symbol_table_t symbol_table;
	for (size_t i = 0; i < n; i++)
	{
		symbol_table.add_variable("x" + to_string(i + 1), x->at(i));
	}
	symbol_table.add_constants();

	expression_t expression;
	expression.register_symbol_table(symbol_table);

	parser_t parser;

	if (!parser.compile(expression_string, expression))
	{
		throw string("Error in expression!");
	}

	return expression;
}

QVector newtone_slove_sys(func_RNt F, Jfunc_RNt J, QVector x0, double epsilon, size_t MAX_STEPS)
{
	size_t n = F(x0).size();
	if (F(x0).size() != J(x0).size())
	{
		throw string("Error! Difference size of Jacobian and Function!");
	}
	QVector x_curr = x0;
	QVector x_old(n);
	QVector p_diff(n);
	size_t k = 0;
	do
	{
		k++;
		swap(x_curr, x_old);

		QMatrix Jfx = J(x_old);
		QVector rFx = (-1)*F(x_old);

		QMatrix L_Jfx(n);
		QMatrix U_Jfx(n);
		Jfx.PLU_decomposition(L_Jfx, U_Jfx, epsilon);
		p_diff = PLU_slove(L_Jfx, U_Jfx, rFx);
		x_curr = x_old + p_diff;

	} while (p_diff.euqlid_norm() >= epsilon && k <= MAX_STEPS);
	cout << "Понадобилось итераций: " << k << " при epsilon = " << epsilon << endl;

	return x_curr;

}

double minimize_golden_section(func_Rt F, double a, double b, double epsilon, size_t max_steps)
{
	if (a > b) swap(a, b);
	if (a == b) return a;

	double gold = (1 + sqrt(5)) / 2;

	size_t step = 0;
	while (abs(F(a) - F(b)) >= epsilon || abs(b - a) >= epsilon)
	{
		double x1 = b - (b - a) / gold;
		double x2 = a + (b - a) / gold;

		if (F(x1) > F(x2))
		{
			a = x1;
		}
		else
		{
			b = x2;
		}
		++step;
		if (step >= max_steps) break;
	}
	cout << "Понадобилось итераций: " << step << " при epsilon = " << epsilon << endl;

	return (a + b) / 2;
}

QVector minimize_fast_grad(func_RNt_double F, func_RNt Grad, QVector x, double epsilon, size_t max_steps, double max_lambda)
{
	size_t n = x.size();
	QVector x_new(n);
	QVector grad_x(n);

	size_t step = 0;

	func_Rt edge_F = [&x, &grad_x, &F](double lambda) -> double
	{
		return F(x - lambda * grad_x);
	};

	do
	{
		grad_x = Grad(x); // оптимизация вычислений
		double lambda = minimize_golden_section(edge_F, 0, max_lambda, epsilon, max_steps);

		x_new = x - lambda * grad_x;
		swap(x_new, x);

		++step;
		if (step >= max_steps) break;
	} while ((x_new - x).euqlid_norm() >= epsilon);

	cout << "\nСводка: понадобилось общих итераций: " << step << " при epsilon = " << epsilon << endl;
	return x;
}

Polynomial interpolate(vector<point>& points)
{
	vector<Polynomial> L;

	L.reserve(points.size());

	vector<double> ent = { 1 };
	for (size_t i = 0; i < points.size(); ++i)
	{
		L.emplace_back(ent);
	}

	for (size_t i = 0; i < points.size(); ++i)
	{
		for (size_t j = 0; j < points.size(); ++j)
		{
			if (i == j) continue;
			L.at(i) *= Polynomial({ -points.at(j).x, 1 }) * (1 / (points.at(i).x - points.at(j).x));
		}
	}

	Polynomial lagrang;
	for (size_t i = 0; i < points.size(); ++i)
	{
		lagrang += L.at(i)*points.at(i).y;
	}

	return lagrang;
}

Polynomial polynomial_approximation(vector<point>& points, Polynomial& p0, double epsilon, size_t max_steps, double max_lambda)
{
	QVector x0(p0.getCoeff());

	auto RSS = [&points](QVector& coeff) -> double
	{
		Polynomial poly(coeff);
		double err = 0;
		for (size_t i = 0; i < points.size(); ++i)
		{
			err += pow(points.at(i).y - poly(points.at(i).x), 2.);
		}

		return err;
	};

	auto RSS_der = [&points](QVector& coeff) -> QVector
	{
		Polynomial poly(coeff);

		QVector grad(coeff.size());
		grad.zero();

		for (size_t k = 0; k < grad.size(); ++k)
		{
			Polynomial poly_der_k = poly.der_by_coeff(k);
			for (size_t i = 0; i < points.size(); ++i)
			{
				grad.at(k) += -2 * (points.at(i).y - poly(points.at(i).x))*poly_der_k(points.at(i).x);
			}
		}

		return grad;
	};

	QVector optimal_coeff = minimize_fast_grad(RSS, RSS_der, x0, epsilon, max_steps, max_lambda);

	return Polynomial(optimal_coeff);
}
