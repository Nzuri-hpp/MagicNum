#pragma once
#include "QVectorT.h"
#include "QMatrix.h"
#include "Polynomial.h"

struct point
{
	double x;
	double y;
	point() : x(0), y(0) {}
	point(double _x, double _y) : x(_x), y(_y) {}
};

typedef function<double(double)> func_Rt;
typedef function<QVector(QVector)> func_RNt;
typedef function<QMatrix(QVector)> Jfunc_RNt;
typedef function<double(QVector)>  func_RNt_double;

QVector JacobiIterations(QMatrix& A, QVector& b, QVector x0, bool std_init = false, double epsilon = 1e-7, size_t MAX_STEPS = 1000000);
QVector ZeidelIterations(QMatrix& A, QVector& b, QVector x0, bool std_init = false, double epsilon = 1e-7, size_t MAX_STEPS = 1000000);
QVector SimpleIterations(QMatrix& B, QVector& c, QVector x0, bool std_init = false, double epsilon = 1e-7, size_t MAX_STEPS = 1000000);

double eval(string expression_string, double x);
exprtk::expression<double> fast_eval(string expression_string, double* x);

double bisection(func_Rt f, double a, double b, double epsilon = 1e-5, size_t MAX_STEPS = 1000000);
double regula_falsi(func_Rt f, double a, double b, double epsilon = 1e-5, size_t MAX_STEPS = 1000000);
double newtone(func_Rt f, func_Rt der_f, double a, double b, double multiplicity = 1, double epsilon = 1e-5, size_t MAX_STEPS = 1000000);
double combnewtone(func_Rt f, func_Rt der_f, double a, double b, double multiplicity = 1, double epsilon = 1e-5, size_t MAX_STEPS = 1000000);
double simpleit(func_Rt f, double a, double b, double epsilon = 1e-5, size_t MAX_STEPS = 1000000);

QVector PLU_slove(QMatrix& L, QMatrix& U, QVector& b);

//Переделать для повышения производительности все eval функции
double eval_RN(string expression_string, QVector x);

exprtk::expression<double> fast_eval_RN(string expression_string, QVector *x);

QVector newtone_slove_sys(func_RNt F, Jfunc_RNt J, QVector x0, double epsilon = 1e-5, size_t MAX_STEPS = 1000000);

double  minimize_golden_section(func_Rt F, double a, double b, double epsilon = 1e-5, size_t max_steps = 1000000);

QVector minimize_fast_grad(func_RNt_double F, func_RNt Grad, QVector x0, double epsilon = 1e-5, size_t max_steps = 1000000, double max_lambda = 50.);

Polynomial interpolate(vector<point>& points);

Polynomial polynomial_approximation(vector<point>& points, Polynomial& p0, double epsilon = 1e-5, size_t max_steps = 1000000, double max_lambda = 50.);