#pragma once
#include "QVectorT.h"
#include <algorithm>
#include <iostream>
#include <sstream>
using namespace std;
class QMatrix
{
#define for_dim(a) for(size_t (a) = 0; (a) < n; (a)++)
	double **src;
	const size_t n;

	QVectorT<size_t> transposition;
	size_t  transp_am = 0;
public:
	struct error
	{
		const size_t code;
		error(size_t code_) : code(code_)
		{
		}
	};
	struct errcode
	{
		enum codes { DIFFERENCE_SIZE, DEGENERATING_MATRIX, ZERO_ON_DIAG, BAD_INDEX };
	};

	QMatrix() = delete;
	QMatrix(size_t n_);
	~QMatrix();
	QMatrix(const QMatrix& A);
	QMatrix(QMatrix&& A);
	QMatrix(double** mem, size_t n_);
	QMatrix(double*** mem, size_t n_);
	QMatrix operator+(const QMatrix& R);
	QMatrix& operator+=(const QMatrix& R);
	QMatrix operator-(const QMatrix& R);
	QMatrix& operator-=(const QMatrix& R);
	QMatrix operator*(const QMatrix& R);
	QMatrix& operator=(const QMatrix& R);
	QMatrix& operator*=(const QMatrix& R);
	QMatrix pow(size_t k);
	QMatrix& spow(size_t k);
	QMatrix transp();
	QMatrix& stransp();
	double det(double epsilon = 1e-5);
	double cubic_norm();
	bool diagonal_dominance();
	double& at(size_t i, size_t j);
	void zero();
	void entity();
	size_t size();
	QMatrix& swap(size_t p, size_t k);
	QMatrix& tr_swap(size_t p, size_t k);
	void forget_transp();
	size_t main_element_tr(size_t k, double epsilon = 1e-5);
	QVector& vector_renumerate(QVector& vec);
	QMatrix direct_gauss_move(double epsilon = 1e-5);
	QMatrix& revitalize();
	QMatrix& not_forget_revitalize();
	QMatrix inverse(double epsilon = 1e-5);
	QMatrix& sinverse(double epsilon = 1e-5);
	QVectorT<size_t> PLU_decomposition(QMatrix& L_ret, QMatrix& U_ret, double epsilon = 1e-5);
	void setTransposition(QVectorT<size_t>& transp);
	QVectorT<size_t> getTransposition();
	QMatrix getTranspositionMatrix();
	size_t getSize();
	string serialize();
	friend ostream& operator<<(ostream& out, const QMatrix& A)
	{
		size_t n = A.n;
		for_dim(i)
		{
			for_dim(j)
			{
				out << setw(10) << left << A.src[i][j] << "  \t";
			}
			out << '\n';
		}
		return out;
	}
	friend istream& operator>> (istream& in, QMatrix& A)
	{
		size_t n = A.n;
		for_dim(i)
		{
			for_dim(j)
			{
				in >> A.src[i][j];
			}
		}

		return in;
	}
	friend QMatrix  operator*(double t, QMatrix& R)
	{
		size_t n = R.n;
		QMatrix P(n);

		for_dim(i)
			for_dim(j)
			P.src[i][j] = t*R.src[i][j];

		return P;
	}
	friend QMatrix  operator*(QMatrix& L, double t)
	{
		size_t n = L.n;
		QMatrix P(n);

		for_dim(i)
			for_dim(j)
			P.src[i][j] = t*L.src[i][j];

		return P;
	}
private:
	inline size_t count_zero(size_t diag, double** A, size_t fix, double epsilon = 1e-5);
	inline void setZero(double** R);
	inline void setEdentity(double** R);
	inline void generate_Pmn(int m_, int n_, double** out);
	inline void prod(double** A, double** B, double** out);
	inline void swap_mn(size_t m_, size_t n_, double** inout);
	void antidegenerating(double** A, double** P, size_t* p_ind, double epsilon = 1e-5);
	void copy_mem(double** dest, double** source);
	double** getMem();
	void freeMem(double** mem);
	inline void move_data(QMatrix&& R);
#undef for_dim
};

QVector operator*(QMatrix& R, QVector& v);