#include "QMatrix.h"
using namespace std;

#define for_dim(a) for(size_t (a) = 0; (a) < n; (a)++)


QMatrix::QMatrix(size_t n_) : n(n_), transposition(n_)
{
	src = getMem();
	for_dim(i)
	{
		transposition.at(i) = i;
	}
}
QMatrix::~QMatrix()
{
	freeMem(src);
}
QMatrix::QMatrix(const QMatrix& A) : QMatrix(A.n)
{
	copy_mem(src, A.src);
	transposition = A.transposition;
	transp_am = A.transp_am;
}
QMatrix::QMatrix(QMatrix&& A) : QMatrix(A.n)
{
	src = A.src;
	A.src = nullptr;
	transposition = A.transposition;
	transp_am = A.transp_am;
}
QMatrix::QMatrix(double** mem, size_t n_) : QMatrix(n_)
{
	copy_mem(src, mem);
}
QMatrix::QMatrix(double*** mem, size_t n_) : n(n_), transposition(n_)
{
	src = *mem;
	*mem = nullptr;
	for_dim(i)
	{
		transposition.at(i) = i;
	}
}
QMatrix QMatrix::operator+(const QMatrix& R)
{
	if (R.n != n)
	{
		throw error(errcode::DIFFERENCE_SIZE);
	}
	QMatrix S(n);
	for_dim(i)
		for_dim(j)
	{
		S.src[i][j] = src[i][j] + R.src[i][j];
	}
	return S;
}
QMatrix& QMatrix::operator+=(const QMatrix& R)
{
	if (R.n != n)
	{
		throw error(errcode::DIFFERENCE_SIZE);
	}
	for_dim(i)
		for_dim(j)
	{
		src[i][j] += R.src[i][j];
	}
	return *this;
}
QMatrix QMatrix::operator-(const QMatrix& R)
{
	if (R.n != n)
	{
		throw error(errcode::DIFFERENCE_SIZE);
	}
	QMatrix S(n);
	for_dim(i)
		for_dim(j)
	{
		S.src[i][j] = src[i][j] - R.src[i][j];
	}
	return S;
}
QMatrix& QMatrix::operator-=(const QMatrix& R)
{
	if (R.n != n)
	{
		throw error(errcode::DIFFERENCE_SIZE);
	}
	for_dim(i)
		for_dim(j)
	{
		src[i][j] -= R.src[i][j];
	}
	return *this;
}
void QMatrix::zero()
{
	setZero(src);
}
void QMatrix::entity()
{
	setEdentity(src);
}
QMatrix QMatrix::operator*(const QMatrix& R)
{
	if (R.n != n)
	{
		throw error(errcode::DIFFERENCE_SIZE);
	}
	QMatrix P(n);
	P.zero();

	for_dim(i)
		for_dim(j)
		for_dim(k)
		P.src[i][j] += src[i][k] * R.src[k][j];
	return P;
}
QMatrix& QMatrix::operator=(const QMatrix& R)
{
	if (R.n != n)
	{
		throw error(errcode::DIFFERENCE_SIZE);
	}
	copy_mem(src, R.src);
	transposition = R.transposition;
	transp_am = R.transp_am;
	return *this;
}
QMatrix& QMatrix::operator*=(const QMatrix& R)
{
	if (R.n != n)
	{
		throw error(errcode::DIFFERENCE_SIZE);
	}
	QMatrix P(n);
	P.zero();

	for_dim(i)
		for_dim(j)
		for_dim(k)
		P.src[i][j] += src[i][k] * R.src[k][j];

	move_data(move(P));;

	return *this;
}
QMatrix QMatrix::pow(size_t k)
{
	QMatrix P(n);
	P.entity();
	for (size_t i = 0; i < k; i++)
	{
		P *= (*this);
	}

	return P;
}
QMatrix& QMatrix::spow(size_t k)
{
	move_data(move(this->pow(k)));
	return *(this);
}
QMatrix QMatrix::transp()
{
	QMatrix T(n);
	for_dim(i)
		for (size_t j = i + 1; j < n; j++)
		{
			T.src[i][j] = src[j][i];
			T.src[j][i] = src[i][j];
		}

	for_dim(d)
	{
		T.src[d][d] = src[d][d];
	}

	return T;
}
QMatrix& QMatrix::stransp()
{
	for_dim(i)
		for (size_t j = i + 1; j < n; j++)
		{
			double t = src[i][j];
			src[i][j] = src[j][i];
			src[j][i] = t;
		}

	return *(this);
}
double QMatrix::det(double epsilon)
{
	QMatrix Self = *this;
	Self.forget_transp();

	try
	{
		Self.direct_gauss_move(epsilon);
	}
	catch (error e)
	{
		if (e.code == errcode::ZERO_ON_DIAG || e.code == errcode::DEGENERATING_MATRIX) return 0.;
		else throw e;
	}

	double determinate = 1;
	if (Self.transp_am % 2 == 1) determinate = -1;

	for (size_t i = 0; i < Self.n; i++)
	{
		determinate *= Self.src[i][i];
	}
	return determinate;
}
double QMatrix::cubic_norm()
{
	double old_acc = 0;
	for_dim(i)
	{
		double acc = 0;
		for_dim(k) acc += abs(src[i][k]);
		if (acc > old_acc) old_acc = acc;
	}

	return old_acc;
}
bool QMatrix::diagonal_dominance()
{
	for_dim(i)
	{
		double summ = 0;
		for_dim(j)
		{
			if (i == j) continue;
			summ += abs(src[i][j]);
		}

		if (abs(src[i][i]) <= summ) return false;
	}
	return true;
}
double& QMatrix::at(size_t i, size_t j)
{
	if (i >= n || j >= n)
	{
		throw error(errcode::BAD_INDEX);
	}
	return src[i][j];
}
size_t QMatrix::size()
{
	return n;
}
QMatrix& QMatrix::swap(size_t p, size_t k)
{
	if (p == k) return *this;

	double* taddr;
	taddr = src[p];
	src[p] = src[k];
	src[k] = taddr;

	return *this;
}
QMatrix& QMatrix::tr_swap(size_t p, size_t k)
{
	if (p == k) return *this;

	this->swap(p, k);
	++transp_am;
	transposition.swap(p, k);

	return *this;
}
void QMatrix::forget_transp()
{
	for_dim(i) transposition.at(i) = i;
	transp_am = 0;
}
size_t QMatrix::main_element_tr(size_t k, double epsilon)
{
	if (k >= n)
	{
		throw error(errcode::BAD_INDEX);
	}

	auto eq = [k](double* a, double* b) -> bool {return abs(a[k]) < abs(b[k]); };

	double** main_line = max_element(src + k, src + n, eq);

	if (abs(main_line[0][k]) < epsilon) throw error(errcode::DEGENERATING_MATRIX);

	double* taddr = src[k];
	src[k] = *main_line;
	*main_line = taddr;

	if (k != main_line - src)
	{
		transposition.swap(k, main_line - src);
		++transp_am;
	}

	return main_line - src;
}
QVector& QMatrix::vector_renumerate(QVector& vec)
{
	if (vec.size() != n)
	{
		throw error(errcode::DIFFERENCE_SIZE);
	}
	QVector old = vec;
	for_dim(k)
	{
		vec.at(k) = old.at(transposition.at(k));
	}

	return vec;
}
QMatrix QMatrix::direct_gauss_move(double epsilon)
{
	QMatrix& A = *this;
	QMatrix Right(n);
	Right.entity();

	for_dim(diag)
	{
		size_t m = A.main_element_tr(diag);
		Right.tr_swap(diag, m);

		double diagVal = A.at(diag, diag);

		if (abs(diagVal) < epsilon) throw error(errcode::ZERO_ON_DIAG);

		for (size_t i = diag + 1; i < n; i++)
		{
			double lineCoeff = A.at(i, diag) / diagVal;
			for (int j = 0; j < n; j++)
			{
				A.at(i, j) -= lineCoeff*A.at(diag, j);
				Right.at(i, j) -= lineCoeff*Right.at(diag, j);
			}
		}
	}

	return move(Right);
}
QMatrix& QMatrix::revitalize()
{
	double** src_old = src;

	src = new double*[n];
	for_dim(i)
	{
		src[i] = src_old[transposition.at(i)];
	}
	delete[] src_old;

	for_dim(i) transposition.at(i) = i;
	transp_am = 0;

	return *this;
}
QMatrix& QMatrix::not_forget_revitalize()
{
	double** src_old = src;

	src = new double*[n];
	for_dim(i)
	{
		src[i] = src_old[transposition.at(i)];
	}
	delete[] src_old;

	return *this;
}
QMatrix QMatrix::inverse(double epsilon)
{
	QMatrix A(*this);

	QMatrix Right = A.direct_gauss_move(epsilon);

	for (int diag = n - 1; diag >= 0; diag--)
	{
		double diagVal = A.at(diag, diag);

		for_dim(i)
		{
			A.at(diag, i) /= diagVal;
			Right.at(diag, i) /= diagVal;
		}

		for (int i = diag - 1; i >= 0; i--)
		{
			double lineCoeff = A.at(i, diag);
			for (int j = n - 1; j >= 0; j--)
			{
				Right.at(i, j) -= lineCoeff*Right.at(diag, j);
			}
		}
	}

	return Right;
}
QMatrix& QMatrix::sinverse(double epsilon)
{
	move_data(move(this->inverse(epsilon)));
	return *this;
}
QVectorT<size_t> QMatrix::PLU_decomposition(QMatrix& L_ret, QMatrix& U_ret, double epsilon)
{
	U_ret = *this;
	U_ret.forget_transp();
	L_ret = U_ret.direct_gauss_move(epsilon);
	L_ret.sinverse(epsilon);
	L_ret.not_forget_revitalize();

	U_ret.forget_transp();

	return L_ret.transposition;
}
void QMatrix::setTransposition(QVectorT<size_t>& transp)
{
	transposition = transp;
}
QVectorT<size_t> QMatrix::getTransposition()
{
	return transposition;
}
QMatrix QMatrix::getTranspositionMatrix()
{
	QMatrix P(n);
	P.zero();
	for_dim(i)
	{
		P.at(i, transposition.at(i)) = 1;
	}

	return P;
}
size_t QMatrix::getSize()
{
	return n;
}
string QMatrix::serialize()
{
	stringstream converter;
	string accum = "[";
	for_dim(i)
	{
		accum += "[";
		for_dim(j)
		{
			string elem;
			converter << src[i][j];
			converter >> elem;
			if (j != n - 1) accum += elem + ", ";
			else accum += elem + "]";
			converter.str("");
			converter.clear();
		}
	}
	accum += "]";
	return accum;
}
size_t QMatrix::count_zero(size_t diag, double** A, size_t fix, double epsilon)
{
	int cnt = 0;
	for (size_t i = diag; i < n; i++)
	{
		if (abs(A[fix][i]) < epsilon) cnt++;
	}

	return cnt;
}
void QMatrix::setZero(double** R)
{
	for_dim(i)
	{
		memset(R[i], 0, sizeof(double)*n);
	}
}
void QMatrix::setEdentity(double** R)
{
	setZero(R);

	for_dim(i)
		R[i][i] = 1.;
}
void QMatrix::generate_Pmn(int m_, int n_, double** out)
{
	setEdentity(out);

	out[m_][m_] = 0.0;
	out[m_][n_] = 1.0;
	out[n_][n_] = 0.0;
	out[n_][m_] = 1.0;
}
void QMatrix::prod(double** A, double** B, double** out)
{
	setZero(out);

	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			for (size_t k = 0; k < n; k++)
			{
				out[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}
void QMatrix::swap_mn(size_t m_, size_t n_, double** inout)
{
	double* taddr;
	taddr = inout[m_];
	inout[m_] = inout[n_];
	inout[n_] = taddr;
}
void QMatrix::antidegenerating(double** A, double** P, size_t* p_ind, double epsilon)
{
	size_t c = 0;
	for_dim(diag)
	{
		if (abs(A[diag][diag]) < epsilon) c++;
	}

	setEdentity(P);

	for (int i = 0; i < n; i++) p_ind[i] = i;

	if (c == 0) return;

	double** P_current = getMem();
	double** temp = getMem();

	for_dim(diag)
	{
		int c_zero = -1;
		size_t k = 0;

		for (int i = diag; i < n; i++)
		{
			if (abs(A[i][diag]) > epsilon)
			{
				int cz_current = count_zero(diag, A, i, epsilon);
				if (cz_current > c_zero)
				{
					c_zero = cz_current;
					k = i;
				}
			}
		}
		if (k == diag)
		{
			continue;
		}
		if (c_zero == -1)
		{
			throw error(errcode::DEGENERATING_MATRIX);
		}

		generate_Pmn(diag, k, P_current);
		prod(P_current, P, temp);
		copy_mem(P, temp);

		swap_mn(diag, k, A);

		size_t t = p_ind[diag];
		p_ind[diag] = p_ind[k];
		p_ind[k] = t;
	}

	freeMem(P_current);
	freeMem(temp);
}
void QMatrix::copy_mem(double** dest, double** source)
{
	for_dim(i)
	{
		memcpy(dest[i], source[i], sizeof(double)*n);
	}
}
double** QMatrix::getMem()
{
	double** mem = new double*[n];
	for_dim(i)
	{
		mem[i] = new double[n];
	}
	return mem;
}
void QMatrix::freeMem(double** mem)
{
	if (mem == nullptr) return;
	for_dim(i)
	{
		delete[] mem[i];
	}
	delete[] mem;
}
void QMatrix::move_data(QMatrix&& R)
{
	double** data = R.src;
	R.src = nullptr;
	freeMem(src);
	this->src = data;
}
QVector operator*(QMatrix& R, QVector& v)
{
	size_t n = v.size();
	QVector p(n);
	p.zero();

	for_dim(i)
		for_dim(k)
		p.at(i) += R.at(i, k)*v.at(k);

	return p;
#undef for_dim
}