#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#define for_dim(a) for(size_t (a) = 0; (a) < n; (a)++)

using namespace std;
template<class T>
class QVectorT
{
	T *src;
	const size_t n;
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
	QVectorT(size_t n_);
	~QVectorT();
	QVectorT(const QVectorT& v);
	QVectorT(double* vsrc, size_t n_);
	QVectorT(const std::vector<T>& vec);
	QVectorT operator+(QVectorT& r);
	QVectorT& operator=(const QVectorT& r);
	QVectorT operator-(QVectorT& r);
	QVectorT& operator+=(QVectorT& r);
	QVectorT& operator-=(QVectorT& r);
	T cubic_norm();
	T euqlid_norm();
	void reindex(size_t* const p_ind);
	QVectorT operator*(double t);
	T& at(size_t i);
	size_t size();
	void zero();
	QVectorT& operator=(QVectorT&& r);
	QVectorT& swap(size_t p, size_t k);
	string serialize();
	friend QVectorT<T> operator*(T t, QVectorT<T>& r) 
	{
		size_t n = r.n;
		QVectorT<T> p(n);

		for_dim(i) p.src[i] = t*r.src[i];

		return p;
	}
	friend void swap(QVectorT<T>& l, QVectorT<T>& r) 
	{
		if (l.size() != r.size())
		{
			throw QVectorT<T>::error(QVectorT<T>::errcode::DIFFERENCE_SIZE);
		}
		T* taddr = l.src;
		l.src = r.src;
		r.src = taddr;
	}
	friend ostream& operator<<(ostream& out, const QVectorT<T>& v)
	{
		size_t n = v.n;

		for_dim(i)
		{
			out << setw(15) << left << v.src[i] << endl;
		}

		return out;
	}
	friend istream& operator>>(istream& in, QVectorT<T>& v)
	{
		size_t n = v.n;
		for_dim(i)
		{
			in >> v.src[i];
		}

		return in;
	}
private:
	void setZero(double* s);
	void copy_mem(T* dest, T* source);
};
template<class T>
QVectorT<T>::QVectorT(size_t n_) : n(n_)
{
	src = new T[n];
}
template<class T>
QVectorT<T>::~QVectorT()
{
	delete[] src;
}
template<class T>
QVectorT<T>::QVectorT(const QVectorT& v) : QVectorT(v.n)
{
	copy_mem(src, v.src);
}
template<class T>
QVectorT<T>::QVectorT(double* vsrc, size_t n_) : QVectorT(n_)
{
	copy_mem(src, vsrc);
}
template<class T>
QVectorT<T>::QVectorT(const std::vector<T>& vec) : QVectorT(vec.size())
{
	for_dim(i)
	{
		src[i] = vec.at(i);
	}
	
}
template<class T>
QVectorT<T> QVectorT<T>::operator+(QVectorT& r)
{
	if (r.n != n)
	{
		throw error(errcode::DIFFERENCE_SIZE);
	}

	QVectorT<T> s(n);

	for_dim(i) s.src[i] = src[i] + r.src[i];

	return s;
}
template<class T>
QVectorT<T>& QVectorT<T>::operator=(const QVectorT& r)
{
	if (r.n != n)
	{
		throw error(errcode::DIFFERENCE_SIZE);
	}
	copy_mem(src, r.src);

	return *this;
}
template<class T>
QVectorT<T> QVectorT<T>::operator-(QVectorT& r)
{
	if (r.n != n)
	{
		throw error(errcode::DIFFERENCE_SIZE);
	}

	QVectorT<T> s(n);

	for_dim(i) s.src[i] = src[i] - r.src[i];

	return s;
}
template<class T>
QVectorT<T>& QVectorT<T>::operator+=(QVectorT& r)
{
	if (r.n != n)
	{
		throw error(errcode::DIFFERENCE_SIZE);
	}

	for_dim(i) src[i] += r.src[i];

	return *this;
}
template<class T>
QVectorT<T>& QVectorT<T>::operator-=(QVectorT& r)
{
	if (r.n != n)
	{
		throw error(errcode::DIFFERENCE_SIZE);
	}

	for_dim(i) src[i] -= r.src[i];

	return *this;
}
template<class T>
T QVectorT<T>::cubic_norm()
{
	T maxi = 0;
	for_dim(i)
	{
		if (abs(src[i]) > maxi) maxi = abs(src[i]);
	}

	return maxi;
}
template<class T>
T QVectorT<T>::euqlid_norm()
{
	T sum = 0;
	for_dim(i)
	{
		sum += src[i] * src[i];
	}

	return sqrt(sum);
}
template<class T>
void QVectorT<T>::reindex(size_t* const p_ind)
{
	T* temp = new T[n];
	copy_mem(temp, src);

	for_dim(i)
	{
		src[i] = temp[p_ind[i]];
	}

	delete[] temp;
}
template<class T>
QVectorT<T> QVectorT<T>::operator*(double t)
{
	QVectorT<T> p(n);

	for_dim(i) p.src[i] = t*src[i];

	return p;
}
template<class T>
T& QVectorT<T>::at(size_t i)
{
	if (i >= n)
	{
		throw error(errcode::BAD_INDEX);
	}
	return src[i];
}
template<class T>
size_t QVectorT<T>::size()
{
	return n;
}
template<class T>
void QVectorT<T>::zero()
{
	setZero(src);
}
template<class T>
QVectorT<T>& QVectorT<T>::operator=(QVectorT<T>&& r)
{
	copy_mem(src, r.src);

	return *this;
}
template<class T>
QVectorT<T>& QVectorT<T>::swap(size_t p, size_t k)
{
	T t = src[p];
	src[p] = src[k];
	src[k] = t;

	return *this;
}
template<class T>
string QVectorT<T>::serialize()
{
	stringstream converter;
	string accum = "[";
	for_dim(i)
	{
		string elem;
		converter << src[i];
		converter >> elem;
		if(i != n-1) accum += elem + ", ";
		else accum += elem + "]";
		converter.str("");
		converter.clear();
	}
	return accum;
}
template<class T>
void QVectorT<T>::setZero(double* s)
{
	memset(s, (T)0, sizeof(T)*n);
}
template<class T>
void QVectorT<T>::copy_mem(T* dest, T* source)
{
	memcpy(src, source, sizeof(T)*n);
}
#undef for_dim

typedef QVectorT<double> QVector;