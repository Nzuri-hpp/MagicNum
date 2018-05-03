#pragma once
#include "Polynomial.h"

namespace integration
{

	enum QUADRATURE
	{
		HOMOGEN,
		VARIOUS
	};

	enum CONSTRUCTOR
	{
		HERMITE,
		LAGGER,
		__ALG
	};

	template<int TYPE>
	class Quadrature
	{
		typedef std::vector<double> vec_d;
		typedef double(*target_func)(double);

		const vec_d & x_nodes;
		const vec_d & weighting_factors;

		double (Quadrature::*convolution_alg)(const double*, size_t);
		double*(Quadrature::*computing_alg)(target_func);
	public:
		struct error
		{
			const size_t code;
			error(size_t code_) : code(code_) {}
		};
		struct errcode
		{
			enum codes { DIFFERENCE_SIZE };
		};
		Quadrature(const vec_d & _x_nodes, const vec_d & _weighting_factors)
			: x_nodes(_x_nodes), weighting_factors(_weighting_factors)
		{
			convolution_alg = &Quadrature<QUADRATURE::VARIOUS>::convolution_simple;
			computing_alg   = &Quadrature<QUADRATURE::VARIOUS>::computing_values_simple;
			if (x_nodes.size() != weighting_factors.size())
				throw error(errcode::DIFFERENCE_SIZE);
		}
		Quadrature(const std::pair<vec_d, vec_d> & kernel)
			: x_nodes(std::get<0>(kernel)), weighting_factors(std::get<1>(kernel))
		{
			convolution_alg = &Quadrature<QUADRATURE::VARIOUS>::convolution_simple;
			computing_alg   = &Quadrature<QUADRATURE::VARIOUS>::computing_values_simple;
			if (x_nodes.size() != weighting_factors.size())
				throw error(errcode::DIFFERENCE_SIZE);
		}
		double operator()(target_func F)
		{
			auto values = (this->*computing_alg)(F);
			return convolution_simple(values, x_nodes.size());
		}

	private:
		double* computing_values_simple(target_func F)
		{
			size_t L = x_nodes.size();
			const auto x_nodes_raw = x_nodes.data();

			double *values = new double[L];
			for (size_t i = 0; i < L; ++i)
			{
				values[i] = F(x_nodes_raw[i]);
			}

			return values;
		}
		double convolution_simple(const double* values, size_t L)
		{
			const auto weighting_factors_raw = weighting_factors.data();

			double res = 0;
			for (size_t i = 0; i < L; ++i)
			{
				res += values[i] * weighting_factors_raw[i];
			}

			return res;
		}
	};

	template<>
	class Quadrature<QUADRATURE::HOMOGEN>
	{
		typedef std::vector<double> vec_d;
		typedef double(*target_func)(double);

		const vec_d & x_nodes;
		const double weighting_factor;

		double (Quadrature<QUADRATURE::HOMOGEN>::*convolution_alg)(const double*, size_t);
		double*(Quadrature<QUADRATURE::HOMOGEN>::*computing_alg)(target_func);
	public:
		Quadrature(const vec_d & _x_nodes, double _weighting_factor)
			: x_nodes(_x_nodes), weighting_factor(_weighting_factor)
		{
			convolution_alg = &Quadrature<QUADRATURE::HOMOGEN>::convolution_simple;
			computing_alg   = &Quadrature<QUADRATURE::HOMOGEN>::computing_values_simple;
		}
		Quadrature(const std::pair<vec_d, double> & kernel)
			: x_nodes(std::get<0>(kernel)), weighting_factor(std::get<1>(kernel))
		{
			convolution_alg = &Quadrature<QUADRATURE::HOMOGEN>::convolution_simple;
			computing_alg   = &Quadrature<QUADRATURE::HOMOGEN>::computing_values_simple;
		}
		double operator()(target_func F)
		{
			auto values = (this->*computing_alg)(F);
			return convolution_simple(values, x_nodes.size());
		}

	private:
		double* computing_values_simple(target_func F)
		{
			size_t L = x_nodes.size();
			double *values = new double[L];
			for (size_t i = 0; i < L; ++i)
			{
				values[i] = F(x_nodes[i]);
			}

			return values;
		}
		double convolution_simple(const double* values, size_t L)
		{
			double res = 0;
			for (size_t i = 0; i < L; ++i)
			{
				res += values[i];
			}

			return weighting_factor * res;
		}
	};

	typedef Quadrature<QUADRATURE::HOMOGEN> Quadrature_H;
	typedef Quadrature<QUADRATURE::VARIOUS> Quadrature_V;

	template<int TYPE>
	class Constructor
	{

	};

	template<>
	class Constructor<CONSTRUCTOR::__ALG>
	{
		std::vector<double> factorial_buffer;
		std::vector<double> sum_log_buffer;
	public:
		Constructor() 
		{
			factorial_buffer.reserve(100);
			factorial_buffer.push_back(1.);

			sum_log_buffer.reserve(100);
			sum_log_buffer.push_back(0.);
		}

		double stirling(double n)
		{
			return sqrt(2.*M_PI*n)*pow(n / M_E, n);
		}

		double factorial(double n)
		{
			if (n < 1) return 1.;

			double k = 1;
			double fl = 1;

			if (n <= factorial_buffer.size())
			{
				return factorial_buffer[static_cast<size_t>(n) - 1];
			}
			if (n > factorial_buffer.size())
			{
				k  = static_cast<double>(factorial_buffer.size()) + 1.;
				fl = factorial_buffer.back();
			}

			for (double i = k; i <= n; i += 1.)
			{
				fl *= i;
				factorial_buffer.push_back(fl);
			}

			return fl;
		}

		double sum_log(double n)
		{
			if (n < 1) return 0.;

			double k = 1;
			double fl = 0;

			if (n <= sum_log_buffer.size())
			{
				return sum_log_buffer[static_cast<size_t>(n) - 1];
			}
			if (n > sum_log_buffer.size())
			{
				k = static_cast<double>(sum_log_buffer.size()) + 1.;
				fl = sum_log_buffer.back();
			}

			for (double i = k; i <= n; i += 1.)
			{
				fl += std::log(i);
				sum_log_buffer.push_back(fl);
			}

			return fl;
		}

		double sum_log(double shift, double n)
		{
			return sum_log(n + shift) - sum_log(shift);
		}

		inline void clear_buffers()
		{
			factorial_buffer.clear();
			sum_log_buffer.clear();
		}

		inline void shrink_buffers()
		{
			factorial_buffer.shrink_to_fit();
			sum_log_buffer.shrink_to_fit();
		}

		template<size_t N> 
		struct ctime_sum_log
		{
			const double value = std::log(N) + sum_log<N - 1>().value;

			inline static double get()
			{
				return ctime_sum_log().value;
			}
		};

		template<>
		struct ctime_sum_log<1>
		{
			const double value = 0;

			inline static double get()
			{
				return ctime_sum_log().value;
			}
		};

		Polynomial lagger(size_t n_size)
		{
			std::vector<double> L_coeffs(n_size + 1);

			size_t k = 0;
			for (auto& a : L_coeffs)
			{
				size_t s = std::min(k, n_size - k);
				size_t o = std::max(k, n_size - k);

				double p1 = sum_log(o, s);
				double p2 = sum_log(s);
				double p3 = sum_log(k) - p2;

				double pi_k = p1 - 2 * p2 - p3;

				double mult = k % 2 ? 1. : -1.;

				L_coeffs[k] = mult * std::exp(pi_k);

				++k;
			}

			return Polynomial(move(L_coeffs));
		}

		std::vector<double> graeffes_iteration(std::vector<double> & p_coeff)
		{
			std::vector<double> p_next(p_coeff.size());

			int n_size = p_coeff.size();
			auto i = [=](int k)-> int {return n_size - k - 1;};

			p_next[i(0)] = p_coeff[i(0)];
			for (int k = 1; k < n_size; ++k)
			{
				p_next[i(k)] = (k % 2 ? 1. : -1.)*p_coeff[i(k)] * p_coeff[i(k)];

				double acc = 0;
				for (int j = 0; j <= k - 1; ++j)
				{
					if (2 * k - j >= n_size) continue;
					acc += (j % 2 ? 1. : -1.) * p_coeff[i(j)] * p_coeff[i(2*k - j)];
				}

				p_next[i(k)] += 2 * acc;

				p_next[i(k)] /= p_next[i(0)];
			}

			return move(p_next);
		}

		std::vector<double> graeffes_roots(Polynomial p, double eps = DBL_EPSILON, size_t max_step = INT_MAX)
		{
			auto criterion = [](const vector<double> & pc_next, const vector<double> & pc_curr)->double
			{
				double max_val = 0;

				for (size_t i = 0; i < pc_next.size(); ++i)
				{
					double val = std::abs(pc_next[i] / (pc_curr[i] * pc_curr[i]));
					if (std::abs(val - 1.) > max_val)
						max_val = std::abs(val - 1.);
				}

				return max_val;
			};

			auto abs_min_coeff = [](const vector<double> & pc_next) -> double
			{
				return abs(*std::min_element(pc_next.begin(), pc_next.end(), [](double a, double b) {return abs(a) < abs(b); }));
			};

			auto abs_max_coeff = [](const vector<double> & pc_next) -> double
			{
				return abs(*std::max_element(pc_next.begin(), pc_next.end(), [](double a, double b) {return abs(a) < abs(b); }));
			};

			vector<double> curr_coeff = p.getCoeff();
			vector<double> next_coeff = graeffes_iteration(curr_coeff);

			size_t step = 1;
			bool aa = criterion(next_coeff, curr_coeff) > eps;
			bool bb = abs_min_coeff(next_coeff) > 1e-90;
			while (criterion(next_coeff, curr_coeff) > eps && step < max_step && abs_min_coeff(next_coeff) > 1e-50 && abs_max_coeff(next_coeff) < 1e+50)
			{
				curr_coeff = move(next_coeff);
				next_coeff = graeffes_iteration(curr_coeff);
				++step;
			}

			vector<double> x_roots;
			x_roots.reserve(p.degree());
			
			auto i = [&](int k)-> int {return p.degree() - k; };

			for (size_t k = 1; k <= p.degree(); ++k)
			{
				double buff = abs(next_coeff[i(k)] / next_coeff[i(k - 1)]);
				x_roots.push_back( std::pow(buff, 1. / std::pow(2., step)) );
			}

			std::sort(x_roots.begin(), x_roots.end());
			return move(x_roots);
		}

		Polynomial radial_normalization(Polynomial & p, double radius)
		{
			auto cfs = p.getCoeff();
			size_t n = cfs.size() - 1;
			for (size_t k = 1; k < n; ++k)
			{
				cfs[n - k] *= factorial(n - k) * factorial(k) / (factorial(n) * std::pow(radius, k));
			}

			return Polynomial(move(cfs));
		}
	};

	template<>
	class Constructor<CONSTRUCTOR::HERMITE>
	{
		typedef std::map<size_t, std::pair<std::vector<double>, double> > cache_t;
		cache_t cache;
	public:
		Quadrature_H construct(size_t n_size)
		{
			double n_size_d = static_cast<double>(n_size);
			auto p = cache.emplace(n_size, 
								   make_pair(move(std::vector<double>()), 
								   M_PI/n_size_d ) );
			auto it = std::get<0>(p);

			if (!std::get<1>(p))
			{
				return Quadrature_H(it->second);
			}

			auto& x_points = std::get<0>(it->second);

			x_points.reserve(n_size);

			double mult_const = M_PI / (2.*n_size_d);
			for (double i = 0; i < n_size_d; i+=1.)
			{
				double xi = std::cos((2 * i - 1)*mult_const);
				x_points.push_back(xi);
			}

			return Quadrature_H(it->second);
		}
	};

	template<>
	class Constructor<CONSTRUCTOR::LAGGER> 
		: private Constructor<CONSTRUCTOR::__ALG>
	{
		typedef std::map<size_t, std::pair<std::vector<double>, std::vector<double> > > cache_t;
		cache_t cache;
	public:
		Quadrature_V construct(size_t n_size)
		{
			double n_size_d = static_cast<double>(n_size);
			auto p = cache.emplace(n_size,
								   make_pair(move(std::vector<double>()),
								   M_PI / n_size_d));
			auto it = std::get<0>(p);

			if (!std::get<1>(p))
			{
				return Quadrature_V(it->second);
			}
			
			auto& x_points = std::get<0>(it->second);
			auto& a_coeffs = std::get<1>(it->second);

			Polynomial Lagger = lagger(n_size);
			Lagger.take_normal();

			double scale = 1;
			auto fake_x = graeffes_roots(Lagger);
			double max_fake_x = *std::max_element(fake_x.begin(), fake_x.end());

			x_points = move(Lagger.pricisie_all_roots(0, scale * max_fake_x));
			Lagger.sturms_method(0, scale * max_fake_x);
			if (x_points.size() != n_size)
			{
				double __INF = std::numeric_limits<double>::infinity();
				while (Lagger.sturms_method(scale * max_fake_x, __INF) != 0)
				{
					scale += 1.;
				}

				x_points = move(Lagger.pricisie_all_roots(0, scale * max_fake_x));
			}

			n_size   = x_points.size() < n_size ? x_points.size() : n_size;
			n_size_d = static_cast<double>(n_size);

			a_coeffs.resize(n_size);

			Polynomial l_poly = lagger(n_size + 1);

			double const_coeff = 1. / ((n_size + 1)*(n_size + 1));

			for (size_t i = 0; i < n_size; ++i)
			{
				double li = l_poly.horners_method(x_points[i]);
				a_coeffs[i] = const_coeff*x_points[i]/(li*li);
			}

			return Quadrature_V(it->second);
		}

		static size_t precision_test(size_t lim, double eps = 0.001)
		{
			Constructor<integration::CONSTRUCTOR::LAGGER> master_L;
			auto c = [](double x) {return 1.; };
			for (size_t i = 1; i <= lim; ++i)
			{
				auto qr = master_L.construct(i);
				if (std::abs(qr(c) - 1) > eps) return i - 1;
			}

			return lim;
		}
	};

}
