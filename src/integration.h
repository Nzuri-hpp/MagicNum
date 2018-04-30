#pragma once

namespace integration
{

	enum QUADRATURE
	{
		HOMOGEN,
		VARIOUS
	};

	enum CONSTRUCTOR
	{
		HERMITE
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

			return it->second;
		}
	};

}
