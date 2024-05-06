#include <iostream>
#include <array>
#include <cmath>
const int n = 50;
double alpha = 4.0;
double beta = 1.0;
double gamma = 2.0;
double h = 1.0 / double(n);
double EPS = h * h * h;

double p(double x)
{
	return 1 + std::pow(x, gamma);
}
double p_(double x)
{
	return gamma * (1, std::pow(x, gamma - 1));
}
double q(double x)
{
	return x + 1;
}
double u(double x)
{
	return std::pow(x, alpha) * std::pow(1 - x, beta);
}
double u_(double x)
{
	return alpha * beta
		* (std::pow(x, alpha - 1))
		* (std::pow(1 - x, beta - 1))
		- beta * std::pow(x, alpha)
		* (std::pow(1 - x, beta - 1));
}
double u__(double x)
{
	return alpha
		* (alpha - 1)
		* (std::pow(x, alpha - 2))
		* std::pow(1 - x, beta)
		- 2 * alpha * beta
		* (std::pow(x, alpha - 1))
		* (std::pow(1 - x, beta - 1))
		+ beta * (beta - 1) * std::pow(x, alpha)
		* (std::pow(1 - x, beta - 2));
}
double f(double x)
{
	return -(p_(x) * u_(x) + p(x) * u__(x)) + q(x) * u(x);
}
double f(int i)
{
	return f(double(i) * h);
}
double a(int i)
{
	return p(double(i) * h);
}
double g(int i)
{
	return q(double(i) * h);
}
class Vector
{
	std::array<double, n> v;
public:
	Vector() { for (int i = 0; i < n; i++) v[i] = 0.0; }
	double& operator[](int i) { return v[i]; }
	const double& operator[](int i) const { return v[i]; }
	std::string printVert() const {
		std::string out;
		for (int i = 0; i < n; i++) {
			std::cout << v[i] << '\n'; out += v[i];
		}
		return out;
	}
	friend std::ostream& operator<<(std::ostream& os, const Vector& v)
	{
		for (int i = 0; i < n; i++)
			os << v[i] << ' ';
		//os << '\n';
		return os;
	}
};
class Matrix
{
	std::array<Vector, n> mtx;
public:
	Vector& operator[](int i) { return mtx[i]; }
	const Vector& operator[](int i) const { return mtx[i]; }
	friend std::ostream& operator<<(std::ostream& os, const Matrix& mtx)
	{
		os << "Matrix<" << n << ">:\n";
		for (int i = 0; i < n - 1; i++)
			os << mtx[i] << '\n';
		os << mtx[n - 1];
		return os;
	}
};
double operator*(const Vector& v1, const Vector& v2)
{
	double res = 0.0;
	for (int i = 0; i < n; i++)
		res += v1[i] * v2[i];
	return res;
}
Vector operator*(double val, Vector v)
{
	for (int i = 0; i < n; i++)
		v[i] = v[i] * val;
	return v;
}
Vector operator*(const Matrix& m, const Vector& v)
{
	Vector res;
	for (int i = 0; i < n; i++)
	{
		res[i] = m[i] * v;
	}
	return res;
}
Vector operator-(Vector v1, const Vector& v2)
{
	for (int i = 0; i < n; i++)
		v1[i] = v1[i] - v2[i];
	return v1;
}
Vector abs(Vector v)
{
	for (int i = 0; i < n; i++)
		v[i] = std::abs(v[i]);
	return v;
}
bool operator>(const Vector& v, double val)
{
	for (int i = 1; i < n - 1; i++)
		if (v[i] < val) return false;
	return true;
}

class Tridiagonal
{
protected:
	Matrix mtx;
	Vector y;
	Vector alphas;
	Vector betas;

public:

	void InitializeMatrix()
	{
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				mtx[i][j] = 0.0;
		double ai = 0.0;
		double ci = a(0) + a(1) + h * h * g(0);
		double bi = a(1);
		if (n > 2)
		{
			alphas[0] = bi / ci;
			betas[0] = 0.0;
			mtx[0][0] = ai;
			mtx[0][1] = ci;
			mtx[0][2] = bi;
		}

		for (int i = 1; i < n - 1; i++)
		{
			ai = alphas[i - 1] * a(i);
			ci = a(i) + a(i + 1) + h * h * g(i);
			bi = a(i + 1);

			mtx[i][i - 1] = ai;
			mtx[i][i] = ci;
			mtx[i][i + 1] = bi;

			alphas[i] = bi / (ci - ai);
			betas[i] = (f(i) * h * h - betas[i - 1] * a(i)) / (ci - ai);
		}

		ai = alphas[n - 2] * a(n - 1);
		ci = a(n - 1) + a(n) + h * h * g(n - 1);
		bi = 0.0;
		alphas[n - 1] = bi / (ci - ai);
		betas[n - 1] = (f(n - 1) * h * h - betas[n - 2] * a(n - 1)) / (ci - ai);
		mtx[n - 1][n - 3] = ai;
		mtx[n - 1][n - 2] = ci;
		mtx[n - 1][n - 1] = bi;
	}

	virtual const Matrix& getMatrix() const = 0;

	virtual const Vector& getAlphas() const = 0;
	virtual const Vector& getBettas() const = 0;

	virtual const Vector& getSolutions() = 0;
};
class Progonka : public Tridiagonal
{
public:
	Progonka()
	{
		InitializeMatrix();
	}

	const Matrix& getMatrix() const override
	{
		return mtx;
	};

	const Vector& getAlphas() const override
	{
		return alphas;
	}
	const Vector& getBettas() const override
	{
		return betas;
	}

	const Vector& getSolutions() override
	{
		y[n - 1] = 0.0;
		y[0] = 0.0;
		for (int i = n - 2; i > 0; i--) {
			y[i] = alphas[i] * y[i + 1] + betas[i];
		}
		return y;
	}
};

class Jacobi : Progonka
{
public:
	Jacobi()
	{
		InitializeMatrix();
	}

	const Matrix& getMatrix() const override
	{
		return mtx;
	}

	const Vector& getAlphas() const override
	{
		return alphas;
	}

	const Vector& getBettas() const override
	{
		return betas;
	}

	const Vector& getSolutions() override
	{
		Vector r;
		Vector temp;
		double norm = 0.0;
		int k = 0;
		do {
			for (int i = 1; i < n - 1; i++)
			{
				double ai = mtx[i][i - 1];
				double ci = mtx[i][i];
				double bi = mtx[i][i + 1];
				y[i] = (a(i) * temp[i - 1] + a(i + 1) * temp[i + 1]
					+ f(i) * h * h) / (mtx[i][i]);
			}
			r = temp - y;
			temp = y;
			k++;
		} while (abs(r) > EPS);
		std::cout << "Количество итераций: " << k << '\n';
		return y;
	}
};

class Zeidel : Progonka
{
public:
	Zeidel()
	{
		InitializeMatrix();
	}

	const Matrix& getMatrix() const override
	{
		return mtx;
	};

	const Vector& getAlphas() const override
	{
		return alphas;
	}
	const Vector& getBettas() const override
	{
		return betas;
	}

	const Vector& getSolutions() override
	{
		Vector r;
		double t = 0.0;
		int k = 0;
		do {
			for (int i = 1; i < n - 1; i++)
			{
				y[i] = y[i] - t * r[i];
				r[i] = -a(i) * y[i - 1] + mtx[i][i] * y[i] - a(i) * y[i + 1] - f(i) * h * h;
			}
			t = r * r / ((mtx * r) * r);

			k++;
		} while (abs(r) > EPS);
		std::cout << "Количество итераций: " << k << '\n';
		return y;
	}
};

class UpRelaxation : Progonka
{
public:
	UpRelaxation()
	{
		InitializeMatrix();
	}
	const Matrix& getMatrix() const override
	{
		return mtx;
	};

	const Vector& getAlphas() const override
	{
		return alphas;
	}
	const Vector& getBettas() const override
	{
		return betas;
	}

	const Vector& getSolutions() override
	{
		Vector best;
		int pre_k = INT32_MAX;
		double step = 1.0 / 10.0;
		for (double w = 1.0; w < 2.0 + step; w += step)
		{
			Vector temp;
			Vector r;
			int k = 0;
			do
			{
				for (int i = 1; i < n - 1; i++)
				{
					y[i] = (1 - w) * temp[i] + w * (a(i) * y[i - 1]
						+ a(i + 1) * temp[i + 1] + f(i) * h * h) / mtx[i][i];
				}
				r = temp - y;
				temp = y;
				k++;
			} while (abs(r) > EPS);
			if (k < pre_k && w>0) {
				pre_k = k;
				best = y;
			}
			std::cout << w << " \t" << k << '\n';
		}
		std::cout << '\n';
		return best;
	}
};
class GradientDescent : Tridiagonal
{
public:
	GradientDescent()
	{
		InitializeMatrix();
	}
	const Matrix& getMatrix() const override
	{
		return mtx;
	};

	const Vector& getAlphas() const override
	{
		return alphas;
	}
	const Vector& getBettas() const override
	{
		return betas;
	}

	const Vector& getSolutions() override
	{
		Vector temp;
		Vector r;
		double t = 0.0;
		int k = 0;
		do
		{
			for (int i = 1; i < n - 1; i++)
			{
				double ai = mtx[i][i - 1];
				double ci = mtx[i][i];
				double bi = mtx[i][i + 1];

				r[i] = -ai * y[i - 1] + ci * y[i] - bi * y[i + 1] - f(i) * h * h;
			}
			t = r * r / ((mtx * r) * r);
			y = y - t * r;
			k++;
		} while (abs(r) > EPS);
		std::cout << "Количество итераций: " << k << '\n';
		return y;
	}
};

int main()
{
	setlocale(LC_ALL, "Russian");
	Progonka method1;
	Matrix mtx = method1.getMatrix();
	std::cout << mtx << '\n';
	std::cout << "Метод прогонки: \n";
	Vector y1 = method1.getSolutions();
	std::cout << "ih \t\t yi \t\t u(ih) \t\t |yi-u(ih)|\n";
	for (int i = 0; i < n; i++)
		std::cout << i * h << "\t\t" << y1[i] <<
		"\t\t" << u(i * h) << "\t\t" << fabs(y1[i] - u(i * h)) << '\n';
	std::cout << "\n\n";

	std::cout << "Метод Якоби: \n";
	Jacobi method6;
	Vector y6 = method6.getSolutions();
	std::cout << "ih \t\t yi \t\t u(ih) \t\t |yi-u(ih)|\n";
	for (int i = 0; i < n; i++)
		std::cout << i * h << "\t\t" << y6[i] <<
		"\t\t" << u(i * h) << "\t\t" << fabs(y6[i] - u(i * h)) << '\n';
	std::cout << "\n\n";

	std::cout << "Метод верхней релаксации:\n";
	UpRelaxation method4;
	std::cout << "w \t k \n";
	Vector y4 = method4.getSolutions();
	std::cout << "\n\n";

	std::cout << "Метод наискорейшего спуска: \n";
	GradientDescent method7;
	Vector y7 = method7.getSolutions();
	std::cout << "ih \t\t yi \t\t progonka \t\t |yi-progonka|\n";
	for (int i = 0; i < n; i++)
		std::cout << i * h << "\t\t" << y7[i] << "\t" <<
		y1[i]/*u(i * h)*/ << "\t\t" <<
		fabs(y7[i] - y1[i]/*u(i * h)*/) << '\n';
	std::cout << "\n\n";

	return 0;
}
