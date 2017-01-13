#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f(double t)
{
  double val = tanh(log(t));

  return sqrt(1 + val * val) * sqrt(t)/(1 + t*t);
}

double ftrans(double x)
{
  double val = tanh(x);

  return sqrt(1 + val * val) * 0.5 / (exp(0.5*x)*cosh(x));
}

double sinc(double x)
{
  double val = 1.0;

  if (x != 0) {
    val = sin(M_PI * x) / (M_PI * x);
  }

  return val;
}

double S(int k, double h, double x)
{
  return sinc((x/h) - k);
}

double sinc_approx(int n, double d, double alpha, double beta, double t)
{
  int k, M, N;
  double mu= fmin(alpha, beta);
  double h = sqrt(M_PI * d / (mu*n));
  if (alpha <= beta) {
    M = n;
    N = (int)ceil(alpha * n / beta);
  } else {
    M = (int)ceil(beta * n / alpha);
    N = n;
  }
  double x = log(t);
  double val1 = 0;
  double val2 = 0;

  for (k = -M; k < 0; k++) {
    val1 += ftrans(k*h) * S(k,h,x);
  }
  for (k = N; k >= 0; k--) {
    val2 += ftrans(k*h) * S(k,h,x);
  }

  return val1 + val2;
}

double err_bound(int N, double K, double d, double alpha, double beta)
{
  double mu  = fmin(alpha, beta);
  double pdm = sqrt(M_PI * d * mu);
  double sqn = sqrt(N);
  double val = 1.0;
  val += 2.0 / (pdm * (1 - exp(-2 * pdm)) * pow(cos(d),(alpha+beta)*0.5));
  val *= 2.0 * K / pdm;

  return val * sqn * exp(- pdm * sqn);
}

int main()
{
  int i, N;
  double alpha = 0.5;
  double beta  = 1.5;
  double K = 1.5;
  double d = M_PI_4;
  double s, t;
  double err, maxerr;

  for (N = 2; N <= 140; N += 5) {
    s = 1;
    t = 1;
    maxerr = fabs(f(t) - sinc_approx(N, d, alpha, beta, t));

    for (i = 1; i <= 100; i++) {
      s *= M_SQRT1_2;
      t *= M_SQRT2;

      err = fabs(f(s) - sinc_approx(N, d, alpha, beta, s));
      if (maxerr < err)
        maxerr = err;

      err = fabs(f(t) - sinc_approx(N, d, alpha, beta, t));
      if (maxerr < err)
        maxerr = err;
    }

    printf("%d\t%e\t%e\n", N, maxerr, err_bound(N, K, d, alpha, beta));
  }

  return EXIT_SUCCESS;
}
