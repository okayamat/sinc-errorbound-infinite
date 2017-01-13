#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f(double t)
{
  double val = tanh(log(t));

  return sqrt(1 + val * val) * sqrt(t)/(1 + t*t);
}

double ftrans(double y)
{
  double x = M_PI_2*sinh(y);
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
  double h = log(4 * d * n / mu) / n;
  if (alpha <= beta) {
    M = n;
    N = n - (int)floor(log(beta/alpha)/h);
  } else {
    M = n - (int)floor(log(alpha/beta)/h);
    N = n;
  }
  double x = asinh(M_2_PI * log(t));
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
  double mu = fmin(alpha, beta);
  double nu = fmax(alpha, beta);
  double val = 4.0 / (M_PI * (1 - exp(- M_PI_2 * mu * M_E)));
  val /= (pow(cos(M_PI_2 * sin(d)),(alpha+beta)*0.5) * cos(d));
  val += mu * exp(M_PI_4*nu);
  val *= 2.0 * K / (M_PI*d*mu);

  return val * exp(- M_PI*d*N/log(4*d*N/mu));
}

int main()
{
  int i, N;
  double alpha = 0.5;
  double beta  = 1.5;
  double K = 1.5;
  double d = M_PI / 6.0;
  double s, t;
  double err, maxerr;

  for (N = 2; N <= 120; N += 5) {
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
