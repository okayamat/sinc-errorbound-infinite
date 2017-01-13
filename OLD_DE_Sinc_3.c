#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f(double t)
{
  return pow(t, M_PI_4) * exp(- t);
}

double ftrans(double y)
{
  double x = M_PI_2*sinh(y);
  double t;
  if (x >= 0) {
    t = x + log1p(exp(- x));
  } else {
    t = log1p(exp(x));
  }

  return pow(t, M_PI_4) / (1 + exp(x));
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

double sinc_approx(int N, double d, double mu, double t)
{
  int k;
  double h = log(4 * d * N / mu) / N;
  double x = asinh(M_2_PI * log(expm1(t)));
  double val = 0;

  for (k = N; k > 0; k--) {
    val += ftrans( k*h) * S( k,h,x);
    val += ftrans(-k*h) * S(-k,h,x);
  } val += ftrans( 0*h) * S( 0,h,x);

  return val;
}
/*** This error bound is NOT correct! ***/
double err_bound(int N, double K, double d, double mu)
{
  double val = 4.0 / (M_PI * (1 - exp(- M_PI * mu * M_E)));
  val /= (pow(cos(M_PI_2 * sin(d)),2*mu) * pow(cos(d),1+mu));
  val += mu * exp(M_PI_2*mu + 1) * pow(2, 1-mu);
  val *= K / (pow(M_PI,1-mu)*d*mu);

  return val * exp(- M_PI*d*N/log(2*d*N/mu));
}

int main()
{
  int i, N;
  double K = 1.0;
  double d = 1.4;
  double mu = M_PI_4;
  double s, t;
  double err, maxerr;

  for (N = 2; N <= 50; N += 5) {
    s = 1;
    t = 1;
    maxerr = fabs(f(t) - sinc_approx(N, d, mu, t));

    for (i = 1; i <= 100; i++) {
      s *= M_SQRT1_2;
      t *= M_SQRT2;

      err = fabs(f( s) - sinc_approx(N, d, mu, s));
      if (maxerr < err)
        maxerr = err;

      err = fabs(f( t) - sinc_approx(N, d, mu, t));
      if (maxerr < err)
        maxerr = err;
    }

    printf("%d\t%e\t%e\n", N, maxerr, err_bound(N, K, d, mu));
  }

  return EXIT_SUCCESS;
}
