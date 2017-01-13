#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f(double t)
{
  double val = tanh(asinh(t));

  return sqrt(1 + val * val) / (1 + t*t);
}

double ftrans(double y)
{
  double x = M_PI_2*sinh(y);
  double t = tanh(x);
  double c = cosh(x);

  return (sqrt(1 + t * t) / c) / c;
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
  double x = asinh(M_2_PI * asinh(t));
  double val = 0;

  for (k = N; k > 0; k--) {
    val += ftrans( k*h) * S( k,h,x);
    val += ftrans(-k*h) * S(-k,h,x);
  } val += ftrans( 0*h) * S( 0,h,x);

  return val;
}

double err_bound(int N, double K, double d, double mu)
{
  double val = 4.0 / (M_PI * (1 - exp(-0.5 * M_PI * mu * M_E)));
  val /= (pow(cos(M_PI_2 * sin(d)),mu) * cos(d));
  val += mu * exp(M_PI_4*mu);
  val *= pow(2,mu+1) * K / (M_PI*d*mu);

  return val * exp(- M_PI*d*N/log(4*d*N/mu));
}

int main()
{
  int i, N;
  double K = 1.5;
  double d = M_PI / 6.0;
  double mu= 2.0;
  double s, t;
  double err, maxerr;

  for (N = 2; N <= 95; N += 5) {
    t = 0;
    maxerr = fabs(f(t) - sinc_approx(N, d, mu, t));

    s = 1;
    t = 1;
    err = fabs(f( t) - sinc_approx(N, d, mu, t));
    if (maxerr < err)
      maxerr = err;

    err = fabs(f(-t) - sinc_approx(N, d, mu,-t));
    if (maxerr < err)
      maxerr = err;

    for (i = 1; i <= 100; i++) {
      s *= M_SQRT1_2;
      t *= M_SQRT2;

      err = fabs(f( s) - sinc_approx(N, d, mu, s));
      if (maxerr < err)
        maxerr = err;

      err = fabs(f(-s) - sinc_approx(N, d, mu,-s));
      if (maxerr < err)
        maxerr = err;

      err = fabs(f( t) - sinc_approx(N, d, mu, t));
      if (maxerr < err)
        maxerr = err;

      err = fabs(f(-t) - sinc_approx(N, d, mu,-t));
      if (maxerr < err)
        maxerr = err;
    }

    printf("%d\t%e\t%e\n", N, maxerr, err_bound(N, K, d, mu));
  }

  return EXIT_SUCCESS;
}
