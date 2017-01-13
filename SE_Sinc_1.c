#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f(double t)
{
  double val = tanh(asinh(t));

  return sqrt(1 + val * val) / (1 + t*t);
}

double ftrans(double x)
{
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
  double h = sqrt(M_PI * d / (mu*N));
  double x = asinh(t);
  double val = 0;

  for (k = N; k > 0; k--) {
    val += ftrans( k*h) * S( k,h,x);
    val += ftrans(-k*h) * S(-k,h,x);
  } val += ftrans( 0*h) * S( 0,h,x);

  return val;
}

double err_bound(int N, double K, double d, double mu)
{
  double pdm = sqrt(M_PI * d * mu);
  double sqn = sqrt(N);
  double val = 1.0;
  val += 2.0 / (pdm * (1 - exp(-2 * pdm)) * pow(cos(d),mu));
  val *= pow(2,mu+1) * K / pdm;

  return val * sqn * exp(- pdm * sqn);
}

int main()
{
  int i, N;
  double K = 1.5;
  double d = M_PI_4;
  double mu= 2.0;
  double s, t;
  double err, maxerr;

  for (N = 2; N <= 140; N += 5) {
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
