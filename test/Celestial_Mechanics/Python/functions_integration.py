import math
import numpy as np

mu = 398600
def rv_from_r0v0(R0, V0, t):
  #  ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃
  # This function computes the state vector (R,V) from the
  # initial state vector (R0,V0) and the elapsed time.
  #
  # mu - gravitational parameter (kmˆ3/sˆ2)
  # R0 - initial position vector (km)
  # V0 - initial velocity vector (km/s)
  # t - elapsed time (s)
  # R - final position vector (km)
  # V - final velocity vector (km/s)
  #
  # User M-functions required: kepler_U, f_and_g, fDot_and_gDot
  # ------------------------------------------------------------

  global mu
  #...Magnitudes of R0 and V0:
  r0 = np.linalg.norm(R0);
  v0 = np.linalg.norm(V0);
  #...Initial radial velocity:
  
  vr0 = np.dot(R0[0], V0[0])/r0;
  #l'espressione corretta, trattandosi di un prodotto scalare, dovrebbe essere la seguente
  #vr0 = np.dot(R0, V0)/r0;
  
  #...Reciprocal of the semimajor axis (from the energy equation):
  alpha = 2/r0 - v0**2/mu;
  #...Compute the universal anomaly:
  x = kepler_U(t, r0, vr0, alpha);
  #...Compute the f and g functions:
  f, g = f_and_g(x, t, r0, alpha);
  #print(f,g,R0,V0)
  #...Compute the final position vector:
  R = f*R0 + g*V0;
  #...Compute the magnitude of R:
  r = np.linalg.norm(R);
  #...Compute the derivatives of f and g:
  fdot, gdot = fDot_and_gDot(x, r, r0, alpha);
  #...Compute the final velocity:
  V = fdot*R0 + gdot*V0;
  return R,V




def kepler_U(dt, ro, vro, a):
  #  ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃
  #
  # This function uses Newton’s method to solve the universal
  # Kepler equation for the universal anomaly.
  #
  # mu - gravitational parameter (kmˆ3/sˆ2)
  # x - the universal anomaly (kmˆ0.5)
  # dt - time since x = 0 (s)
  # ro - radial position (km) when x = 0
  # vro - radial velocity (km/s) when x = 0
  # a - reciprocal of the semimajor axis (1/km)
  # z - auxiliary variable (z = a*xˆ2)
  # C - value of Stumpff function C(z)
  # S - value of Stumpff function S(z)
  # n - number of iterations for convergence
  # nMax - maximum allowable number of iterations
  #
  # User M-functions required: stumpC, stumpS
  # ------------------------------------------------------------
  global mu
  #...Set an error tolerance and a limit on the number of
  # iterations:
  error = 1.e-8;
  nMax = 1000;
  #..Starting value for x:
  x = np.sqrt(mu)*np.abs(a)*dt;
  #...Iterate on Equation 3.62 until convergence occurs within
  #...the error tolerance:
  n = 0;
  ratio = 1;
  while abs(ratio) > error and n <= nMax:
    n = n + 1;
    C = stumpC(a*x**2);
    S = stumpS(a*x**2);
    F = ro*vro/np.sqrt(mu)*x**2*C + (1 - a*ro)*x**3*S + ro*x-np.sqrt(mu)*dt;
    dFdx = ro*vro/np.sqrt(mu)*x*(1 - a*x**2*S)+(1 - a*ro)*x**2*C+ro;
    ratio = F/dFdx;
    x = x - ratio;

  #...Deliver a value for x, but report that nMax was reached:
  if n > nMax:
    print('\n **No. iterations of Kepler''s equation')
    print(' = %g', n)
    print('\n F/dFdx = %g\n', F/dFdx)

  return x

def f_and_g(x, t, ro, a):
  #  ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃
  #
  # This function calculates the Lagrange f and g coefficients.
  #
  # mu - the gravitational parameter (kmˆ3/sˆ2)
  # a - reciprocal of the semimajor axis (1/km)
  # ro - the radial position at time t (km)
  # t - the time elapsed since t (s)
  # x - the universal anomaly after time t (kmˆ0.5)
  # f - the Lagrange f coefficient (dimensionless)
  # g - the Lagrange g coefficient (s)
  #
  # User M-functions required: stumpC, stumpS
  # ------------------------------------------------------------
  global mu
  z = a*x**2;
  #...Equation 3.66a:
  f = 1 - x**2/ro*stumpC(z);
  #...Equation 3.66b:
  g = t - 1/np.sqrt(mu)*x**3*stumpS(z);
  return f,g


def fDot_and_gDot(x, r, ro, a):
  #  ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃
  #
  # This function calculates the time derivatives of the
  # Lagrange f and g coefficients.
  #
  # mu - the gravitational parameter (kmˆ3/sˆ2)
  # a - reciprocal of the semimajor axis (1/km)
  # ro - the radial position at time t (km)
  # t - the time elapsed since initial state vector (s)
  # r - the radial position after time t (km)
  # x - the universal anomaly after time t (kmˆ0.5)
  # fDot - time derivative of the Lagrange f coefficient (1/s)
  # gDot - time derivative of the Lagrange g coefficient
  # (dimensionless)
  #
  # User M-functions required: stumpC, stumpS
  # ------------------------------------------------------------
  global mu
  z = a*x**2;
  #...Equation 3.66c:
  fdot = np.sqrt(mu)/r/ro*(z*stumpS(z) - 1)*x;
  #...Equation 3.66d:
  gdot = 1 - x**2/r*stumpC(z);
  return fdot,gdot


def stumpS(z):
  #  ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃
  #
  # This function evaluates the Stumpff function S(z) according
  # to Equation 3.49.
  #
  # z - input argument
  # s - value of S(z)
  #
  #User M-functions required: none
  # ------------------------------------------------------------
  if z > 0:
    s = (np.sqrt(z) - np.sin(np.sqrt(z)))/(np.sqrt(z))**3;
  elif z < 0:
    s = (np.sinh(np.sqrt(-z)) - np.sqrt(-z))/(np.sqrt(-z))**3;
  else:
    s = 1/6;
  return s


def stumpC(z):
  #  ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃
  #
  # This function evaluates the Stumpff function C(z) according
  # to Equation 3.50.
  #
  # z - input argument
  # c - value of C(z)
  #
  # User M-functions required: none
  # ------------------------------------------------------------
  if z > 0:
    c = (1 - np.cos(np.sqrt(z)))/z;
  elif z < 0:
    c = (np.cosh(np.sqrt(-z)) - 1)/(-z);
  else:
    c = 1/2;
  return c
