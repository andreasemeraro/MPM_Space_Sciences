function [ delta, e, t0, y0, tstop ] = kepler_perturbed_parameters ( ...
  delta_user, e_user, t0_user, y0_user, tstop_user )

%*****************************************************************************80
%
%% kepler_perturbed_parameters(): parameters for the perturbed Kepler ODE.
%
%  Discussion:
%
%    If input values are specified, this resets the default parameters.
%    Otherwise, the output will be the current defaults.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    23 April 2021
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    real DELTA_USER: ?;
%
%    real E_USER: the orbital eccentricity;
%
%    real T0_USER: the initial time;
%
%    real Y0_USER(4): the initial values.
%
%    real TSTOP_USER: the final time.
%
%  Output:
%
%    real DELTA: ?;
%
%    real E: the orbital eccentricity;
%
%    real T0: the initial time;
%
%    real Y0(4): the initial values.
%
%    real TSTOP: the final time.
%
  persistent delta_default;
  persistent e_default;
  persistent t0_default;
  persistent y0_default;
  persistent tstop_default;
%
%  Initialize defaults.
%
  if ( isempty ( delta_default ) )
    delta_default = 0.015;
  end

  if ( isempty ( e_default ) )
    e_default = 0.6;
  end

  if ( isempty ( t0_default ) )
    t0_default = 0.0;
  end

  if ( isempty ( y0_default ) )
    p0 = 1.0 - e_default;
    p1 = 0.0;
    q0 = 0.0;
    q1 = sqrt ( ( 1.0 + e_default ) / ( 1.0 - e_default ) );
    y0_default = [ p0, p1, q0, q1 ];
  end

  if ( isempty ( tstop_default ) )
    tstop_default = 120.0;
  end
%
%  Update defaults if input was supplied.
%
  if ( 1 <= nargin )
    delta_default = delta_user;
  end

  if ( 2 <= nargin )
    e_default = e_user;
  end

  if ( 3 <= nargin )
    t0_default = t0_user;
  end

  if ( 4 <= nargin )
    y0_default = y0_user;
  end

  if ( 5 <= nargin )
    tstop_default = tstop_user;
  end
%
%  Return values.
%
  delta = delta_default;
  e = e_default;
  t0 = t0_default;
  y0 = y0_default;
  tstop = tstop_default;

  return
end

