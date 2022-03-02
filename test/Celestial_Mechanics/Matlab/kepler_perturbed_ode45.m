function kepler_perturbed_ode45 ( )

%*****************************************************************************80
%
%% kepler_perturbed_ode45() uses ode45() on the perturbed Kepler ODE.
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
  fprintf ( 1, '\n' );
  fprintf ( 1, 'kepler_perturbed_ode45():\n' );
  fprintf ( 1, '  Use ode45() to solve the perturbed Kepler ODE.\n' );

  [ delta, e, t0, y0, tstop ] = kepler_perturbed_parameters ( )

  f = @ kepler_perturbed_deriv
  tspan = [ t0, tstop ];

  [ t, y ] = ode45 ( f, tspan, y0 );

  n = size ( t, 1 );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Number of variable size steps = %d\n', n );
%
%  Time plot.
%
  figure ( );
  clf ( );
  plot ( t, y, 'linewidth', 2 );
  grid ( 'on' );
  xlabel ( '<---  t  --->' );
  ylabel ( '<---  y(1:4)  --->' );
  title ( 'Perturbed Kepler ode45: time plot' )
  filename = 'kepler_perturbed_ode45_plot.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Phase plot.
%
  figure ( );
  clf ( );
  plot ( y(:,1), y(:,2) );
  grid ( 'on' );
  xlabel ( '<---  y1  --->' );
  ylabel ( '<---  y2  --->' );
  title ( 'Perturbed Kepler ode45: phase' )
  filename = 'kepler_perturbed_ode45_phase.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );

  return
end
