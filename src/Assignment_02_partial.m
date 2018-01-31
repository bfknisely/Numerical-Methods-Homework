%%%%%%%%
% ME523 Spring 2018
% Dan Haworth, January 2018
%%%%%%%%

% Modified by Brian Knisely, January 29, 2018


clear; close all; format compact; more off;

% This script solves a steady, linear, one-dimensional convection-diffusion
%   equation on the domain xmin to xmax with constant properties and with
%   Dirichlet boundary conditions at both ends:
%       ?*u*??/?x = ?*?^2(?)/?x^2
% The analytic solution is compared with the numerical solution obtained using
%   a FDM with UDS or CDS for convection and with CDS for diffusion.

%%%%% define global variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% variables that are defined as global can be accessed and changed
%   in any function
% see more under the function definitions below, and in Octave Documentation

% scalars: physical quantities
global xmin;    % minimum value of independent variable x
global xmax;    % maximum value of independent variable x
global den;     % fluid mass density
global vel;     % fluid x-direction velocity component
global dif;     % fluid diffusivity
global pe;      % Peclet number based on domain length (xmax-xmin)
global phil;    % value of phi at left-hand boundary  (at xmin)
global phir;    % value of phi at right-hand boundary (at xmax)

% scalars: grid-related quantities
global n;       % number of intervals into which the 1d domain is subdivided
global nm1;     % n-1
global np1;     % n+1 (number of grid points)
global dxrat;   % grid expansion ratio - see below
global dxmin;   % minimum spacing between two adjacent grid points
global dxmax;   % maximum spacing between two adjacent grid points
global pedxmin; % Peclet number based on dxmin
global pedxmax; % Peclet number based on dxmax

% arrays (vectors) of length np1 holding values at np1 grid points
global x;       % x position
global phiuds;  % computed phi value for UDS convection
global phicds;  % computed phi value for CDS convection
global phian;   % analytic solution for phi
global erruds;  % error for UDS convection wrt/analytic solution
global errcds;  % error for CDS convection wrt/analytic solution

% arrays (vectors) of length np1 used in solving the linear system A*phi=q
global aw ap ae;% three nonzero coefficients for each row of the linear system
global q;       % right-hand side of the linear system
global phi;     % solution of the linear system

%%%%% define the physical problem and grid specifications %%%%%%%%%%%%%%%%%%%%%%

xmin = 0.;
xmax = 1.;
den  = 1.;
vel  = 10.;
dif  = 0.2;
phil = 1.;
phir = 2.;

% dxrat is the grid expansion ratio (dxrat > 0.)
%   dxrat=1.    for a uniform grid
%   0.<dxrat<1. for a finer grid  (smaller dx) with increasing x
%   dxrat>1.    for a coarser grid (larger dx) with increasing x

n     = 25;
dxrat = 1.0;

%%%%% compute derived variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grid-related derived variables are computed in function grid below

nm1 = n - 1;
np1 = n + 1;
pe  = den*vel*( xmax - xmin ) / dif;

%%%%% define and initialize arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nodes are indexed from i=1 (corresponding to x=xmin)
%   to i=np1=n+1 (corresponding to x=xmax)
% all array element values are initialized to zero

x      = zeros(np1,1);
phiuds = zeros(np1,1);
phicds = zeros(np1,1);
phian  = zeros(np1,1);
erruds = zeros(np1,1);
errcds = zeros(np1,1);

aw     = zeros(np1,1);
ap     = zeros(np1,1);
ae     = zeros(np1,1);
q      = zeros(np1,1);
phi    = zeros(np1,1);

%%%%% define functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% here the functions are defined - they are not excuted until called
%   by name (see below)

% to access a global variable within a function, it must be declared as
%   a global variable within the function
% any variable that is not declared as a global variable within a function
%   is a local variable, even if its name is the same as that of a global
%   variable
% for example, if dxmin and dxmax were not declared as global variables
%   in function grid_1d below, the values of dxmin and dxmax computed there
%   would not change the values of xmin and xmax outside of function grid_1d

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set up the 1d grid %%%
function grid_1d

  % compute grid-point locations and related quantities

  % global variables needed by this function
  global xmin xmax
  global den vel dif
  global dxrat
  global n np1

  % global variables set by this function
  global dxmin dxmax
  global pedxmin pedxmax
  global x

  % set dxmin and dxmax
  if dxrat == 1.
      dxmin = ( xmax - xmin ) / real(n);
      dxmax = dxmin;
  else
      dx1   = ( xmax - xmin )*( 1. - dxrat ) / ( 1. - dxrat**real(n) );
      dxn   = dx1*dxrat**(n-1);
      dxmin = min(dx1,dxn);
      dxmax = max(dx1,dxn);
  end

  % set pedxmin and pedxmax
  pedxmin = den*vel*dxmin / dif;
  pedxmax = den*vel*dxmax / dif;

  % set x
  i = 1;
  x(i) = xmin;
  if dxrat == 1.
      dx = ( xmax - xmin ) / real(n);
  else
      dx = ( xmax - xmin )*( 1. - dxrat ) / ( 1. - dxrat^real(n) );
      endif
      
      for i=2:np1
          x(i)  = x(i-1) + dx;
          dx    = dx*dxrat;
      end
      
      % make sure that x(np1)=xmax, in spite of potential roundoff errors
      x(np1) = xmax;
    
end


%%% end of grid_1d %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% compute the analytic solution %%%

function analytic_1d
  % compute the analytic solution of the steady, linear, 1D convection-diffusion
  %   problem at each grid point

  % global variables needed by this function
  global xmin xmax np1 phil phir pe x

  % global variables set by this function
  global phian

  % compute the analytic solution
  rx  = 1. / ( xmax - xmin );
  for i=1:np1
      phian(i) = phil + ( phir - phil )*( exp(pe*(x(i)-xmin)*rx) - 1. ) ...
          / ( exp(pe) - 1. );
  end

  % make sure that the BCs are satisfied, in spite of potential roundoff errors
  phian(1)   = phil;
  phian(np1) = phir;

end

%%% end of analytic_1d %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% solve a tridiagonal linear system using TDMA %%%

function tdma
  % solve the tridiagonal linear system A*phi=q using the Thomas algorithm
  %   (tridiagonal matrix algorithm)

  % aw, ap, ae, q, and phi are arrays (vectors) of length np1
  % aw(i), ap(i), and ae(i) are the three nonzero coefficients in row i of A,
  %   with ap(i) being on the diagonal
  % q(i) is the right-hand side value for row i
  % this function computes phi(i) for i=2,...,n (where n=np1-1)
  % in addition to aw, ap, ae, and q, the value of phi(np1) must be set
  %   before calling this function

  % global variables needed by this function
  global n np1 aw ap ae q

  % global variables set by this function
  global phi

  % define and initialize local work arrays
  bpr = zeros(np1,1);
  v   = zeros(np1,1);

  % calculate 1./u_p (bpr) and modified source term (v)
  for i=2:n
      bpr(i) = 1. / ( ap(i) - aw(i)*ae(i-1)*bpr(i-1) );
      v(i)   = q(i) - aw(i)*v(i-1)*bpr(i-1);
  end

  % solve for phi by back substitution
  for i=n:-1:2 % this decrements the value of i by 1 on each pass
      phi(i) = ( v(i) - ae(i)*phi(i+1) )*bpr(i);
  end

end
%%% end of tdma %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% compute the numerical solution using UDS convection and CDS diffusion %%%

function numerical_1d_uds
  % compute the numerical solution of the steady, linear, 1D convection-diffusion
  %   problem at each grid point using a finite-difference method
  %   with UDS convection and CDS diffusion

  % global variables needed by this function
  global phil phir np1 n den vel dif x
  global aw ap ae q phi

  % global variables set by this function
  global phiuds

  % apply the boundary conditions to the discrete solution
  phiuds(1)   = phil;
  phiuds(np1) = phir;

  % set up the discrete system of algebraic equations
  %   for a FDM with Dirichlet BCs at x=xmin and x=xmax
  % we do not need to solve for phiuds(1) or for phiuds(np1)
  % arrays aw(i), ap(i), and ae(i) hold the three nonzero values
  %   of the coefficients for row i of the coefficient matrix
  % array q(i) holds the right-hand-side for row i

  denvel = den*vel;

  %--- assemble the coefficient matrix ---

  % loop over interior grid points (i=2 to i=n)
  for i=2:n
      
      % UDS convection contributions
      awc = -max(denvel,0.) / ( x(i)   - x(i-1) );
      aec =  min(denvel,0.) / ( x(i+1) - x(i)   );
      
      % CDS diffusion contributions
      dxr = 2. / ( x(i+1) - x(i-1) );
      awd = -dif*dxr / ( x(i)   - x(i-1) );
      aed = -dif*dxr / ( x(i+1) - x(i)   );
      
      % fill row i coefficient matrix and right-hand-side values
      aw(i) = awc + awd;
      ae(i) = aec + aed;
      ap(i) = -aw(i) - ae(i);
      q(i)  = 0.;
      
  end

  % take care of the boundary conditions
  i     = 2;
  q(i)  = q(i) - aw(i)*phiuds(i-1);
  aw(i) = 0.;

  i     = n;
  q(i)  = q(i) - ae(i)*phiuds(i+1);
  ae(i) = 0.;

  phi(1)   = phiuds(1);
  phi(np1) = phiuds(np1);

  %--- solve for phi using the TDMA direct solver ---

  tdma;

  % copy the solution from phi to phiuds
  for i=2:n
      phiuds(i) = phi(i);
  end

end

%%% end of numerical_1d_uds %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% compute the numerical solution using CDS convection and CDS diffusion %%%

function numerical_1d_cds
  % compute the numerical solution of the steady, linear, 1D convection-diffusion
  %   problem at each grid point using a finite-difference method
  %   with CDS convection and CDS diffusion

  % Call global variables needed by this function
  global phil phir np1 n den vel dif x
  global aw ap ae q phi

  % global variables set by this function
  global phicds Pe

  % Apply BCs
  phicds(1) = phil;
  phicds(np1) = phir;

  % compute product of density and velocity
  denvel = den*vel;

  % loop over interior grid points (i = 2 to i = n)
  for i = 2:n
      
      % CDS convection contributions
      awc = -denvel / ( x(i+1) - x(i-1) );  % West coefficient
      aec =  denvel / ( x(i+1) - x(i-1)   );  % East coefficient
      apc = 0;  % Coefficient at point
    
      % CDS diffusion contributions (unchanged)
      %dxr = 2.*dif / ( x(i) - x(i-1) );  % 
      awd = -dif / (x(i) - x(i-1))**2;  % West coefficient
      aed = -dif / (x(i) - x(i-1))**2;  % East coefficient
      apd = 2*dif / ( x(i) - x(i-1) )**2;  % Coefficient at point
      
      % fill row i coefficient matrix and right-hand-side values
      aw(i) = awc + awd;
      ae(i) = aec + aed;    
      ap(i) = apc + apd;
      q(i)  = 0.;
  end


  % take care of the boundary conditions
  i     = 2;
  q(i)  = q(i) - aw(i)*phicds(i-1);
  aw(i) = 0.;

  i     = n;
  q(i)  = q(i) - ae(i)*phicds(i+1);
  ae(i) = 0.;

  phi(1)   = phicds(1);
  phi(np1) = phicds(np1);

  %--- solve for phi using the TDMA direct solver ---

  tdma;

  % copy the solution from phi to phicds
  for i=2:n
      phicds(i) = phi(i);
  end

end

%%% end of numerical_1d_cds %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% compute errors of the numerical solutions wrt/the analytic solution %%%

function errors
  % compute differences between two numerical solutions and the analytic solution

  # global inputs to this function
  global phicds phiuds phian

  # global outputs to this function
  global errcds erruds

  erruds = abs(phian - phiuds);  % compute error in uds solution
  errcds = abs(phian - phicds);  % compute error in cds solution

end

%%% end of errors %%%

%%% generate figures %%%

function figures
  % plot the analytic solution and two numerical solutions as functions of x
  % plot errors between two numerical solutions and the analytic solution
  %   as functions of x
  % save both figures to the current working directory as .png files

  global x phian phiuds phicds erruds errcds

  % Make figure for phi vs x for analytical, UDS, and CDS solutions
  figure(1)  % Create Figure "1"
  % Plot Phi_analytic vs x with black line
  % Plot Phi_UDS vs x with dash-dot blue line
  % Plot Phi_CDS vs x with dotted red line
  p = plot(x, phian, 'k', x, phiuds, 'b-.', x, phicds, 'r:');
  set(p, 'linewidth', 3)  % Set line thickness to 3
  xl = xlabel('x'); yl = ylabel('\phi');  % Add axis labels
  % Add legend and specify its location
  l = legend('Analytic', 'UDS1', 'CDS2', 'location', 'northwest');
  % Increase font size of plot elements
  set([gca, l, xl, yl], 'fontsize', 16);
  figure(1, 'position', [50 50 600 450])
  saveas(1,'hw2_figure1.png');  % Save figure to file "hw2_figure1.png"
  close(1)  % Close the figure

  % Make figure for |error| vs x for UDS and CDS solutions
  figure(2)  % Create Figure "2"
  p = plot(x, erruds, 'b-.', x, errcds, 'r:');
  set(p, 'linewidth', 3)  % Set line thickness to 3
  xl = xlabel('x'); yl = ylabel('|Error|');
  % Add legend and specify its location
  l = legend('UDS1', 'CDS2', 'location', 'northwest');
  % Increase font size of plot elements
  set([gca, l, xl, yl], 'fontsize', 16);
  figure(2, 'position', [650 50 600 450])
  saveas(2,'hw2_figure2.png');   % Save figure to file "hw2_figure2.png"
  close(2)  % Close the figure

end

%%% end of figures %%%

%%% write an output file %%%

function outfile
% write an ascii text output file
% echoes input parameters and key derived quantitites,
%   and generates six columns of output for each grid point:
%     x phian phiuds phicds erruds errcds

% Call global variables needed to state at top of text file
global xmin xmax den vel dif phil phir n dxrat dxmin dxmax pe pedxmin pedxmax


% Call global variables needed for table
global x phian phiuds phicds erruds errcds

  % Name output file and open it
  outfile = 'hw2_results.txt'
  fid = fopen(outfile, 'w');  % Open a file to store results

  % Print each variable in nice formatting
  fprintf(fid, 'Input quantities:\n')
  fprintf(fid, '%8s = %12.6f\n', 'xmin', xmin)
  fprintf(fid, '%8s = %12.6f\n', 'xmax', xmax)
  fprintf(fid, '%8s = %12.6f\n', 'den', den)
  fprintf(fid, '%8s = %12.6f\n', 'vel', vel)
  fprintf(fid, '%8s = %12.6f\n', 'dif', dif)
  fprintf(fid, '%8s = %12.6f\n', 'phil', phil)
  fprintf(fid, '%8s = %12.6f\n', 'phir', phir)
  fprintf(fid, '%8s = %12.6f\n', 'n', n)
  fprintf(fid, '%8s = %12.6f\n', 'dxrat', dxrat)
  
  fprintf(fid, '\nDerived quantities:\n')
  fprintf(fid, '%8s = %12.6f\n', 'dxmin', dxmin)
  fprintf(fid, '%8s = %12.6f\n', 'dxmax', dxmax)
  fprintf(fid, '%8s = %12.6f\n', 'pe', pe)
  fprintf(fid, '%8s = %12.6f\n', 'pedxmin', pedxmin)
  fprintf(fid, '%8s = %12.6f\n', 'pedxmax', pedxmax)

  % Print header
  fprintf(fid, '\nTabulated results:\n')
  fprintf(fid, '\n%3s,%12s,%12s,%12s,%12s,%12s\n',...
      'x', 'phian', 'phiuds', 'phicds', 'erruds', 'errcds')

  % Print results one line at a time
  for i = 1:length(phian)
      fprintf(fid, '%3.2f,%12.6e,%12.6e,%12.6e,%12.6e,%12.6e\n',...
      x(i), phian(i), phiuds(i), phicds(i), erruds(i), errcds(i))
  end

  fclose(fid);  % Close the file

end
%%% end of outfile %%%

%%%%% execute the program %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate the grid
grid_1d;

% compute the analytic solution
analytic_1d;

% compute the UDS convection solution
numerical_1d_uds;

% compute the CDS convection solution
numerical_1d_cds;

% compute errors with respect to the analytic solution
errors;

% plot the results
figures;

% write an output file with the results
outfile;

%%%%% all done %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%