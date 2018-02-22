%%%%%%%%
% ME523 Spring 2018
% Dan Haworth, January 2018
%%%%%%%%

% Modified by Brian Knisely, January 29, 2018


clear; close all; format compact; more off;

nn = [10, 20, 40, 160, 320, 640];
dxrat_arr = [1., 1., 1., 1., 1., 1.];
comb = zeros((length(nn))*length(dxrat_arr), 3);

ind = 0;

%errU = zeros((length(nn)+1)*length(dxrat), max(nn));
%errC = zeros((length(nn)+1)*length(dxrat), max(nn));

for nnn = 1:length(nn);
    ind = ind + 1;
    st = sprintf('n%i',nn(ind));


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
dif  = 0.5;
phil = 1.;
phir = 2.;

% dxrat is the grid expansion ratio (dxrat > 0.)
%   dxrat=1.    for a uniform grid
%   0.<dxrat<1. for a finer grid  (smaller dx) with increasing x
%   dxrat>1.    for a coarser grid (larger dx) with increasing x

n = nn(nnn);
dxrat = dxrat_arr(nnn);

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
      awd = -dif / ( (x(i) - x(i-1))*(x(i+1) - x(i-1))/2 );  % West coefficient
      aed = -dif / ( (x(i+1) - x(i-1))*(x(i+1) - x(i))/2 );  % East coefficient
      apd = 2*dif / ( (x(i) - x(i-1))*(x(i+1) - x(i)) );  % Coefficient at point
      
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

  % global inputs to this function
  global phicds phiuds phian

  % global outputs to this function
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
  figure(2, 'position', [650 50 600 450]);
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
%figures;

% write an output file with the results
%outfile;

[a, ind90] = min(abs((x-0.9)));
phian90 = phian(ind90);


comb(ind, 1) = n;
comb(ind, 2) = dxrat;
comb(ind, 3) = ind;

if nnn < 9
  phiU.(st) = phiuds(ind90);
  phiC.(st) = phicds(ind90);
  errU.(st) = abs(erruds(ind90));
  errC.(st) = abs(errcds(ind90));
  dX.(st) = x(2)- x(1);
else
  errU.(st) = erruds;
  errC.(st) = errcds;
  X.(st) = x;
end
end  % loop for values of N


%%%%% all done %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pu = [];
pc = [];

for ii = 1:6;
  st = sprintf('n%i',nn(ii));
  pu(end+1) = errU.(st);
  pc(end+1) = errC.(st);
%  if ii <= 3
%    figure(1); hold on;
%  else
%    figure(2); hold on;
%  end
%  loglog(dX.(st), errU.(st), 'ko', 'markersize', 10, 'linewidth', 5);
%  loglog(dX.(st), errC.(st), 'rx', 'markersize', 10,'linewidth', 5);
%  l = legend('UDS1', 'CDS2', ...
%  'location', 'northwest');
%  xl = xlabel('\Deltax'); yl = ylabel('|Error|');
%  set([gca, l, xl, yl], 'fontsize', 18);
%  figure(gcf, 'position', [650 50 600 450]);
end
%figure(gcf, 'position', [50 50 600 450]);

% Split values of pu and pc for low gridpoints case and high gridpoints case
pu_coarse = pu(1:3);  % For low numbers of gridpoints (10, 20, 40)
pc_coarse = pc(1:3);  % For low numbers of gridpoints (10, 20, 40)
pu_fine = pu(4:6);  % For high numbers of gridpoints (160, 320, 640)
pc_fine = pc(4:6);  % For high numbers of gridpoints (160, 320, 640)

% Calculate convergence rate (P) for coarse grids (10, 20, 40)
% UDS case
P_uds_coarse = 1/log(2) * log((pu_coarse(end-1)-pu_coarse(end-2)) / (pu_coarse(end)-pu_coarse(end-1)));
% CDS case
P_cds_coarse = 1/log(2) * log((pc_coarse(end-1)-pc_coarse(end-2)) / (pc_coarse(end)-pc_coarse(end-1)));

% Calculate convergence rate (P) for fine grids (160, 320, 640)
% UDS case
P_uds_fine = 1/log(2) * log((pu_fine(end-1)-pu_fine(end-2)) / (pu_fine(end)-pu_fine(end-1)));
% CDS case
P_cds_fine = 1/log(2) * log((pc_fine(end-1)-pc_fine(end-2)) / (pc_fine(end)-pc_fine(end-1)));

phiEst_uds_coarse = phiU.n40 + P_uds_coarse * dX.n40;
phiEst_cds_coarse = phiC.n40 + P_cds_coarse * dX.n40;
phiEst_uds_fine = phiU.n640 + P_uds_fine * dX.n640;
phiEst_cds_fine = phiC.n640 + P_cds_fine * dX.n640;

% First use CDS for all parts
% Part a
fprintf('\nPart a:\n')
phiEst_cds_coarse
P_cds_coarse
errC.n40
phian90

% Part b
fprintf('\nPart b:\n')
phiEst_cds_fine
P_cds_fine
errC.n640
phian90

% Part c
dxCoarse = [dX.n10, dX.n20, dX.n40];

errCoarseEst(1) = abs(phiC.n10-phiEst_cds_coarse);
errCoarseEst(2) = abs(phiC.n20-phiEst_cds_coarse);
errCoarseEst(3) = abs(phiC.n40-phiEst_cds_coarse);

errCoarseAn(1) = abs(phiC.n10-phian90);
errCoarseAn(2) = abs(phiC.n20-phian90);
errCoarseAn(3) = abs(phiC.n40-phian90);

dxFine = [dX.n160, dX.n320, dX.n640];

errFineEst(1) = abs(phiC.n160-phiEst_cds_fine);
errFineEst(2) = abs(phiC.n320-phiEst_cds_fine);
errFineEst(3) = abs(phiC.n640-phiEst_cds_fine);

errFineAn(1) = abs(phiC.n160-phian90);
errFineAn(2) = abs(phiC.n320-phian90);
errFineAn(3) = abs(phiC.n640-phian90);

p = loglog(dxCoarse, errCoarseEst, 'ko-', dxCoarse, errCoarseAn, 'rs-.',...
     dxFine, errFineEst, 'bd--', dxFine, errFineAn, 'm:p');
xlabel('\Deltax', 'fontsize', 18);
ylabel('|Error|', 'fontsize', 18);
leg = legend('wrt CDS, coarse, RE', 'wrt Analytic', 'wrt CDS, fine, RE', ...
       'wrt Analytic', 'location', 'southeast');
set(gca, 'fontsize', 18);
set(leg, 'fontsize', 16);
set(p, 'linewidth', 3);

% Part d is to do everything again but with UDS instead of CDS
% Part d - like part a
fprintf('\nPart da:\n')
phiEst_uds_coarse
P_uds_coarse
errU.n40
phian90

% Part d - like part b
fprintf('\nPart db:\n')
phiEst_uds_fine
P_uds_fine
errU.n640
phian90