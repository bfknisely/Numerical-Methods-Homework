%%%%%%%%
% ME523 Spring 2018
% Dan Haworth, February 2018
%%%%%%%%

% This program solves a steady, linear, two-dimensional convection-diffusion
%   equation in the rectangular domain 0<=x<=xmax, 0<=y<=ymax 
%   with constant fluid properties rho and dif:
%       rho*u*dphi/dx + rho*v*dphi/dy = dif*( d^2phi/dx^2 + d^2phi/dy^2 )
% The numerical solution is obtained using a FVM,
%   with either UDS or CDS for convection and with CDS for diffusion.
% The mesh spacings in x (dx) and in y (dy) are each uniform,
%   but dx is not necessarily equal to dy.
% The velocity field is specified, and corresponds to a stagnation-point flow
%   about the origin (requires xmin = 0, ymin = 0):
%       u = x, v = -y
%   This satisfies the incompressible continuity equation du/dx + dv/dy = 0.
% The boundary conditions on phi are as follows:
%   x=xmin=0: phi varies linearly from 0. at y=ymax to 1. at y=ymin=0
%   x=xmax:   dphi/dx=0.
%   y=ymin=0: dphi/dy=0.
%   y=ymax:   phi=0.

% This is adapted from a problem in Chapter 4 of Ferziger & Peric,
%   3rd Edition, 2002.
% The Dirichlet BC implementations are formally only 1st-order accurate here,
%   while the Neumann BC implementations are formally 2nd-order accurate.
% The code executes very slowly. This could be improved significantly
%   by using array operators and something other than a direct linear
%   system solver. That has not been done here, so that the array indexing
%     is evident and so that an "exact" solution to the discrete system
%     is available, which can be used later to investigate accuracy and CPU
%     time for various iterative linear system solvers.

% Mass flow rates, and convective and diffuxive flow rates of phi, across
%   the four boundaries of the computational domain are computed. These are
%   used to check for global conservation.

%%%%% define global variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% scalars: physical quantities
global xmin;    % minimum value of independent variable x (=0)
global xmax;    % maximum value of independent variable x (=1)
global ymin;    % minimum value of independent variable y (=0)
global ymax;    % maximum value of independent variable y (=1)
global den;     % fluid mass density
global dif;     % fluid diffusivity

% scalars: grid-related quantities - see grid_2d below
global n;       % number of subintervals in the x direction
global m;       % number of subintervals in the y direction

global nm;      % total number of control volumes or cells (nm=n*m)
global nx;      % number of grid lines in x direction, plus 1 (nx=n+2)
global ny;      % number of grid lines in y direction, plus 1 (ny=m+2)

global dx;      % grid spacing in x direction (constant)
global dy;      % grid spacing in y direction (constant)

% scalars: numerical parameters
global cflag;   % convection discretization flag
                %   UDS: cflag=1
                %   CDS: cflag=2
                
% scalars: maximum grid Peclet numbers in x and in y directions
global pedx_max % max grid Peclet number magnitude based on x velocity and dx
global pedy_max % max grid Peclet number magnitude based on y velocity and dy                
                
% scalars: net flowrates across the boundaries of the computational domain -
%   see bndryflowrates below
global mdot_xmin % net mass flow rate across x = xmin boundary
global mdot_xmax % net mass flow rate across x = xmax boundary
global mdot_ymin % net mass flow rate across y = ymin boundary
global mdot_ymax % net mass flow rate across y = ymax boundary

global cdot_xmin % net convective flow rate of phi across x = xmin boundary
global cdot_xmax % net convective flow rate of phi across x = xmax boundary
global cdot_ymin % net convective flow rate of phi across y = ymin boundary
global cdot_ymax % net convective flow rate of phi across y = ymax boundary

global ddot_xmin % net diffusive flow rate of phi across x = xmin boundary
global ddot_xmax % net diffusive flow rate of phi across x = xmax boundary
global ddot_ymin % net diffusive flow rate of phi across y = ymin boundary
global ddot_ymax % net diffusive flow rate of phi across y = ymax boundary

global mdot_tot   % net mass flow rate across all four boundaries
global cdot_tot   % net convective flow rate of phi across all four boundaries
global ddot_tot   % net diffusive flow rate of phi across all four boundaries
global phidot_tot % net flowrate of phi acroess all four boundaries

% scalars: CPU times
global lufactor_cpu % CPU time to factor the matrix (forward elimination)
global lusolve_cpu  % CPU time to solve the factored system (back substitution)

% arrays (vectors) holding spatial coordinates - see grid_2d below
global xc;      % control volume (cell-center) locations in x direction (nx x 1)
global yc;      % control volume (cell-center) locations in y direction (ny x 1)
global xf;      % cell-face-center locations in x direction (nx x 1)
global yf;      % cell-face-center locations in y direction (ny x 1)

% arrays (vectors) holding cell-face geometric interpolation factors -
%   see ginterp below
global gfacx;   % face interpolation factors in x direction (nx x 1)
global gfacy;   % face interpolation factors in y direction (ny x 1)

% arrays holding cell-face mass flowrates - see facemdots below
global mdotuf;  % face x-direction mass flux (nx x ny)
global mdotvf;  % face y-direction mass flux (nx x ny)

% array holding computed nodal values of the dependent variable
global phi;     % control volume (cell-center) values of dependent variable
                %   (nx x ny)

% arrays holding nonzero elements of coefficient matrix (all nx x ny)
global awe;     % coefficient for control volume west  of control volume P
global aso;     % coefficient for control volume south of control volume P
global ap;      % coefficient on diagonal for control volume P
global ano;     % coefficient for control volume north of control volume P
global aea;     % coefficient for control volume east  of control volume P

% array holding right-hand-side values
global q;       % right-hand-side value for control volume P (nx x ny)

% arrays for coefficient matrix, right-hand-side vector, and dependent variable
%   vector using 1D indexing (indexing cells from 1 to nm)
global a1d;     % coefficient matrix (nm x nm)
global q1d;     % right-hand-side vector (nm x 1)
global phi1d;   % dependent variable (nm x 1)

% work arrays (vectors) used in the course of matrix assembly and solution
global ipvt;    % pivot element pointer in LU decomposition (nm x 1)
global dw;      % work array used in matrix assembly (ny x 1)

%%%%% define the physical problem, grid size, and numerical parameters %%%%%%%%%

xmin  = 0.;
xmax  = 1.;
ymin  = 0.;
ymax  = 1.;

den   = 100.;
dif   = 0.1;

n     = 25;
m     = 25;

cflag = 2;

%%%%% compute derived scalar variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nm = n*m;
nx = n + 2;
ny = m + 2;
dx = ( xmax - xmin ) / n;
dy = ( ymax - ymin ) / m;

%%%%% define and initialize arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% all array element values are initialized to zero

xc     = zeros(nx,1);
yc     = zeros(ny,1);
xf     = zeros(nx,1);
yf     = zeros(ny,1);

gfacx  = zeros(nx,1);
gfacy  = zeros(ny,1);

mdotuf = zeros(nx,ny);
mdotvf = zeros(nx,ny);

phi    = zeros(nx,ny);

awe    = zeros(nx,ny);
aso    = zeros(nx,ny);
ap     = zeros(nx,ny);
ano    = zeros(nx,ny);
aea    = zeros(nx,ny);
q      = zeros(nx,ny);

a1d    = zeros(nm,nm);
q1d    = zeros(nm,1);
phi1d  = zeros(nm,1);

ipvt   = zeros(nm,1);
dw     = zeros(ny,1);

%%%%% define functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% set up the 2d grid %%%

function grid_2d
% compute cell-center and cell-face-center locations

% global variables needed by this function
  global xmin xmax ymin ymax dx dy
  global nx ny
  
% global variables set by this function
  global xf yf xc yc

% --- indexing conventions ---
%   in the x-direction:
%     n is the number of cells or control volumes in the x direction
%     nx  = n + 2
%     xf(i) is the x-location of grid-line (or face) i, 
%       with an extra entry at xmax
%         xf(1)   = xmin
%         xf(2)   = xmin + dx
%         . . .
%         xf(n)   = xmax - dx
%         xf(n+1) = xmax
%         xf(n+2) = xmax
%     xc(i) is the x-location of control volume (or cell center) i,
%       with an extra entry at xmin and an extra entry at xmax
%         xc(1)   = xmin
%         xc(2)   = xmin + dx/2
%         xc(3)   = xmin + 3*dx/2
%         . . .
%         xc(n)   = xmax - 3*dx/2
%         xc(n+1) = xmax - dx/2
%         xc(n+2) = xmax
%   in the y-direction:
%     m is the number of cells or control volumes in the y direction
%     ny  = m + 2
%     yf(j) is the y-location of grid-line (or face) j,
%       with an extra entry at ymax
%         yf(1)   = ymin
%         yf(2)   = ymin + dy
%         . . .
%         yf(m)   = ymax - dy
%         yf(m+1) = ymax
%         yf(m+2) = ymax
%     yc(j) is the y-location of control volume (or cell center) j,
%       with an extra entry at ymin and an extra entry at ymax
%         yc(1)   = ymin
%         yc(2)   = ymin + dy/2
%         yc(3)   = ymin + 3*dy/2
%         . . .
%         yc(m)   = ymax - 3*dy/2
%         yc(m+1) = ymax - dy/2
%         yc(m+2) = ymax
% --- indexing conventions ---

  for i=1:nx-1
    xf(i) = xmin + (i-1)*dx;
  end
  i     = nx;
  xf(i) = xf(i-1);
    
  i     = 1;
  xc(i) = xmin;
  for i=2:nx-1
    xc(i) = 0.5*( xf(i) + xf(i-1) );
  end
  i     = nx;
  xc(i) = xmax;
    
  for j=1:ny-1
    yf(j) = ymin + (j-1)*dy;
  end
  j     = ny;
  yf(j) = yf(j-1);
    
  j     = 1;
  yc(j) = ymin;
  for j=2:ny-1
    yc(j) = 0.5*( yf(j) + yf(j-1) );
  end  
  j     = ny;
  yc(j) = ymax;

endfunction

%%% end of grid_2d %%%

%%% calculate geometric interpolation factors %%%

function ginterp
% compute face-center linear interpolation coefficients

% global variables needed by this function
  global xf yf xc yc
  global n m
  
% global variables set by this function
  global gfacx gfacy

  i        = 1;
  gfacx(i) = 0.;
  for i=2:n+1
    gfacx(i) = ( xf(i) - xc(i) ) / ( xc(i+1) - xc(i) );
  end
  i        = n + 2;
  gfacx(i) = 0.;
        
  j        = 1;
  gfacy(j) = 0.;
  for j=2:m+1
    gfacy(j) = ( yf(j) - yc(j) ) / ( yc(j+1) - yc(j) );
  end  
  j        = m + 2;
  gfacy(j) = 0.;

endfunction    
    
%%% end of ginterp %%%

%%% compute cell-face mass flowrates %%%

function facemdots
% compute face mass flow rates corresponding to the given value of the
%   fluid mass density and the presribed velocity field

% global variables needed by this function
  global den xf yf
  global n m
  
% global variables set by this function
  global mdotuf mdotvf

% mdotuf(i,j) is the mass per unit time leaving cell(i,j) across the east face
%   mdotuf = den*u*dy, with u=x
% mdotvf(i,j) is the mass per unit time leaving cell(i,j) across the north face
%   mdotvf = den*v*dx, with v=-y

  for i=1:n+1
    for j=2:m+1
      ufac        = xf(i);    % face-centered x velocity component
      denvel      = den*ufac; % face-centered mass flux out across east face
      mdotuf(i,j) = denvel*( yf(j) - yf(j-1) );
    end
  end  
        
  for j=1:m+1
    for i=2:n+1
      vfac        = -yf(j);   % face-centered y velocity component
      denvel      = den*vfac; % face-centered mass flux out across north face
      mdotvf(i,j) = denvel*( xf(i) - xf(i-1) );
    end
  end

endfunction  

%%% end of facemdots %%%

%%% compute maximum cell Peclet number magnitudes in computational domain %%%

function cellPecletnumbers
% compute maximum cell Peclet number magnitudes in x (based on local x velocity 
%   component and dx) and in y (based on local y velocity component and dy)

% global variables needed by this function
  global den dif dx dy
  global xf yf
  global n m
  
% global variables set by this function
  global pedx_max pedy_max

  pedx_max = 0.;
  pedy_max = 0.;
  
  for i=2:n+1
    ucel     = 0.5*( xf(i-1) + xf(i) ) ; % average x velocity for the cell
    pedx     = abs( den*ucel*dx / dif );
    pedx_max = max( pedx_max , pedx );
    for j=2:m+1
      vcel     = -0.5*( yf(j-1) + yf(j) ); % average y velocity for the cell
      pedy     = abs( den*vcel*dy / dif );
      pedy_max = max( pedy_max , pedy );
    end
  end  

endfunction  

%%% end of cellPecletnumbers %%%

%%% impose Dirichlet BC's on the dependent variable %%%

function dirichletbcs
% set phi(i,j) to presribed boundary values on Dirichlet boundaries

% global variables needed by this function
  global ymin ymax yc
  global n m
  
% global variables set by this function
  global phi
  
% --- indexing conventions ---
% phi(i,j) for 2<=i<=n+1 and 2<=j<=m+1 are the cell-centered values of phi
%   for which we are solving
% phi(1,j) for 2<=j<=m+1 are the x=xmin cell-face-center values of phi,
%   which are prescribed (Dirichlet BC)
% phi(n+2,j) for 2<=j<=m+1 are the x=xmax cell-face-center values of phi;
%   here these are not used, as a Neumann BC is prescribed at x=xmax
% phi(i,1) for 2<=i<=n+1 are the y=ymin cell-face-center values of phi;
%   here these are not used, as a Neumann BC is prescribed at y=ymin
% phi(i,m+2) for 2<=i<=n+1 are the y=ymax cell-face-center values of phi,
%   which are prescribed (Dirichlet BC)
% --- indexing conventions ---

% Dirichlet BC at x=xmin
%   phi varies linearly from 0. at y=ymax to 1. at y=ymin
  i = 1;
  for j=2:m+1
    phi(i,j) = 1. - ( yc(j) - ymin ) / ( ymax - ymin );
  end
  phi(1,1) = 1.;
        
% Dirichlet BC at y=ymax
%   phi=0.
  j = m + 2;
  for i=2:n+1
    phi(i,j) = 0.;
  end    

endfunction

%%% end of dirichletbcs %%%

%%% assemble the coefficient matrix, without consideration of BCs %%%   
        
function assemble_matrix%%% end of assemble_matrix %%%
% set up to five non-zero coefficients (awe, aso, ap, ano, aea) and right-hand
% side values (q) for each interior finite-volume cell (i,j),%%% modify coefficient matrix and rhs to account for BCs %%%
% without considering the BCs
% use UDS for convection for cflag=1function apply_bcs
% use CDS for convection for cflag=2% modify matrix coeffients (awe, aso, ap, ano, aea) and right-hand-side
% use CDS for diffusion in all cases

% global variables needed by this function
global n m
global cflagendfunction  
global dif xf yf xc yc mdotuf mdotvf gfacx gfacy
%%% end of apply_bcs %%%
% global work array used by this function
global dw%%% convert from 2D indexing to 1D indexing %%%

% global variables set by this functionfunction from_2d_to_1d
global awe aso ap ano aea q% convert from 2D indexing to 1D indexing in preparation for solving
%   the linear system
% initialize work array dw(j)
% dw(j) is used to hold east/west face diffusive mass flow rates% global variables needed by this function
% note that xc(2) - xc(1) = dx/2  global n m

for j=2:m+1
dw(j) = -dif*( yf(j) - yf(j-1) ) / ( xc(2) - xc(1) );  
end
  
for i=2:n+1
% note that yc(2) - yc(1) = dy/2  for i=2:n+1
asdif = -dif*( xf(i) - xf(i-1) ) / ( yc(2) - yc(1) );    

for j=2:m+1   
awdif = dw(j);        
aedif = -dif*( yf(j) - yf(j-1) ) / ( xc(i+1) - xc(i) ); 
andif = -dif*( xf(i) - xf(i-1) ) / ( yc(j+1) - yc(j) );

% select UDS or CDS for convection      
if (cflag==1) % UDS 
awcon = -max( mdotuf(i-1,j) , 0. );      
ascon = -max( mdotvf(i,j-1) , 0. );      
ancon = min( mdotvf(i,j) , 0. );      
aecon = min( mdotuf(i,j) , 0. );      
elseif (cflag==2) % CDS      
awcon = -mdotuf(i-1,j)*( 1. - gfacx(i-1) );
ascon = -mdotvf(i,j-1)*( 1. - gfacy(j-1) );
ancon = mdotvf(i,j)*gfacy(j); 
aecon = mdotuf(i,j)*gfacx(i);            
endif 
        
% set coefficients and right-hand-side vector 
awe(i,j) = awdif + awcon;
aso(i,j) = asdif + ascon;        
ano(i,j) = andif + ancon;
aea(i,j) = aedif + aecon;

ap(i,j) = -( awe(i,j) + aso(i,j) + ano(i,j) + aea(i,j) );

q(i,j) = 0.; % physical source term is zero here

dw(j) = aedif;
asdif = andif;

end % end of loop over j

end % end of loop over i

endfunction      
%%% end of assemble_matrix %%%% global work array used by this function

  global ipvt
%%% modify coefficient matrix and rhs to account for BCs %%%  
% global variables changed by this function
function apply_bcs  global a1d
% modify matrix coeffients (awe, aso, ap, ano, aea) and right-hand-side
% values (q) for each finite-volume cell (i,j) adjacent to one or more  info  = 0;
% boundaries, to account for the BCs  nmm1  = nm - 1;
% the treatment for Neumann BCs is formally second-order accurate:    
% phi_e = phi_P + 0.5*dx*dphi/dx along east face at x = xmax  for i=1:nmm1
% phi_s = phi_P - 0.5*dy*dphi/dy along south face at y = ymin    mi  = i;
% the treatment for Dirichlet BCs is formally first-order accurate:    ip1 = i + 1;
% dphi/dx = 2*( phi_P - phi_w ) / dx along west face at x = xmin
% dphi/dy = 2*( phi_n - phi_P ) / dy along north face at y = ymax    for k=ip1:nm
      if (abs(a1d(k,i))>abs(a1d(mi,i)))
% global variables needed by this function        mi = k;
global n m      endif 
global phi    end 

% global variables changed by this function    ipvt(i) = mi;
global awe aso ap ano aea q    
    if (a1d(mi,i)==0.)
% along x=xmin=0: Dirichlet BC      info = i;
% boundary values are in phi(1,j)    else
% no adjustment to ap(2,j) is required since u_w = 0 here (zero convective      if (mi!=i)
% flux across the boundary) and the 1/2-cell dx was accounted for        for k=i:nm
% in computing the original coefficient for diffusion          t         = a1d(i,k);
for j=2:m+1          a1d(i,k)  = a1d(mi,k);
q(2,j) = q(2,j) - awe(2,j)*phi(1,j);          a1d(mi,k) = t;
awe(2,j) = 0.;        end
      endif  
      rp = -1. / a1d(i,i);
      for k=ip1:nm
        a1d(k,i) = a1d(k,i)*rp;
      end
      for j=ip1:nm
        c = a1d(i,j);
        if (c!=0.);
          for k=ip1:nm
            a1d(k,j) = a1d(k,j) + a1d(k,i)*c;
          end  
        endif
      end
    endif
  end % end of loop over i  
      
  ipvt(nm) = nm;
  if (a1d(nm,nm)==0.)
    info = nm;
  endif

% a non-zero value of info indicates that a zero was encountered along the
%   diagonal in the course of the Gauss elimination; this is not necessarily
%   a problem here
  info;  
    
endfunction

%%% end of lufactor %%%

%%% solve the linear system via back substitution %%%

function lusolve
% on entry, a1d holds the LU decomposition of the original coefficient matrix,
%   q1d holds the right-hand-side vector, and ipvt holds the pivot vector
%   from the LU decomposition
% on exit, q1d holds the solution to the linear system

% global variables needed by this function
  global nm
  global a1d
  
% global work array used by this function
  global ipvt
     
% global variables changed by this function
  global q1d

  nmm1 = nm - 1;

  for i=1:nmm1
    mi = ipvt(i);
    c  = q1d(mi);
    if (mi!=i)
      q1d(mi) = q1d(i);
      q1d(i)  = c;
    endif  
    if (c!=0.)
      ip1 = i + 1;
      for k=ip1:nm
        q1d(k) = q1d(k) + a1d(k,i)*c;
      end
    endif
  end % end of loop over i  

  for i=1:nmm1
    j      = nm - i + 1;
    q1d(j) = q1d(j) / a1d(j,j);
    c      = q1d(j);
    if (c!=0.)
      jm1 = j - 1;
      for k=1:jm1
        q1d(k) = q1d(k) - a1d(k,j)*c;
      end
    endif
  end % end of loop over i  

  q1d(1) = q1d(1) / a1d(1,1);

endfunction

%%% end of lusolve %%%

%%% convert the solution from 1D indexing back to 2D indexing %%% 

function from_1d_to_2d
% on entry, q1d holds the solution to the linear system using 1D indexing
% on exit, phi holds the solution to the linear system using 2D indexing

% global variables needed by this function
  global n m
  global q1d phi1d
    
% global variables set by this function
  global phi

  phi1d = q1d;

  for i=2:n+1
    for j=2:m+1
      k       = ( i - 2 )*m + j - 1;
      phi(i,j) = phi1d(k);
    end
  end  

endfunction
        
%%% end of from_1d_to_2d %%%

%%% impose BCs on the computed dependent variable array %%%

function postbc
% reimpose BCs on the computed dependent variable array
% this is mainly for plotting purposes

% global variables needed by this function
  global n m
    
% global variables changed by this function
  global phi

% dphi/dy=0. along y = ymin
  for i=2:n+1
    phi(i,1) = phi(i,2);
  end  
    
% dphi/dx=0. along x = xmax
  for j=2:m+1
    phi(n+2,j) = phi(n+1,j);
  end
    
% Dirichlet BC along x = xmin (lower left-hand corner)
  phi(1,1) = 1.;

% Dirichlet BC along y = ymax (upper right-hand corner)
  phi(n+2,m+2) = 0.;

endfunction

%%% end of postbc %%%

%%% calculate flowrates across boundaries of the computational domain %%%

function bndryflowrates
% compute mdot_xmin, cdot_xmin, ddot_xmin, mdot_xmax, cdot_xmax, ddot_xmax,
%         mdot_ymin, cdot_ymin, ddot_ymin, mdot_ymax, cdot_ymax, ddot_ymax,
%         mdot_tot,  cdot_tot,  ddot_tot,  phidot_tot

% mdot_BNDRY is the total mass flow rate across the boundary
%   the units of mdot_BNDRY are the same as those of mdotuf and mdotvf
% cdot_BNDRY is the convective flow rate of phi across the boundary
%   the units of cdot_BNDRY are the units of mdot_BNDRY times phi
% ddot_BNDRY is the diffusive flow rate of phi across the boundary
%   the units of ddot_BNDRY are the units of mdot_BNDRY times phi

% mdot_tot is the sum of the values of mdot_BNDRY over all boundaries
% cdot_tot is the sum of the values of cdot_BNDRY over all boundaries
% ddot_tot is the sum of the values of ddot_BNDRY over all boundaries
% phidot_tot is the sume of cdot_tot and ddot_tot

% flows out of the domain are positive, flows into the domain are negative

% it is essential to use discrete approximations that are consistent with those
%   that were used in the discrete algebraic system:
%     midpoint rule for surface integrals
%     first-order approximation for diffusive fluxes across Dirichlet boundaries
%     first- (UDS convection) or second- (CDS convection) order approximations
%       for convective fluxes across Neumann boundaries

endfunction
        
%%% end of bndryflowrates %%%

%%% generate a 2d contour plot of phi %%%

function phi2dcontourplot
% generate and save a 2d contour plot of phi
% this version works only for n=m 

% global variables needed by this function
  global xc yc phi
  
% global variables set by this function
% none

  figure;
  contourf(xc,yc,phi);
  colorbar;
  xlabel('x');
  ylabel('y');
  title('Phi versus x and y');
  saveas(gcf,'phi_2dcontours.png');

endfunction

%%% end of phi2dcontourplot %%%

%%% write an output file %%%

function outfile
% write an ascii text output file
% echoes input parameters and key derived quantitites

% global variables needed by this function
  global xmin xmax ymin ymax den dif
  global n m
  global cflag
  global pedx_max pedy_max
  global lufactor_cpu lusolve_cpu
  global mdot_xmin mdot_xmax mdot_ymin mdot_ymax mdot_tot
  global cdot_xmin cdot_xmax cdot_ymin cdot_ymax cdot_tot
  global ddot_xmin ddot_xmax ddot_ymin ddot_ymax ddot_tot
  global phidot_tot
    
% global variables set by this function
% none
  
% open the output file
  fileID = fopen('output.txt','w');  
  
% echo input parameters
  fprintf(fileID,'xmin,xmax          = %12f %12f \n', ...
                  xmin,xmax);
  fprintf(fileID,'ymin,ymax          = %12f %12f \n', ...
                  ymin,ymax);
  fprintf(fileID,'den,dif            = %12f %12f \n', ...
                  den,dif);
  fprintf(fileID,'n,m                = %12i %12i \n', ...
                  n,m);
  fprintf(fileID,'cflag              = %12i      \n', ...
                  cflag);
  fprintf(fileID,'\n');
  
% write key derived quantities
  fprintf(fileID,'pedx_max,pedy_max  = %12f %12f \n', ...
                  pedx_max,pedy_max);
  fprintf(fileID,'\n');
  
% write CPU times (s) for key operations
  fprintf(fileID,'lufactor_cpu,lusolve_cpu = %12f %12f \n', ...
                  lufactor_cpu,lusolve_cpu);
  fprintf(fileID,'\n');
                  
% write boundary mass flow rates
  fprintf(fileID,'mdot_xmin,mdot_xmax,mdot_ymin,mdot_ymax = %+15e %+15e %+15e %+15e \n', ...              
                  mdot_xmin,mdot_xmax,mdot_ymin,mdot_ymax);
                  
% write boundary convective flow rates of phi
  fprintf(fileID,'cdot_xmin,cdot_xmax,cdot_ymin,cdot_ymax = %+15e %+15e %+15e %+15e \n', ...              
                  cdot_xmin,cdot_xmax,cdot_ymin,cdot_ymax);
                  
% write boundary diffusive flow rates of phi
  fprintf(fileID,'ddot_xmin,ddot_xmax,ddot_ymin,ddot_ymax = %+15e %+15e %+15e %+15e \n', ...              
                  ddot_xmin,ddot_xmax,ddot_ymin,ddot_ymax);
                  
% write net flow rates across all boundaries
  fprintf(fileID,'mdot_tot,cdot_tot,cdot_tot,phidot_tot   = %+15e %+15e %+15e %+15e \n', ...              
                  mdot_tot,cdot_tot,cdot_tot,phidot_tot);
                                   
% close the output file  
  fclose(fileID);

endfunction

%%% end of outfile %%%
     
%%%%% execute the program %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate the grid
grid_2d;

% calculate geometric interpolation coefficients
ginterp;

% compute face mass flow rates
facemdots;

% compute maximum cell Peclet number magnitudes in the computational domain
cellPecletnumbers;

% impose Dirichlet BCs on dependent variable
dirichletbcs;

% assemble the coefficient matrix, without consideration of BCs
assemble_matrix;
%awe
%aso
%ap
%ano
%aea
%q

% modify the coefficient matrix and right-hand side to account for BCs
apply_bcs;
%awe
%aso
%ap
%ano
%aea
%q

% convert from 2D indexing to 1D indexing
from_2d_to_1d;

% factor the coefficient matrix, and track the CPU time required
t0 = clock();
lufactor;
lufactor_cpu = etime (clock (), t0);

% solve the linear system, and track the CPU time required
t0 = clock();
lusolve;
lusolve_cpu = etime (clock (), t0);

% convert from 1D indexing to 2D indexing
from_1d_to_2d;

% impose BCs on the computed dependent variable array
postbc;

% generate a 2D contour plot of phi
phi2dcontourplot;

% compute flowrates across boundaries of the computational domain
bndryflowrates;

% write an output file
outfile;
 
%%%%% all done %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%