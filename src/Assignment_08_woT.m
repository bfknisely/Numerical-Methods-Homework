%%%%%%%%
% ME523 Spring 2018
% Dan Haworth, March 2018
%%%%%%%%

% This program solves the 2D steady or unsteady Navier-Stokes equations
%   in a rectangular domain, using a colocated pressure-based
%   finite-volume method and a Cartesian grid with constant grid spacing
%   in x and constant grid spacing in y.

%  It is set up for lid- and buoyancy-driven flows in closed cavities
%    using a blended UDS/CDS scheme for convective fluxes, 
%    CDS for diffusive fluxes,
%    and implicit Euler time advancement. 

% Adapted from a 2d Navier-Stokes solver "pcol" by Ferziger & Peric,
%   3rd Edition, 2002.

% Echo all user specified quantities to an output file, along with results

% Not tested:
%   buoyancy source term and calct(beta,gravx,gravy)

%%%%% define global variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% scalars: physical quantities
global xmin;    % minimum value of independent variable x (m)
global xmax;    % maximum value of independent variable x (m)
global ymin;    % minimum value of independent variable y (m)
global ymax;    % maximum value of independent variable y (m)

global densit;  % fluid mass density (kg/m^3)
global visc;    % fluid dynamic viscosity (kg/m-s)
global prm;     % fluid Prandtl number (-)
global prr;     % 1/prm
global beta;    % fluid volumetric coefficient of thermal expansion (1/K)

global gravx;   % body force vector in x direction (m/s^2)
global gravy;   % body force vector in y direction (m/s^2)

global th;      % temperature of hot wall at x=xmin (K)
global tc;      % temperature of cold wall at x=xmax (K)
global tref;    % reference temperature (K)
global ulid;    % x-direction lid speed for top wall at y=ymax (m/s)

global uin;     % uniform initial x velocity (m/s)
global vin;     % uniform initial y velocity (m/s)
global pin;     % uniform initial pressure (Pa)
global tin;     % uniform initial temperature (K)

% scalars: grid-related quantities - see grid_2d below
global ncvx;    % number of control volumes in the x direction
global ncvy;    % number of control volumes in the y direction

global ni;      % ni  = ncvx + 2
global nx;      % nx  = ni
global nim;     % nim = ni - 1 = ncvx + 1

global nj;      % nj  = ncvy + 2
global ny;      % ny  = nj
global njm;     % njm = nj - 1 = ncvy + 1

global ncv;     % ncv = ncvx*ncvy
global nij;     % nij = ni*nj
global nxy;     % nxy = nx*ny

global dx;      % grid spacing in x direction (constant)
global dy;      % grid spacing in y direction (constant)

% scalars: numerical solution control parameters
global nphi;    % number of dependent variables for which to solve
global iu;      % dependent variable corresponding to x velocity component
global iv;      % dependent variable corresponding to y velocity component
global ip;      % dependent variable corresponding to pressure
global ien;     % dependent variable corresponding to energy (temperature)

global small;   % a small positive real number
global great;   % a large positive real number

global ltime;   % flag to select steady or unsteady problem
                %   steady:   ltime=0
                %   unsteady: ltime=1
global dt;      % constant computational time step for ltime=1 (s)
global dtr;     % 1/dt
global itst;    % number of computational time steps for ltime=1
                %   itst=1 for ltime=0                
global maxit;   % maximum number of outer iterations
                %   per time step, for ltime=1
                
global ipr;     % i (x) index of node at which the pressure is kept fixed
                %   1<=ipr<=ncvx
global jpr;     % j (y) index of node at which the pressure is kept fixed
                %   1<=jpr<=ncvy
                % ipr,jpr define the reference location for pressure
                
global sormax;  % level of residual norm at which outer iterations are stopped
                %   (a convergence criterion)
global slarge;  % level of residual norm at which iterations are stopped  
                %   because of divergence
                
global alfa;    % parameter for the SIP linear system solver

% scalars: output control
global nprt;    % number of computational time steps between writing output
                %   for ltime=1
global imon;    % i (x) index of monotoring location control volume
                %   1<=imon<=ncvx
global jmon;    % j (y) index of monotoring location control volume
                %   1<=jmon<=ncvy
global ijmon;   % monitoring cell id in 1D indexing

global sipfile  % output file to monitor linear system solver convergence
global resfile; % output file to monitor outer iteration convergence
global monfile; % output file to monitor dependent variable solution               
                
% arrays (vectors) for cell face and cell center locations
global xf;      % cell face (coordinate line) x coordinates (nx x 1)
global yf;      % cell face (coordinate line) y coordinates (ny x 1)
global xc;      % cell center x coordinates (nx x 1)
global yc;      % cell center y coordinates (ny x 1)

% arrays (vectors) for cell-face geometric interpolation factors
global gfacx;   % in x direction (nx x 1)
global gfacy;   % in y direction (ny x 1)

% arrays (vectors) for cell-face mass flow rates
global f1;      % across east  face (nxy x 1)
global f2;      % across north face (nxy x 1)
                           
% arrays (vectors) for numerical solution control parameters (all are nphi x 1)
global lcal;    % flag to solve or not for each dependent variable
                %   set to 1 to solve, set to 0 to bypass
global urf;     % under-relaxation factor for each dependent variable
                %   must be between 0. and 1.
global sor;     % required ratio of reduction of residual norm during
                %   inner iterations for each dependent variable                
global nsw;     % max number of inner iterations for each dependent variable
global gds;     % UDS/CDS blending factor for each dependent variable
                %   0. for UDS
                %   1. for CDS
                %   between 0. and 1. for UDS/CDS blend

% arrays (vectors) holding nodal values of the dependent variables, plus
%   additional values for BCs (all nxy x 1)
global u;       % x velocity component
global v;       % y velocity component
global p;       % pressure
global t;       % temperature

% arrays (vectors) holding end-of-last-time-step values of dependent variables
%   (all nxy x 1)
global u0;      % x velocity component
global v0;      % y velocity component
global p0;      % pressure
global t0;      % temperature

% array (vector) holding pressure correction
global pp;      % (nxy x 1)

% arrays (vectors) holding coefficients and right-hand sides (all nxy x 1)
global su;      % a right-hand side
global sv;      % a right-hand side
global apu;     % a coefficient in momentum equation
global apv;     % a coefficient in momentum equation
global awe;     % coefficient for control volume west  of control volume P
global aso;     % coefficient for control volume south of control volume P
global ap;      % diagonal coefficient for control volume P
global ano;     % coefficient for control volume north of control volume P
global aea;     % coefficient for control volume east  of control volume P
global dpx;     % x-direction pressure gradient
global dpy;     % y-direction pressure gradient

% work arrays (vectors) used in the SIP linear system solver (all nxy x 1)
global lw;      % 
global ls;      % 
global lpr;     % 
global un;      % 
global ue;      % 
global res;     % 
global fi;      % 

global resor;   % 
                
% array (vector) to facilitate 2D-1D indexing conversions
global li;      % li(i) holds pointer to last location in column i-1 (nx x 1)

%%%%% define the physical problem, grid, and numerical control parameters %%%%%%

xmin    = 0.;
xmax    = 1.;
ymin    = 0.;
ymax    = 1.;

densit  = 1.;
visc    = 0.001;
prm     = 0.1;
beta    = 0.0;

gravx   = 0.;
gravy   = 0.;

th      = 400.;
tc      = 300.;
tref    = 300.;

ulid    = 1.;

uin     =   0.;
vin     =   0.;
pin     =   0.; % pressure difference with respect to an arbitrary reference
                %   pressure - not an absolute thermodynamic pressure
tin     = 300.;

nphi    = 4;

ncvx    =  20;
ncvy    =  20;

ltime   = 0;
dt      = 0.1;
itst    = 1;
nprt    = 1;

ipr     = 3;
jpr     = 3;

maxit   = 500;
sormax  = 0.0001;
slarge  = 1000.;

alfa    = 0.92;

imon    = 5;
jmon    = 5;

lcal    = zeros(nphi,1);
lcal(1) = 1;  
lcal(2) = 1;  
lcal(3) = 1;  
lcal(4) = 1;  

urf     = zeros(nphi,1); 
urf(1)  = 0.8;
urf(2)  = 0.8;
urf(3)  = 0.2;
urf(4)  = 0.8;

sor     = zeros(nphi,1); 
sor(1)  = 0.2;
sor(2)  = 0.2;
sor(3)  = 0.2;
sor(4)  = 0.2;

nsw     = zeros(nphi,1);
nsw(1)  = 1;
nsw(2)  = 1;
nsw(3)  = 6;
nsw(4)  = 2;

gds     = zeros(nphi,1);
gds(1)  = 1.;
gds(2)  = 1.;
gds(3)  = 1.;
gds(4)  = 1.;

%%%%% compute derived scalar variables for indexing and array dimensioning %%%%%

ni  = ncvx + 2;
nx  = ni;
nim = ni - 1;
nj  = ncvy + 2;
ny  = nj;
njm = nj -1;

ncv = ncvx*ncvy;
nij = ni*nj;
nxy = nx*ny;

% Re based on ulid, other dimensionless groups - or do this in a function

%%%%% define and initialize arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% all array element values are initialized to zero

xf    = zeros(nx,1);      
yf    = zeros(ny,1);       
xc    = zeros(nx,1);      
yc    = zeros(ny,1);      

gfacx = zeros(nx,1);      
gfacy = zeros(ny,1);     

f1    = zeros(nxy,1);      
f2    = zeros(nxy,1);      
      
u     = zeros(nxy,1);      
v     = zeros(nxy,1);       
p     = zeros(nxy,1);       
t     = zeros(nxy,1);      

u0    = zeros(nxy,1);     
v0    = zeros(nxy,1);      
p0    = zeros(nxy,1);     
t0    = zeros(nxy,1);      

pp    = zeros(nxy,1);     

su    = zeros(nxy,1);      
sv    = zeros(nxy,1);       
apu   = zeros(nxy,1);      
apv   = zeros(nxy,1);     
awe   = zeros(nxy,1);   
aso   = zeros(nxy,1); 
ap    = zeros(nxy,1);    
ano   = zeros(nxy,1);    
aea   = zeros(nxy,1);   
dpx   = zeros(nxy,1);     
dpy   = zeros(nxy,1);     

lw    = zeros(nxy,1);     
ls    = zeros(nxy,1);     
lpr   = zeros(nxy,1);     
un    = zeros(nxy,1);     
ue    = zeros(nxy,1);     
es    = zeros(nxy,1);     
fi    = zeros(nxy,1);
res   = zeros(nxy,1);       

resor = zeros(nphi,1);   
                
li    = zeros(nx,1);      

%%%%% define functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% set up the 2d grid %%%

function grid_2d
% compute cell-center and cell-face-center locations, geometric
%   interpolation coefficients, array for 1D-2D indexing conversions,
%   and ijmon

% here the grid spacings in x (dx) and in y (dy) are each constant,
%   but dx and dy are not necessarily equal

% global variables needed by this function
  global xmin xmax ymin ymax
  global ncvx ncvy nx ny
  global imon jmon
  
% global variables set by this function
  global dx dy
  global xf yf xc yc
  global gfacx gfacy
  global li
  global ijmon

% --- begin indexing conventions ---

%   in the x-direction:
%     ncvx is the number of cells or control volumes in the x direction
%     ni   = ncvx + 2
%     nx   = ni
%     nim  = ni - 1 = ncvx + 1
%     xf(i) is the x-location of grid-line (or face) i, 
%       with an extra entry at xmax
%         xf(1)    = xmin
%         xf(2)    = xmin + dx
%         . . .
%         xf(ni-2) = xmax - dx
%         xf(ni-1) = xmax
%         xf(ni)   = xmax
%     xc(i) is the x-location of control volume (or cell center) i,
%       with an extra entry at xmin and an extra entry at xmax
%         xc(1)    = xmin
%         xc(2)    = xmin + dx/2
%         xc(3)    = xmin + 3*dx/2
%         . . .
%         xc(ni-2) = xmax - 3*dx/2
%         xc(ni-1) = xmax - dx/2
%         xc(ni)   = xmax

%   in the y-direction:
%     ncvy is the number of cells or control volumes in the y direction
%     nj   = ncvy + 2
%     ny   = nj
%     njm  = nj - 1 = ncvy + 1
%     yf(j) is the y-location of grid-line (or face) j,
%       with an extra entry at ymax
%         yf(1)    = ymin
%         yf(2)    = ymin + dy
%         . . .
%         yf(nj-2) = ymax - dy
%         yf(nj-1) = ymax
%         yf(nj)   = ymax
%     yc(j) is the y-location of control volume (or cell center) j,
%       with an extra entry at ymin and an extra entry at ymax
%         yc(1)    = ymin
%         yc(2)    = ymin + dy/2
%         yc(3)    = ymin + 3*dy/2
%         . . .
%         yc(nj-2) = ymax - 3*dy/2
%         yc(nj-1) = ymax - dy/2
%         yc(nj)   = ymax

%   1D indexing
%     arrays are indexed from 1 to (ncvx+2)*(ncvy+2)
%     there are entries for the ncvx*ncvy nodal (or cell or
%       control volume) values of the dependent variables for which we are
%       solving, plus entries to accommodate boundary conditions 

% --- end indexing conventions ---
 
% grid spacings in x and in y

  dx  = ( xmax - xmin ) / ncvx;
  dy  = ( ymax - ymin ) / ncvy;

% xf and xc

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
  
% yf and yc
    
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
  
% gfacx and gfacy  

  i        = 1;
  gfacx(i) = 0.;
  for i=2:ncvx+1
    gfacx(i) = ( xf(i) - xc(i) ) / ( xc(i+1) - xc(i) );
  end
  i        = ncvx + 2;
  gfacx(i) = 0.;
        
  j        = 1;
  gfacy(j) = 0.;
  for j=2:ncvy+1
    gfacy(j) = ( yf(j) - yc(j) ) / ( yc(j+1) - yc(j) );
  end  
  j        = ncvy + 2;
  gfacy(j) = 0.;
  
% li

  for i=1:nx
    li(i) = ( i - 1 )*ny;
  end  
  
% ijmon  
  
  ijmon = li(imon) + jmon;

endfunction

%%% end of grid_2d %%%

%%% set internal solution control parameters and variables %%%

function solnparm
% compute cell-center and cell-face-center locations, geometric

% global variables needed by this function
  global dt prm
  
% global variables set by this function
  global iu iv ip ien
  global small great
  global dtr prr

  iu  = 1; % dependent variable corresponding to x velocity component
  iv  = 2; % dependent variable corresponding to y velocity component
  ip  = 3; % dependent variable corresponding to pressure
  ien = 4; % dependent variable corresponding to energy (temperature)
  
  small = 1.e-15; % a very small real number
  great = 1.e+15; % a very large real number
  
  dtr = 1. / dt;  % inverse of computational time step
  prr = 1. / prm; % inverse of Prandtl number

endfunction

%%% end of solnparm %%%

%%% set boundary conditions for dependent variables %%%

function bcs
% all arrays were initialized to zero, so just set nonzero values here

% global variables needed by this function
    global ni nj li
    global ulid tc th
  
% global variables modified by this function
  global u t

% t = th along x=xmin
  for j=1:nj
    t(j) = th;
  end

% t = tc along x=xmax
  for j=1:nj
    ij    = li(ni) + j;
    t(ij) = tc;
  end    
  
% u = ulid along y=ymax
  for i=1:ni
    ij    = li(i) + nj;
    u(ij) = ulid;
  end  

endfunction

%%% end of bcs %%%

%%% set initial conditions for dependent variables %%%

function ics
% all initial conditions are uniform in the interior of the computational
%   domain
% take care not to overwrite the boundary conditions

% global variables needed by this function
  global ni nj li
  global uin vin pin tin
  
% global variables set or modified by this function
  global u  v  p  t
  global u0 v0 p0 t0

  for i=2:ni-1
    ijstart = li(i) + 2;
    ijend   = li(i) + nj - 1; 
    for ij=ijstart:ijend
      u (ij) = uin;
      v (ij) = vin;
      p (ij) = pin;
      t (ij) = tin;
      
      u0(ij) = uin;
      v0(ij) = vin;
      p0(ij) = pin;
      t0(ij) = tin;
    end  
  end

endfunction

%%% end of ics %%%

%%% solve a linear system using Stone's SIP solver %%%

function sipsol(ifi)
% solve a linear system using Stone's SIP solver

% ifi=1: solve for x velocity component
% ifi=2: solve for y velocity component
% ifi=3: solve for pressure correction
% ifi=4: solve for temperature

% global variables needed by this function
  global nij nim njm nj nxy
  global li
  global nsw sor
  global awe aso ap ano aea su
  global alfa
  global sipfile
  
% global variables set or modified by this function
  global u v pp t
  global lw ls lpr un ue fi
  global ue un res resor
   
% load fi() with dependent variable field for which we are currently solving
  if (ifi==1);
    for ij=1:nij
      fi(ij) = u(ij);
    end
  elseif (ifi==2);
    for ij=1:nij
      fi(ij) = v(ij);
    end
  elseif (ifi==3);
    for ij=1:nij
      fi(ij) = pp(ij);
    end
  elseif (ifi==4);
    for ij=1:nij
      fi(ij) = t(ij);
    end
  endif
   
% initialize arrays
  ue (1:nxy,1) = 0.;  
  un (1:nxy,1) = 0.;
  res(1:nxy,1) = 0.;

% coefficients of upper and lower triangular matrices
  for i=2:nim
    for ij=li(i)+2:li(i)+njm
      lw(ij)  = awe(ij) / ( 1. + alfa*un(ij-nj) );
      ls(ij)  = aso(ij) / ( 1. + alfa*ue(ij-1)  );
      p1      = alfa*lw(ij)*un(ij-nj);
      p2      = alfa*ls(ij)*ue(ij-1);
      lpr(ij) = 1. / ( ap(ij) + p1 + p2 - lw(ij)*ue(ij-nj) - ls(ij)*un(ij-1) );
      un(ij)  = ( ano(ij) - p1 )*lpr(ij);
      ue(ij)  = ( aea(ij) - p2 )*lpr(ij);
    end
  end
  
%%% inner iteration loop
  l   = 0;
  rsm = 2.*sor(ifi);  
  while ( l<nsw(ifi) && rsm>sor(ifi) );
    l    = l + 1;
    resl = 0.;

% calculate residual and overwrite it by intermediate vector      
    for i=2:nim
      for ij=li(i)+2:li(i)+njm
        res(ij) = su(ij) - ano(ij)*fi(ij+1) - aso(ij)*fi(ij-1) ...
                - aea(ij)*fi(ij+nj) - awe(ij)*fi(ij-nj) - ap(ij)*fi(ij);
        resl    = resl + abs( res(ij) );
        res(ij) = ( res(ij) - ls(ij)*res(ij-1) - lw(ij)*res(ij-nj) )*lpr(ij);
      end
    end

% store initial residual sum for checking convergence of outer iterations
    if(l==1);
      resor(ifi) = resl;
    endif  
    rsm = resl / ( resor(ifi) + 1.e-20 );

% back substitution and correction
    for i=nim:-1:2
      for ij=li(i)+njm:-1:li(i)+2
        res(ij) = res(ij) - un(ij)*res(ij+1) - ue(ij)*res(ij+nj);
        fi(ij)  = fi(ij) + res(ij);
      end
    end

  endwhile % end of inner iteration loop

  fprintf(sipfile,'%12i %12i %12i %+15e %+15e \n', ...
                   ifi,l,nsw(ifi),rsm,sor(ifi) );

% fill dependent variable vector for which we are currently solving
  if (ifi==1)
    for ij=1:nij
      u(ij) = fi(ij);
    end
  elseif (ifi==2)
    for ij=1:nij
      v(ij) = fi(ij);
    end
  elseif (ifi==3)
    for ij=1:nij
      pp(ij) = fi(ij);
    end
  elseif (ifi==4)
    for ij=1:nij
      t(ij) = fi(ij);
    end
  endif

endfunction

%%% end of sipsol %%%

%%% set time-dependent BCs %%%

function bctime
% update BCs for current time step

% currently there are no time-varying BCs - lid speed is constant
% here the lid speed could be specified as a function of time, for example

% global variables needed by this function
  global ni nj
  global li
  global ulid
  
% global variables set or modified by this function
  global u
  
% u = ulid along y=ymax
  for i=1:ni
    ij    = li(i) + nj;
    u(ij) = ulid;
  end  
  
endfunction

%%% end of bctime %%%

%%% set boundary pressure or pressure correction %%%

function pbound(pflag)
% set boundary pressure or pressure correction using linear extrapolation
%   from inside

% pflag = 1: set boundary pressure
% pflag = 2: set boundary pressure correction

% global variables needed by this function
  global nij nim njm ni nj
  global li
  global gfacx gfacy
  
% global variables set or modified by this function
  global fi p pp
 
% load fi() with dependent variable field for which we are solving
  if (pflag==1);
    for ij=1:nij
      fi(ij) = p(ij);
    end
  elseif (pflag==2);
    for ij=1:nij
      fi(ij) = pp(ij);
    end
  endif

% ymin and ymax boundaries
  for i=2:nim
    ij     = li(i) + 1;
    fi(ij) = fi(ij+1) + ( fi(ij+1) - fi(ij+2) )*gfacy(2);
    ij     = li(i) + nj;
    fi(ij) = fi(ij-1) + ( fi(ij-1) - fi(ij-2) )*( 1. - gfacy(njm-1) );
  end

% xmin and xmax boundaries
  nj2 = 2*nj;
  for j=2:njm
    ij     = li(1) + j;
    fi(ij) = fi(ij+nj) + ( fi(ij+nj) - fi(ij+nj2) )*gfacx(2);
    ij     = li(ni) + j;
    fi(ij) = fi(ij-nj) + ( fi(ij-nj) - fi(ij-nj2) )*( 1. - gfacx(nim-1) );
  end
  
% load dependent variable field array for which we are solving
  if (pflag==1);
    for ij=1:nij
      p(ij) = fi(ij);
    end
  elseif (pflag==2);
    for ij=1:nij
      pp(ij) = fi(ij);
    end
  endif
 
endfunction

%%% end of pbound %%%

%%% implement boundary conditions for momentum equations %%%

function bcuv
% implement boundary conditions for momentum equations 

% global variables needed by this function
  global nim njm ni nj
  global li
  global xf yf xc yc u v
  global visc
  
% global variables set or modified by this function
  global apu su apv sv
  
% y=ymin: wall shear force in u equation, zero gradient in v equation (dv/dy=0)  
  for i=2:nim
    ij      = li(i) + 2;
    d       = visc*( xf(i) - xf(i-1) ) / ( yc(2) - yc(1) );
    apu(ij) = apu(ij) + d;
    su(ij)  = su(ij) + d*u(ij-1);
  end

% y=ymax: wall shear force in u equation, zero gradient in v equation (dv/dy=0)
  for i=2:nim
    ij      = li(i) + njm;
    d       = visc*( xf(i) - xf(i-1) ) / ( yc(nj) - yc(njm) );
    apu(ij) = apu(ij) + d;
    su(ij)  = su(ij) + d*u(ij+1);
  end

% x=xmin: wall shear force in v equation, zero gradient in u equation (du/dx=0)  
  for j=2:njm
    ij      = li(2) + j;
    d       = visc*( yf(j) - yf(j-1) ) / ( xc(2) - xc(1) );
    apv(ij) = apv(ij) + d;
    sv(ij)  = sv(ij) + d*v(ij-nj);
  end

% x=xmax: wall shear force in v equation, zero gradient in u equation (du/dx=0)  
  for j=2:njm
    ij      = li(nim) + j;
    d       = visc*( yf(j) - yf(j-1) ) / ( xc(ni) - xc(nim) );
    apv(ij) = apv(ij) + d;
    sv(ij)  = sv(ij) + d*v(ij+nj);
  end   
 
endfunction

%%% end of bcuv %%%

%%% implement boundary conditions for energy (temperature) equation %%%

function bct
% implement boundary conditions for temperature equation 

% global variables needed by this function
  global nim njm ni nj
  global li
  global xf yf xc yc
  global visc prr
  
% global variables set or modified by this function
  global t ap su
  
% y=ymin: adiabatic wall (dt/dy=0) - zero flux  
  for i=2:nim
    ij   = li(i) + 1;
    t(ij)= t(ij+1);
  end

% y=ymax: adiabatic wall (dt/dy=0) - zero flux
  for i=2:nim
    ij    = li(i) + nj;
    t(ij) = t(ij-1);
  end

% x=xmin: isothermal wall - nonzero diffusive flux  
  for j=2:njm
    ij     = li(2) + j;
    d      = visc*prr*( yf(j) - yf(j-1) ) / ( xc(2) - xc(1) );
    ap(ij) = ap(ij) + d;
    su(ij) = su(ij) + d*t(ij-nj);
  end

% x=xmax: isothermal wall - nonzero diffusive flux  
  for j=2:njm
    ij     = li(nim) + j;
    d      = visc*prr*( yf(j) - yf(j-1) ) / ( xc(ni) - xc(nim) );
    ap(ij) = ap(ij) + d;
    su(ij) = su(ij) + d*t(ij+nj);
  end   
 
endfunction

%%% end of bct %%%

%%% solve for x and y velocity components (momentum predictor equation) %%%

function calcuv
% solve linearized momentum equation (momentum predictor)

% global variables needed by this function
  global iu iv ien
  global nij nim njm ni nj
  global li
  global urf gds lcal ltime
  global gfacx gfacy
  global f1 f2
  global xf yf xc yc
  global visc beta densit tref gravx gravy
  global dtr gamt
  global p t u0 v0
  
% global variables set or modified by this function
  global aea awe ano aso ap su sv apu apv
  global dpx dpy
  global u v
  
% reciprocals of under-relaxation factors for u and v
  urfu = 1. / urf(iu);
  urfv = 1. / urf(iv);

% set boundary pressure
  pflag = 1;
  pbound(pflag);
  
% initialize temporary variables
  for ij=1:nij
    su(ij)  = 0.;
    sv(ij)  = 0.;
    apu(ij) = 0.;
    apv(ij) = 0.;
  end  
  
%%% surface integrals: east and west CV faces

  for i=2:nim-1
    fxe  = gfacx(i);
    fxp  = 1. - fxe;
    dxpe = xc(i+1) - xc(i);
    
    for j=2:njm
      ij  = li(i) + j;
      ije = ij + nj;
      
% face area
      s = yf(j) - yf(j-1);
      
% diffusive flux coefficient
      d = visc*s / dxpe;

% explicit convective fluxes for uds and cds
      ce = min( f1(ij) , 0. );
      cp = max( f1(ij) , 0. );

      fuuds = cp*u(ij) + ce*u(ije);
      fvuds = cp*v(ij) + ce*v(ije);
      fucds = f1(ij)*( u(ije)*fxe + u(ij)*fxp );
      fvcds = f1(ij)*( v(ije)*fxe + v(ij)*fxp );      

      aea(ij)  =  ce - d;
      awe(ije) = -cp - d;
      
      su(ij)  = su(ij)  + gds(iu)*( fuuds - fucds );
      su(ije) = su(ije) - gds(iu)*( fuuds - fucds );
      sv(ij)  = sv(ij)  + gds(iu)*( fvuds - fvcds );
      sv(ije) = sv(ije) - gds(iu)*( fvuds - fvcds );
      
    end % end of loop over j
  end % end of loop over i

%%% surface integrals: north and south CV faces

  for j=2:njm-1
    fyn  = gfacy(j);
    fyp  = 1. - fyn;
    dypn = yc(j+1) - yc(j);
    
    for i=2:nim
      ij  = li(i) + j;
      ijn = ij + 1;
      
% face area
      s = xf(i) - xf(i-1);
      
% diffusive flux coefficient
      d = visc*s / dypn;

% explicit convective fluxes for uds and cds
      cn = min( f2(ij) , 0. );
      cp = max( f2(ij) , 0. );

      fuuds = cp*u(ij) + cn*u(ijn);
      fvuds = cp*v(ij) + cn*v(ijn);
      fucds = f2(ij)*( u(ijn)*fyn + u(ij)*fyp );
      fvcds = f2(ij)*( v(ijn)*fyn + v(ij)*fyp );      

      ano(ij)  =  cn - d;
      aso(ijn) = -cp - d;
      
      su(ij)  = su(ij)  + gds(iu)*( fuuds - fucds );
      su(ijn) = su(ijn) - gds(iu)*( fuuds - fucds );
      sv(ij)  = sv(ij)  + gds(iu)*( fvuds - fvcds );
      sv(ijn) = sv(ijn) - gds(iu)*( fvuds - fvcds );
      
    end % end of loop over i
  end % end of loop over j

%%% volume integrals

  for i=2:nim
    dx = xf(i) - xf(i-1);
    for j=2:njm
      dy  = yf(j) - yf(j-1);
      vol = dx*dy;
      ij  = li(i) + j;
      
% pressure terms
      pe = p(ij+nj)*gfacx(i)   + p(ij)   *( 1. - gfacx(i)   );
      pw = p(ij)   *gfacx(i-1) + p(ij-nj)*( 1. - gfacx(i-1) );
      pn = p(ij+1) *gfacy(j)   + p(ij)   *( 1. - gfacy(j)   );
      ps = p(ij)   *gfacy(j-1) + p(ij-1) *( 1. - gfacy(j-1) );
      
      dpx(ij) = ( pe - pw ) / dx;
      dpy(ij) = ( pn - ps ) / dy;
      su(ij)  = su(ij) - dpx(ij)*vol;
      sv(ij)  = sv(ij) - dpy(ij)*vol;
      
% buoyancy term
      if (lcal(ien)==1)
        sb     = -beta*densit*vol*( t(ij) - tref );
        su(ij) = su(ij) + gravx*sb;
        sv(ij) = sv(ij) + gravy*sb;
      endif  
        
% unsteady terms (implicit Euler method)
      if (ltime==1)
        apt     = densit*vol*dtr;
        su(ij)  = su(ij) + apt*u0(ij);
        sv(ij)  = sv(ij) + apt*v0(ij);
        apu(ij) = apu(ij) + apt;
        apv(ij) = apv(ij) + apt;
      endif  
      
    end % end of loop over j
  end % end of loop over i

%%% modify coefficients and right-hand side to account for BCs

  bcuv;

%%% apply under-relaxation for u

  for i=2:nim
    for ij=li(i)+2:li(i)+njm
      ap(ij)  = ( -aea(ij) - awe(ij) - ano(ij) - aso(ij) + apu(ij) )*urfu;
      su(ij)  = su(ij) + ( 1. - urf(iu) )*ap(ij)*u(ij);
      apu(ij) = 1. / ap(ij);
    end
  end
  
%%% solve for u

   sipsol(iu);
 
%%% apply under-relaxation for v

  for i=2:nim
    for ij=li(i)+2:li(i)+njm
      ap(ij)  = ( -aea(ij) - awe(ij) - ano(ij) - aso(ij) + apv(ij) )*urfv;
      su(ij)  = sv(ij) + ( 1. - urf(iv) )*ap(ij)*v(ij);
      apv(ij) = 1. / ap(ij);
    end
  end
  
%%% solve for v

   sipsol(iv);   
     
endfunction

%%% end of calcuv %%%

%%% solve for pressure %%%

function calcp
% solve the pressure correction equation, and apply pressure and velocity
%   corrections
% also set cell face mass flow rates

% global variables needed by this function
  global nim njm ni nj
  global ip ipr jpr
  global xf yf xc yc
  global gfacx gfacy li
  global densit
  global dpx dpy apu apv
  global urf

% global variables set or modified by this function
  global u v p pp
  global awe aso ano aea ap su sv
  global f1 f2
  
%%% east and west faces
  for i=2:nim-1
    dxpe = xc(i+1) - xc(i);
    fxe  = gfacx(i);
    fxp  = 1. - fxe;
    
    for j=2:njm
      ij  = li(i) + j;
      ije = ij + nj;
      
      s    = yf(j) - yf(j-1);
      vole = dxpe*s;
      d    = densit*s;
      
      dpxel = 0.5*( dpx(ije) + dpx(ij) );
      uel   = u(ije)*fxe   + u(ij)*fxp;
      apue  = apu(ije)*fxe + apu(ij)*fxp;
      
      dpxe   = ( p(ije) - p(ij) ) / dxpe;
      ue     = uel - apue*vole*( dpxe - dpxel );
      f1(ij) = d*ue;
      
      aea(ij)  = -d*apue*s;
      awe(ije) = aea(ij);
    end
  end  
      
%%% north and south faces
  for j=2:njm-1
    dypn = yc(j+1) - yc(j);
    fyn  = gfacy(j);
    fyp  = 1. - fyn;
    
    for i=2:nim
      ij  = li(i) + j;
      ijn = ij + 1;
      
      s    = xf(i) - xf(i-1);
      voln = s*dypn;
      d    = densit*s;
      
      dpynl = 0.5*( dpy(ijn) + dpy(ij) );
      vnl   = v(ijn)*fyn   + v(ij)*fyp;
      apvn  = apv(ijn)*fyn + apv(ij)*fyp;
      
      dpyn   = ( p(ijn) - p(ij) ) / dypn;
      vn     = vnl - apvn*voln*( dpyn - dpynl );
      f2(ij) = d*vn;
      
      ano(ij)  = -d*apvn*s;
      aso(ijn) = ano(ij);
    end
  end  
  
%%% source term and coefficient of node P
  sum = 0.;
  for i=2:nim
    for ij=li(i)+2:li(i)+njm
      su(ij)  = f1(ij-nj) - f1(ij) + f2(ij-1) - f2(ij);
      ap(ij)  = -( aea(ij) + awe(ij) + ano(ij) + aso(ij) );
      sum     = sum + su(ij);
      pp(ij)  = 0.;
    end
  end
  
%%% solve for pressure correction
  sipsol(ip);
  
%%% calculate pressure correction at boundaries
  pflag = 2;
  pbound(pflag);
  
%%% subtract value of pressure correction at reference location from all
  ijpref = li(ipr) + jpr;
  ppo    = pp(ijpref);
  
  for i=2:nim-1
    for ij=li(i)+2:li(i)+njm
      f1(ij) = f1(ij) + aea(ij)*( pp(ij+nj) - pp(ij) );
    end
  end

  for i=2:nim
    for ij=li(i)+2:li(i)+njm-1
      f2(ij) = f2(ij) + ano(ij)*( pp(ij+1) - pp(ij) );
    end
  end

  for i=2:nim
    dx = xf(i) - xf(i-1);
    
    for j=2:njm
      ij = li(i) + j;
      dy = yf(j) - yf(j-1);
      
      ppe = pp(ij+nj)*gfacx(i) + pp(ij)*( 1. - gfacx(i) );
      ppw = pp(ij)*gfacx(i-1)  + pp(ij-nj)*( 1. - gfacx(i-1) );
      ppn = pp(ij+1)*gfacy(j)  + pp(ij)*( 1. - gfacy(j) );
      pps = pp(ij)*gfacy(j-1)  + pp(ij-1)*( 1. - gfacy(j-1) );
      
      u(ij) = u(ij) - ( ppe - ppw )*dy*apu(ij);
      v(ij) = v(ij) - ( ppn - pps )*dx*apv(ij);
      p(ij) = p(ij) + urf(ip)*( pp(ij) - ppo );
      
    end
  end  
    
endfunction

%%% end of calcp %%%

%%% solve for temperature %%%

function calct
% solve for temperature

% global variables needed by this function
  global nij nim njm nj
  global urf ien
  global gfacx gfacy
  global xf yf xc yc
  global li
  global visc prr densit dtr
  global f1 f2
  global gds
  global ltime
  global t0

% global variables set or modified by this function
  global su ap awe aso ano aea
  global t
  
%%% initialize su and ap to zero

  
%%% compute convection and diffusion contributions for east and west faces
% set contributions to aea, awe, and su
% use CDS for diffusion
% use CDS/UDS blend for convection
% keep UDS convection on left-hand side of the linear system (implicit terms),
%   and put the rest of the convection contribution on the right-hand side
%   as a deferred correction; see calcuv for examples
  

%%% compute convection and diffusion contributions for north and south faces
% set contributions to ano, aso, and su
% use CDS for diffusion
% use CDS/UDS blend for convection
% keep UDS convection on left-hand side of the linear system (implicit terms),
%   and put the rest of the convection contribution on the right-hand side
%   as a deferred correction; see calcuv for examples
  

%%% compute volume integrals
% here the only contribution is from the unsteady term


%%% apply bcs

  
%%% apply under-relaxation to ap and su
% see calcuv for examples
  

%%% solve for temperature

  
endfunction

%%% end of calct %%%

%%% generate 2d contour plots of dependent variables %%%

function contourplots
% generate 2d contout plots of computed dependent variable fields
% simple version works only for equal number of cells in x and in y

% global variables needed by this function
  global nx ny li
  global xc yc u v p t
  global iu iv ip ien lcal
  
%  global iu iv ip ien lcal

  phi2d = zeros(nx,ny);
  
% global variables set or modified by this function

  if(lcal(iu)==1)
    for i=1:nx
      for j=1:ny
        ij         = li(i) + j;
        phi2d(j,i) = u(ij);
      end
    end
    figure;
    contourf(xc,yc,phi2d);
    colorbar;
    xlabel('x');
    ylabel('y');
    title('x velocity component versus x and y');
    saveas(gcf,'ucon_steady.png');
  endif
  if(lcal(iv)==1)
    for i=1:nx
      for j=1:ny
        ij         = li(i) + j;
        phi2d(j,i) = v(ij);
      end
    end
    figure;
    contourf(xc,yc,phi2d);
    colorbar;
    xlabel('x');
    ylabel('y');
    title('y velocity component versus x and y');
    saveas(gcf,'vcon_steady.png');
  endif
  if(lcal(ip)==1)
    for i=1:nx
      for j=1:ny
        ij         = li(i) + j;
        phi2d(j,i) = p(ij);
      end
    end  
    figure;
    contourf(xc,yc,phi2d);
    colorbar;
    xlabel('x');
    ylabel('y');
    title('Pressure versus x and y');
    saveas(gcf,'pcon_steady.png');
  endif
  if(lcal(ien)==1)
    for i=1:nx
      for j=1:ny
        ij = li(i) + j;
        phi2d(j,i) = t(ij);
      end
    end
    figure;
    contourf(xc,yc,phi2d);
    colorbar;
    xlabel('x');
    ylabel('y');
    title('Temperature versus x and y');
    saveas(gcf,'tcon_steady.png');
  endif  
  
endfunction

%%% end of contourplots %%%

%%%%% execute the program %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate the grid and set related geometric parameters
grid_2d;

% set internal solution control parameters and variables
solnparm;

% set boundary conditions for dependent variables
%   here the BCs are constant in time
bcs;

% set initial conditions for dependent variables
ics;

% open output files
sipfile = fopen('sip_steady.txt','w');
resfile = fopen('res_steady.txt','w');
monfile = fopen('mon_steady.txt','w');

%%% march solution in time from the initial condition for itst time steps %%%
% note that itst=1 for a steady-state problem (with ltime=0)

time = 0;
for itim=1:itst
  time = time + dt;
  
% save values from end of previous time step
  if(ltime==1);
    for ij=1:nij
      u0(ij) = u(ij);
      v0(ij) = v(ij);
      t0(ij) = t(ij);
    end
  endif 

% update for time-dependent BCs
  if(ltime==1);
    bctime;
  endif
  
%%% outer iteration loop %%%
  iter   = 0;
  source = 2.*sormax;
  while ( iter<=maxit && source>sormax );
    iter = iter + 1;
    
    if( lcal(iu)==1 );
      calcuv;
    endif
    if( lcal(ip)==1 );
      calcp;
    endif
    if( lcal(ien)==1 );
      calct;
    endif
    
% check for convergence of outer iterations
    source = max( resor );
    if( source>slarge );
      printf('solution diverging');
      break;
    endif

% write end-of-outer-iteration output
  fprintf(resfile,'%+15e %15i %+15e %+15e %+15e %+15e %+15e %+15e %+15e %+15e \n', ...
                   time,iter,u(ijmon),v(ijmon),p(ijmon),t(ijmon),resor(1),resor(2),resor(3),resor(4));
     
  endwhile % end of outer iteration loop
   
% write end-of-timestep output
    fprintf(resfile,'%+15e %15i %+15e %+15e %+15e %+15e %+15e %+15e %+15e %+15e \n', ...
                   time,iter,u(ijmon),v(ijmon),p(ijmon),t(ijmon),resor(1),resor(2),resor(3),resor(4));
    fprintf(monfile,'%+15e %+15e %+15e %+15e %+15e \n', ...
                     time,u(ijmon),v(ijmon),p(ijmon),t(ijmon));
    
% plot end-of-timestep fields

end % end of time loop

% close output files
fclose(sipfile);
fclose(resfile);
fclose(monfile);

% plot final fields
contourplots;

%%%%% all done %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
