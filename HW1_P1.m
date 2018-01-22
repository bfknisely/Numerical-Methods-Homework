%{
Brian Knisely
ME523, HW1, P1

The purpose of this code is to determine the regions in space in which a 
PDE is characterized as elliptic, parabolic, or hyperbolic

Given PDE is
u_xx + x*u_xy + y*u_yy = 0
Typical format is
a*u_xx + b*u_xy + c*u_yy = f

Character determined by value of b^2 - 4*a*c
%}

clear; close all; format compact; home;

x = -10:0.01:10;  % Range of x
y = -10:0.01:10;  % Range of y
[b, c] = meshgrid(x, y);
% make mesh grid for x, y locations and define coefficients b and c

a = 1;  % Set first coefficient equal to 1

ch = b.^2 - 4.*a.*c;  % Compute value of character at every x-y location

v = [-100,0,200];  % Set contour levels (so the plot has distinct regions)
contourf(b, c, ch, v, 'linewidth',6,'linecolor','r'); colormap gray
% Plot filled contour in b, c (x, y) space with character as z-values
xlabel('x'); ylabel('y');  % Label axes

text(0, 5, 'Elliptic', 'horizontalAlignment', 'center', 'color', 'w');
text(0, -0.5, 'Parabolic (along curve)', ...
    'horizontalAlignment', 'center', 'color', 'r');
text(0, -6, 'Hyperbolic', 'horizontalAlignment', 'center', 'color', 'k');
% Add text to plot to show regions