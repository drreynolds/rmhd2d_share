function [rho,u,v,w,Bx,By,Bz,p,divB,J,e,x,y,z] = read_output2d(fname,Nx,Ny)
%
% [rho,u,v,w,Bx,By,Bz,p,divB,J,e,x,y] = read_output2d(fname,Nx,Ny)
%
% Extracts 2D fields from gnuplot output file.
%
% Inputs:  
%     fname - file name for input
%        Nx - dimension of fields in x-direction
%        Ny - dimension of fields in y-direction
%
% Outputs:  
%       rho - density field (2D)
%         u - x-velocity field (2D)
%         v - y-velocity field (2D)
%         w - z-velocity field (2D)
%        Bx - x-magnetic field (2D)
%        By - y-magnetic field (2D)
%        Bz - z-magnetic field (2D)
%         p - pressure field (2D)
%      divB - div(B) field [for constraint check] (2D)
%         J - toroidal current (2D)
%         e - total energy field (2D)
%         x - x-coordinates (1D)
%         y - y-coordinates (1D)
%
% Daniel Reynolds
% Mathematics @ SMU
%

% open file for input
fin = fopen(fname,'r');
if fin < 0
   error([ 'Could not open ',fname,' for input']);
end

% check whether there is a label line (marked by #)
nhead = 0;
for i=1:100
   buffer = fgetl(fin);
   [text,buffer] = strtok(buffer);
   if (text == '#') 
      nhead = nhead + 1;
   else
      break;
   end
end
frewind(fin);
   
% discard header text in preparation of reading data
for i=1:nhead,  buffer=fgetl(fin);  end

% initialize output arrays
rho = zeros(Nx,Ny);
u = zeros(Nx,Ny);
v = zeros(Nx,Ny);
w = zeros(Nx,Ny);
Bx = zeros(Nx,Ny);
By = zeros(Nx,Ny);
Bz = zeros(Nx,Ny);
p = zeros(Nx,Ny);
divB = zeros(Nx,Ny);
J = zeros(Nx,Ny);
e = zeros(Nx,Ny);
x = zeros(Nx,1);
y = zeros(Ny,1);

% read in data
for j=1:Ny
   for i=1:Nx

      % read in line of text
      buffer = fgetl(fin);
      
      % extract fields into appropriate output arrays
      [x(i), y(j), rho(i,j), u(i,j), v(i,j), w(i,j), Bx(i,j), ...
	     By(i,j), Bz(i,j), p(i,j), divB(i,j), J(i,j), e(i,j)] ...
	  = strread(buffer,'%f %f %f %f %f %f %f %f %f %f %f %f %f');
   
   end
   
   % skip a line
   buffer = fgetl(fin);
   
end


% close file
fclose(fin);
