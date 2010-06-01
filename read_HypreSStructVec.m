function [v,nx,ny,nz] = read_HypreSStructVec(vecfile)

% Usage: [v,nx,ny,nz] = read_HypreSStructVec('vecfile')
% 
% reads in the file 'vecfile' and determines the grid structure and
% resulting the data array that determines the associated vector values.
% The vector data is then converted to a standard vector format, in which
% unknowns are ordered with the x-dimension first, then y, then z; e.g. a
% 2x2x2 problem will index the entries as  
%       (1,1,1) -> 1
%       (2,1,1) -> 2
%       (1,2,1) -> 3
%       (2,2,1) -> 4
%       (1,1,2) -> 5
%       (2,1,2) -> 6
%       (1,2,2) -> 7
%       (2,2,2) -> 8
%
% Outputs:
%       v - vector containing the specified HYPRE SStruct vector
%       nx - x-dimensional size of grid
%       ny - y-dimensional size of grid
%       nz - z-dimensional size of grid
%
% Daniel R. Reynolds
% drreynolds@ucsd.edu

% open file
fid = fopen(vecfile,'r');

% get the grid dimensions
line = fgets(fid);  % should be 'StructVector'
line = fgets(fid);  % should be blank
line = fgets(fid);  % should be 'Grid:'
ndim = fscanf(fid,'%i',1);  
line = fgets(fid);  % rest of line
if ((ndim < 1) | (ndim > 3))
   disp(sprintf('error: illegal ndim = %i',ndim));
   return
end
line = fgets(fid);  % should contain the number 1 (don't know what it's for)
extents = fscanf(fid,'0:  (%i, %i, %i) x (%i, %i, %i)',6);  
line = fgets(fid);  % rest of line
nx = extents(4)-extents(1)+1;
ny = extents(5)-extents(2)+1;
nz = extents(6)-extents(3)+1;
if (ndim < 3)
   if (nz ~= 1)
      disp(sprintf('error, ndim = %i and nz = %i do not match',ndim,nz))
      return
   end
end
if (ndim == 1)
   if (ny ~= 1)
      disp(sprintf('error, ndim = 1 and ny = %i do not match',ny))
      return
   end
end


% read the data and put into vector
v = zeros(nx*ny*nz,1);
line = fgets(fid);  % should be blank
line = fgets(fid);  % should be 'Data:'
for i=1:nx*ny*nz
   vecrow = fscanf(fid,'%i: (%i, %i, %i; %i) %g',6);
   line = fgets(fid);  % rest of line
   ix = vecrow(2)-extents(1)+1;  % shift to start at 1
   iy = vecrow(3)-extents(2)+1;  % shift to start at 1
   iz = vecrow(4)-extents(3)+1;  % shift to start at 1
   
   irow = ((iz-1)*ny + (iy-1))*nx + ix;
   v(irow) = vecrow(6);
end


% close hypre matrix file
fclose(fid);


% end of routine