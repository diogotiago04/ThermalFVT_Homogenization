
function ThermalFVT_energy_based (nx, ny, frac, k_m, k_i)

R = sqrt(frac*nx*ny/pi);                     % radius of circular material heterogeneity 

% Subvolume faces
[i,j] = meshgrid(1:nx,1:ny);
s = (i+(j-1)*nx)';
faces = [s(:),s(:)+nx*ny+1,s(:)+nx,s(:)+nx*ny];
faces(end-nx+1:end,3) = faces(1:nx,1);
faces(nx:nx:end,2) = faces(1:nx:end-nx+1,4);

% Degrees of freedom
ndof = max(faces(:));
dofIS = unique([faces(1:nx, 1); faces(end-nx+1:end, 3)]);
dofIS = [dofIS(1),dofIS(end)];
dofDE = unique([faces(nx:nx:end, 2); faces(1:nx:end-nx+1, 4)]);
dofDE = [dofDE(1),dofDE(end)];

% Degrees of freedom: fixed and free
fixed = [dofIS dofDE];
free = setdiff(1:ndof,fixed);

% Sparse mapping indices
iK = reshape(kron(faces,ones(4,1))',16*nx*ny,1);
jK = reshape(kron(faces,ones(1,4))',16*nx*ny,1);
iF = repmat((faces)',2,1);                
jF = [ones(4,nx*ny); 2*ones(4,nx*ny)];


%____________________________________LOCAL CONDUCTIVITY MATRIX AND HEAT FLUX VECTOR

k0 = eye(2);

% Auxiliary matrices
a = ones(4,1);
N1 = [0,-1]; N2 = [1,0]; N3 = [0,1]; N4 = [-1,0];
N = [N1,zeros(1,6); zeros(1,2),N2,zeros(1,4);
    zeros(1,4),N3,zeros(1,2); zeros(1,6),N4];
A = [0 -1/2 0 1/4; 1/2 0 1/4 0;
    0 1/2 0 1/4; -1/2 0 1/4 0];
E = [0 0 0 0; 0 -1 0 3/2; -1 0 -3/2 0; 0 0 0 0;
    0 0 0 0; 0 -1 0 -3/2; -1 0 3/2 0; 0 0 0 0];
B = N*eye(8)*E;
ab = (B*(A\a))\(B/A);
Ab = A\(eye(4)-a*ab);

% Local conductivy matrix
K0 = B*Ab;

% Local heat flux vector
H0 = [N1;N2;N3;N4]*k0;

x = InitialMaterialDesign(nx,ny,R,k_m,k_i);

%_____________________________________________________PREALLOCATE VARIABLES

Tf = zeros(ndof,2);        % surface-averaged fluctuating temperatures
T = zeros(nx*ny,4,2);      % total surface-averaged temperatures (macroscopic + fluctuating)
C = zeros(2);              % effective thermal conductivity matrix

% Compute macroscopic temperatures
T0 = cell(nx*ny,1);
for j = 1:ny
    for i = 1:nx, s = i+(j-1)*nx;     
T0{s} = [ (1/2+(i-1)), (j-1);...
          (1  +(i-1)), (j-1/2);...
          (1/2+(i-1)),  j;...
          (i-1), (j-1/2) ];
    end
end

% FINITE-VOLUME THEORY ANALYSIS

% Material Interpolation
sK = K0(:)*x(:)';
sF = H0(:)*x(:)';

% Assembly global conductivy matrix
K = sparse(iK, jK, sK, ndof, ndof); K = (K+K')/2;

% Assembly of heat flux vectors corresponding to two unit temperature gradient tests
Q0 = sparse(iF(:), jF(:), sF, ndof, 2);

% Compute fluctuating temperatures for two unit temperature gradient tests
Tf(free,:) = K(free,free) \ Q0(free,:);

% MATERIALS' HOMOGENIZATION

% Compute total temperatures
for s = 1:(nx*ny)
            T(s,:,:) = T0{s}+Tf(faces(s,:),:);
        end
        for i = 1:2, Ti = T(:,:,i);
            for j = 1:2, Tj = T(:,:,j);
                sumE = reshape(sum((Ti*K0).*Tj,2),nx,ny)/(nx*ny);
                C(i,j) = sum(sum(x.*sumE));     
            end
        end
        disp(-C);
end

%___________________________________________________INITIAL MATERIAL DESIGN
function x = InitialMaterialDesign(nx,ny,R,k_m,k_i)
x = k_m*ones(nx,ny);
for j = 1:ny
    for i = 1:nx
        if sqrt((i-nx/2-0.5)^2 + (j- ny/2-0.5)^2) < R
            x(i,j) = k_i;
        end
    end
end

end

