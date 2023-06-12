clc, clear, close all; 

%% NAMES OF GROUP MEMBERS
% Name 1 - Piper Morris
% Name 2 - Jay Matter

%% ABOUT
% 2D Static Solid Mechanics FEA Code

% Q4 Isoparametric Elements
%   - quadrilateral elements
%   - linear interpolation functions
%   - linear shape functions
%   - quadratic gauss quadrature

%% INPUTS

% SECTION 1:
Lx = 2; % length of geometry in the x-direction [m]
Ly = 2; % length of geometry in the y-direction [m]
Nx = 4; % number of elements in the x-direction
Ny = 4; % number of elements in the y-direction

% SECTION 2:
E_mp = 1E9; % Young's modulus (material property) [Pa]
v_mp = 0.3; % Poisson's ratio (material property)
h_ps = 0.01; % 2D plate height/thickness [m]

% SECTION 3:
DispDir = 'y'; % direction of displacement given by 'x' or 'y'
DispVal = 0.1; % value of displacement [m]

%% OUTPUTS

% K_e - 8x8 element stiffness matrices
% K_G - global (unmodified) stiffness matrix

%% SECTION 1 - Define the elements and nodes

% Calculate the number of nodes in each direction:

nx = Nx+1; %Number of nodes in the x-direction
ny = Ny+1;%Number of nodes in the y-direction

% Calculate the total number of elements and total number of nodes:

m = Nx*Ny; %Total number of elements
n = nx*ny;%Total number of nodes

% Generate the node number matrix:

node_post = reshape( 1:n, ny, nx);

% Define the global x and y coordinates for the nodes (in their undeformed configuration):

A = [0:1:Ny]; %Initialize vectors based on number of elements
B = [0:1:Nx];

for i = 1:nx; %Populate matrices with x and y values of the nodes
    for j = 1:ny;
        x0(j,:) = B;
        y0(:,i) = A;
    end
end

x = (Lx/Nx)*x0; %Scale coordinate values based on inputs
y = (Ly/Ny)*y0;

% Generate element information - **[el_connect, el_x, el_y]**:

row1 = 1; %Initialize values for loops below
row2 = 2;
column1 = 1;
column2 = 2;
i = 0;
%Create output matrices with rows representing values from each element in
%proper order: node_post -> el_connect , x -> el_x , y -> el_y
for j = 1:m
    el_connect(j,1:4) = reshape(node_post(row1:row2,column1:column2)',[1,4]);
    el_x(j,1:4) = reshape(x(row1:row2,column1:column2)',[1,4]);
    el_y(j,1:4) = reshape(y(row1:row2,column1:column2)',[1,4]);
    row1 = row1 + 1;
    row2 = row2 + 1;
    i = i +1;
    if (i==Ny)
        column1 = column1 + 1;
        column2 = column2 + 1;
        row1 = 1;
        row2 = 2;
        i = 0;
    end
end

%% SECTION 2 - Calculate the element and assemble global stiffness matrices

% Calculate the properties c11, c12, c22, and c66 for a PLANE
% STRESS ISOTROPIC material

matrix = (E_mp / (1-(v_mp^2))) * [1 v_mp 0;
    v_mp 1 0;
    0 0 (1-v_mp)/2;];
%Use the given values to calculate the weights
C_11 = matrix(1,1);
C_12 = matrix(2,1);
C_21 = matrix(1,2);
C_22 = matrix(2,2);
C_66 = matrix(3,3);

%Initialize each portion of the ke matrix
for j = 1:m
k_11 = zeros(4,4);
k_12 = zeros(4,4);
k_21 = zeros(4,4);
k_22 = zeros(4,4);

for i = 1:4
%Define xi and eta
xi_1 = -1/sqrt(3); eta_1 = 1/sqrt(3);
xi_2 = 1/sqrt(3); eta_2 = 1/sqrt(3);
xi_3 = -1/sqrt(3); eta_3 = -1/sqrt(3);
xi_4 = 1/sqrt(3); eta_4 = -1/sqrt(3);

%Assemble xi and eta values in vectors
xi = [xi_1 xi_2 xi_3 xi_4];
eta = [eta_1 eta_2 eta_3 eta_4];

% Define the derivatives of the interpolation/shape functions
dN_dxi_1 = (1/4)*(-1 + eta(i)); dN_deta_1 = (1/4)*(-1 + xi(i));
dN_dxi_2 = (1/4)*(1 -eta(i));  dN_deta_2 = (1/4)*(-1 - xi(i));
dN_dxi_3 = (1/4)*(-1 -eta(i));  dN_deta_3 = (1/4)*(1 -xi(i));
dN_dxi_4 = (1/4)*(1 + eta(i)); dN_deta_4 = (1/4)*(1 + xi(i));

%Assemble the derivatives in vectors
dN_dxi = [dN_dxi_1;  dN_dxi_2; dN_dxi_3; dN_dxi_4;];
dN_deta = [dN_deta_1; dN_deta_2; dN_deta_3; dN_deta_4;];

%Extract x and y values from Section 1
x1 = el_x(i,1:4)';
x2 = el_y(i,1:4)';

%Define the components of the Jacobian
    J_11 = dN_dxi'  * x1;
    J_12 = dN_deta' * x1;
    J_21 = dN_dxi'  * x2;
    J_22 = dN_deta' * x2;

%Assemble above values into a 2 x 2 matrix
    J = [J_11 J_12;
       J_21 J_22;];

%Compute the inverse of the Jacobian
    J_inv = inv([J]);

%Define the terms of the approximation function for each ke matrix
    term_one = (dN_dxi * J_inv(1,1)) + (dN_deta * J_inv(1,2));
    term_two = (dN_dxi' * J_inv(1,1)) + (dN_deta' * J_inv(1,2));

    term_three = (dN_dxi * J_inv(2,1)) + (dN_deta * J_inv(2,2));
    term_four = (dN_dxi' * J_inv(2,1)) + (dN_deta' * J_inv(2,2));

%Compute each component of the ke matrix
New_k_11 = h_ps * det([J]) *( ( C_11 * term_one * term_two) + (C_66 * term_three * term_four) );
k_11 = k_11 + New_k_11;

New_k_12 = h_ps * det([J]) *( ( C_12 * term_one * term_four) + (C_66 * term_three * term_two) );
k_12 = k_12 + New_k_12;

New_k_21 = h_ps * det([J]) *( ( C_21 * term_three * term_two) + (C_66 * term_one * term_four) );
k_21 = k_21 + New_k_21;

New_k_22 = h_ps * det([J]) *( ( C_22 * term_three * term_four) + (C_66 * term_one * term_two) );
k_22 = k_22 + New_k_22;

%Assemble the ke matrix
ke = [k_11 k_12;
    k_21 k_22;];
        
end 

end

%Initialize the global stiffness matrix
k_G = zeros(n*2,n*2);

%Create a vector that defines the indices of k_G
New_el_connect = el_connect + n;
New_new_el_connect = [el_connect New_el_connect];

%Extract the values of the ke to populate the k_G
for h = 1:m
index = New_new_el_connect(h,:);
New_KG = zeros(n*2,n*2);

for k = 1:8 

matrix = ke(k,:);
New_KG(index(k),index) = matrix;

end 

k_G = k_G + New_KG;

end

%% SECTION 3 - Calculate the displacements and forces at each of the nodes

% Modify Global Stiffness Matrix [K_Gmod] Based on Fixed Boundary Conditions

k_Gmod = k_G;
%creat a matrix for the y components of each node to use in U
new_nodes = node_post + n ;


nodes = node_post(:,1); % The nodes corresponding to the x disp of the left side in global matrix
nodes_2 = new_nodes(:,1); % The nodes corresponding to the y disp of the left side in global matrix 

% create a vector with the node values of the leftmost column 
index_2 = [nodes; nodes_2;]';

% Modify Global Stiffness Matrix Based on Displaced Boundary Conditions
for i = 1:ny*2

k_Gmod(index_2(i),:) = 0;
k_Gmod(index_2(i),index_2(i)) = 1;

end

%Initialize the modified Force vector
F_mod = zeros(n*2,1);

%Extract last column of node location matrices
x_nodes = node_post(:,nx)';  
y_nodes = new_nodes(:,nx)';

for l = 1:ny % number of nodes in the x 
% Displace the right side in either the x or the y
if DispDir ==  'x'
    
      k_Gmod(x_nodes(l),:) = 0;
      k_Gmod(x_nodes(l),x_nodes(l)) = 1;

      F_mod(x_nodes,1) = DispVal;

elseif DispDir == 'y'
        
       k_Gmod(y_nodes(l),:) = 0;
       k_Gmod(y_nodes(l),y_nodes(l)) = 1;

       F_mod(y_nodes,1) = DispVal

end 

end

%Display the modified global stiffness matrix
k_Gmod

% Calculate the displacement at each of the nodes using the equation: U = K_Gmod\F_mod
U = k_Gmod\F_mod

% Calculate the force at each of the nodes using the equation: F = K_G*U
F = k_G*U

% Generate a plot of the undeformed and deformed configurations
x2 = x(:);
y2 = y(:);
undef = [x2;y2];
plot(x2,y2,'.','markersize',30,'color','b')
hold on
grid on
def = undef + U;
x3 = def(1:n);
y3 = def(n+1:2*n);
plot(x3,y3,'.','markersize',30,'color','r')
legend('undeformed','deformed','location','northoutside','orientation','horizontal');
