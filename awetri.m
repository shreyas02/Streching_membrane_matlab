%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    2D Poisson equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all

%%%%% Definition of Discretized Domain %%%%%

% Length of the domain
Lx = 1 ;
Ly = 1 ;
% Number of elements
nElemX = 6 ;
nElemY = 6 ;
nElem = 2*nElemX*nElemY ;
nCrds = (nElemX+1)*(nElemY+1) ;
% Number of nodes per element
nne = 3 ;
% Source term
f = 0 ;
%Speed of sound 
c = 0.1;
% Boundary nodes
BCBottom = [1:1:nElemX+1]';
BCLeft = [1:nElemX+1:nCrds]';
BCRight = [nElemX+1:nElemX+1:nCrds]';
BCTop = [nCrds-nElemX:1:nCrds]';
BCNodes = unique([BCBottom; BCLeft; BCRight; BCTop]);
BCValues = zeros(size(BCNodes,1),1);
for i=1:size(BCTop,1)
    bb = ismember(BCNodes, BCTop(i,1));
    BCValues(bb~=0) = 0 ;
end

% Coordinates of the nodes
x = linspace(0,Lx,nElemX+1);
y = linspace(0,Ly,nElemY+1);
[X,Y] = meshgrid(x,y);
crd = [reshape(X',nCrds,1) reshape(Y',nCrds,1)];

% Connectivity matrix of the elements
conn = zeros(nElem,nne); 
nn = 1;
a0 = 0 ;
a1 = 1 ;
for i=1:nElem
      if (mod(i,2) ~=0)
            conn(i,1) = i-a0 ;
            conn(i,2) = conn(i,1)+1 ;
            conn(i,3) = i-a0+nElemX+1 ;
            if (mod(i+1,2*nElemX)~=0)
                a0 = a0+1 ;
            end
        elseif (mod(i,2) ==0)
            conn(i,1) = i-a1+nElemX+1 ;
            conn(i,2) = conn(i-1,2) ;
            conn(i,3) = conn(i,1)+1 ;
            if (mod(i,2*nElemX)~=0)
                a1 = a1+1 ;
            end
      end
end

%%%%% Definition of Gauss Quadrature and Shape Function Space %%%%%

% Location of Gauss points
gP = [1/3, 4/3, 1/3;
      1/3, 1/3, 4/3];
% Weights of Gauss points
gW = [2/3,  2/3,  2/3] ;
% Number of Gauss points
nQuad = length(gW) ;
% Shape functions
N(1,:) = 0.5.*(2-gP(1,:)-gP(2,:)) ;
N(2,:) = 0.5.*(gP(1,:)) ;
N(3,:) = 0.5.*(gP(2,:)) ;
% Gradient of shape functions
Nx(1,:) = -0.5.*ones(1,3) ;
Nx(2,:) =  0.5.*ones(1,3) ;
Nx(3,:) =  0.*ones(1,3) ;
Ny(1,:) = -0.5.*ones(1,3) ;
Ny(2,:) =  0.*ones(1,3) ;
Ny(3,:) =  0.5.*ones(1,3) ;

%%%%% Formation of Local to Global DOF Mapping %%%%%

ii = zeros(nne^2*nElem,1); 
jj = zeros(nne^2*nElem,1);
index = 0;
for i = 1:nne
   for j = 1:nne
      ii(index+1:index+nElem) = double(conn(:,i)); 
      jj(index+1:index+nElem) = double(conn(:,j));  
      index = index + nElem;
   end
end

%%%%% Initialization of Solution Vector %%%%%
ndof = size(crd,1) ;
Sol.u = zeros(ndof,1);
Sol.iu = zeros(ndof,1);
Sol.du = zeros(ndof,1);
Sol.ddu = zeros(ndof,1);
Sol.nu = zeros(ndof,1);
Sol.ndu = zeros(ndof,1);
Sol.nddu = zeros(ndof,1);

%Initial condition 
Sol.u = crd(:,1).*crd(:,2).*(Lx - crd(:,1)).*(Ly - crd(:,2));

% Satisfy boundary conditions
Sol.u(BCNodes) = BCValues ;
Sol.iu(BCNodes) = BCValues ;
Sol.nu(BCNodes) = BCValues ;

%%%%% Post-processing %%%%%
figure(1);
trimesh(conn,X',Y',Sol.u);
set(gca,'TickLabelInterpreter','latex','FontSize',30);
g1 = colorbar;
colormap('jet');
set(g1,'TickLabelInterpreter','latex','FontSize',30);
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$u$','Interpreter','latex');
title("Initial condition");

%%%%% Initialization for Element-level Calculations %%%%%

xx = zeros(size(conn));
yy = zeros(size(conn));

for i=1:nne
   xx(:,i) = crd(conn(:,i),1);
   yy(:,i) = crd(conn(:,i),2);
end

%%%%% Element-level Evaluation of Matrices %%%%%

% Gauss quadrature loop
for p = 1:nQuad  
    
    % Jacobian evaluation and its absolute value
    J = [xx*[Nx(:,p)], yy*[Nx(:,p)],...
         xx*[Ny(:,p)], yy*[Ny(:,p)]];
    absJac =abs( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
    
    % Evaluation of gradients of shape functions 
    DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,2))*Ny(:,p)')./repmat(absJac,1,nne);
    DNDy = ((-J(:,3))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(absJac,1,nne);
            
    index = 0;
    for i = 1:nne
        for j = 1:nne
            
            % Galerkin diffusion term
            Kij_1 = gW(p)*(DNDx(:,i).*DNDx(:,j));
            Kij_2 = gW(p)*(DNDy(:,i).*DNDy(:,j));
            Kij_1 = -1 * Kij_1.*absJac;
            Kij_2 = -1 * Kij_2.*absJac;
            sK_1(index+1:index+nElem,p) = Kij_1;
            sK_2(index+1:index+nElem,p) = Kij_2;

            %Mass matrix term
            Mij = gW(p)*(N(p,i)*N(p,j));
            Mij = (1/(c*c)).*Mij.*absJac;
            sM(index+1:index+nElem,p) = Mij;
                    
            % Galerkin source term
            Fij = gW(p)*(N(p,i)*N(p,j));
            Fij = Fij.*absJac ;
            sF(index+1:index+nElem,p) = Fij;

            index = index + nElem;
        end
    end
end

% Summation of All Quadrature Data for Numerical Integration
sK_1 = sum(sK_1,2);
sK_2 = sum(sK_2,2);
sM = sum(sM,2);
sF = sum(sF,2);

%%%%% Assembly of Local Matrices to Global Matrix %%%%%
K1_global = sparse(ii,jj,sK_1,ndof,ndof); 
K2_global = sparse(ii,jj,sK_2,ndof,ndof); 
M_global = sparse(ii,jj,sM,ndof,ndof);
F_global = sparse(ii,jj,sF,ndof,ndof); 

%%%%% Defining Matrices %%%%%
K = K1_global + K2_global ;
M = M_global;
Force = (f.*ones(ndof,1));
G = F_global*Force ;

%%%%% Time step integration  %%%%%

% Selecting unknown degrees of freedom for solution
unKnowns = setdiff([1:ndof]',BCNodes) ;

%Setting parameters for generalised alpha method 
rho = 1;
alpham = 0.5*((3 - rho)/(1+rho));
alpha = 1/(1+rho);
gamma = 0.5 + alpham - alpha ;
beta = 0.25*(1-alpham + alpha);
dt = 01; 
nt =3;

%Time integration loop 
for i = 1:nt
        i
        LHS = ((alpham)/ (dt*dt*beta)).*M + (alpha * 0.5).*K ; 
        RHS = G + ((alpham)/ (dt*beta)).*M*Sol.du + ((alpham - 2*beta)/(2*beta)).*M*Sol.ddu + (((alpham)/ (dt*dt*beta)).*M - (1-alpha).*K)*Sol.u;

        %%%%% Incorporating Boundary Conditions to the RHS %%%%%
        RHS = RHS - LHS*Sol.iu ;

        % Solve for the solution
        Sol.nu(unKnowns) = LHS(unKnowns,unKnowns)\RHS(unKnowns) ;

        %update acceleration 
        Sol.nddu = ((Sol.nu - (Sol.u + dt.*Sol.du + (0.5 - beta)*dt*dt.*Sol.ddu))./(beta*dt*dt));
        
        %update velocity 
        Sol.ndu = Sol.du + dt.*((1-gamma).*Sol.ddu - gamma.*Sol.nddu);

        %resetting variables 
        Sol.u = Sol.nu;
        Sol.du = Sol.ndu;
        Sol.ddu = Sol.nddu;

        %%%%% Post-processing %%%%%
        figure(2);
        trimesh(conn,X',Y',Sol.u);
        set(gca,'TickLabelInterpreter','latex','FontSize',30);
        g1 = colorbar;
        colormap('jet');
        set(g1,'TickLabelInterpreter','latex','FontSize',30);
        xlabel('$x$','Interpreter','latex');
        ylabel('$y$','Interpreter','latex');
        zlabel('$u$','Interpreter','latex');
        title("Time marching " );
end


