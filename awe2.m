%%AWE Weak form implementation
%%SUPG stabalisation 

clc
clear all

%%%%% Definition of Discretized Domain %%%%%

%speed of sound 
c = 100;

% Length of the domain
Lx = 1 ;
Ly = 1 ;
% Number of elements
nElemX = 5;
nElemY = 5;
nElem = nElemX*nElemY ;
nCrds = (nElemX+1)*(nElemY+1) ;
% Number of nodes per element
nne = 4 ;
% Number of elements in 1D 
nne1 = 2;
% Source term
f = 0 ;
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
a1 = 0 ;
for i=1:nElem
    conn(i,1) = i+a0 ;
    conn(i,2) = conn(i,1)+1 ;
    conn(i,3) = i+nElemX+2+a1 ;
    conn(i,4) = conn(i,3)-1 ;
    if (mod(i,nElemX)==0)
        a0 = a0+1 ;
        a1 = a1+1 ;
    end
end

%%%%% Definition of Gauss Quadrature and Shape Function Space %%%%%

% Location of Gauss points
gP = [-1/sqrt(3),  1/sqrt(3), 1/sqrt(3), -1/sqrt(3);
      -1/sqrt(3), -1/sqrt(3), 1/sqrt(3),  1/sqrt(3)];
% Weights of Gauss points
gW = [1,  1,  1,  1] ;
% Number of Gauss points
nQuad = length(gW) ;

% 2DShape functions
N(1,:) = 0.25.*(1-gP(1,:)).*(1-gP(2,:)) ;
N(2,:) = 0.25.*(1+gP(1,:)).*(1-gP(2,:)) ;
N(3,:) = 0.25.*(1+gP(1,:)).*(1+gP(2,:)) ;
N(4,:) = 0.25.*(1-gP(1,:)).*(1+gP(2,:)) ;

% Gradient of shape functions
Nx(1,:) = -0.25.*(1-gP(2,:)) ;
Nx(2,:) =  0.25.*(1-gP(2,:)) ;
Nx(3,:) =  0.25.*(1+gP(2,:)) ;
Nx(4,:) = -0.25.*(1+gP(2,:)) ;
Ny(1,:) = -0.25.*(1-gP(1,:)) ;
Ny(2,:) = -0.25.*(1+gP(1,:)) ;
Ny(3,:) =  0.25.*(1+gP(1,:)) ;
Ny(4,:) =  0.25.*(1-gP(1,:)) ; 

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

%%%%Gaussian intial function 
ndof = size(crd,1) ;
[X, Y] = meshgrid(1:nElemX + 1, 1: nElemY + 1);
%G = (1/(pi*2*20)).*exp(-1/(20^2).*((Y-(nElemY+1)/2).^2 + (X-(nElemX+1)/2).^2));
G = (X-1).*(Y-1).*(X-6).*(Y-6);
G  = reshape(G,ndof,1);
%G(BCNodes) = 0;
G12 = reshape(G,nElemX+1,nElemY+1);

%Initial condition 
Sol.u = crd(:,1).*crd(:,2).*(Lx - crd(:,1)).*(Ly - crd(:,2));

%%%%% Initialization of Solution Vector %%%%%

Sol.ddu = zeros(ndof,1);
Sol.dxu = zeros(ndof,1);
Sol.du = zeros(ndof,1);
Sol.iu =zeros(ndof,1);

 figure(1);
        mesh(X',Y',reshape(Sol.u,(nElemX+1),(nElemY+1)));
        set(gca,'TickLabelInterpreter','latex','FontSize',30);
        g1 = colorbar;
        colormap('jet');
        set(g1,'TickLabelInterpreter','latex','FontSize',30);
        xlabel('$x$','Interpreter','latex');
        ylabel('$y$','Interpreter','latex');
        zlabel('$Pa$','Interpreter','latex');
        title("Initial conditions ");


% Satisfy Dirichlet boundary conditions
Sol.u(BCNodes) = BCValues ;
Sol.iu(BCNodes) = BCValues ;

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
            Kij_1g = -gW(p)*(DNDx(:,i).*DNDx(:,j));
            Kij_2g = -gW(p)*(DNDy(:,i).*DNDy(:,j));
            Kij_1g = Kij_1g.*absJac;
            Kij_2g = Kij_2g.*absJac;
            sK_1g(index+1:index+nElem,p) = Kij_1g;
            sK_2g(index+1:index+nElem,p) = Kij_2g;

            %Mass matrix terms 
             Mijg = gW(p)*(N(p,j)*N(p,i));
             Mijg = (1/(c*c))*Mijg.*absJac ;
             sMg(index+1:index+nElem,p) = Mijg;
                    
            % Galerkin source term
            Fijg = gW(p)*(N(p,i)*N(p,j));
            Fijg = Fijg.*absJac ;
            sFg(index+1:index+nElem,p) = Fijg;

            index = index + nElem;
        end
    end
end


% Summation of All Quadrature Data for Numerical Integration
sK_1g = sum(sK_1g,2);
sK_2g = sum(sK_2g,2);
sFg = sum(sFg,2);
sMg = sum(sMg,2);

%%%%% Assembly of Local Matrices to Global Matrix %%%%%
K1_g = sparse(ii,jj,sK_1g,ndof,ndof); 
K2_g = sparse(ii,jj,sK_2g,ndof,ndof); 
F_g = sparse(ii,jj,sFg,ndof,ndof); 
M_g = sparse(ii,jj,sMg,ndof,ndof);

%setting unknowns 
unKnowns = setdiff([1:ndof]',BCNodes) ;

%%Time Discritisation 

Sol.ndu = Sol.du ; 
Sol.nu = Sol.u;
Sol.nddu = Sol.ddu;

%Taking inputs of params 
rho = 1;
alpham = 0.5*((3 - rho)/(1+rho));
alpha = 1/(1+rho);
gamma = 0.5 + alpham - alpha ;
beta = 0.25*(1-alpham + alpha);

%time step input 
dt = 01;
%Number of time steps 
nt = 200;
%starting time iterations 
for ni = 1:nt
           %%%%% Defining of K Matrix and RHS Vector %%%%%
           K1 = K1_g + K2_g ;
           force = f.*ones(ndof,1);
           RHS1 = (F_g)*force ;

           %%%%% Incorporating Boundary Conditions to the RHS %%%%%
          RHS1 = RHS1 - K1*Sol.u;

          %%%final matrices 
         M = M_g ;
         K=K1;
         %damping terms 
         a = 10;
         b = 10;
         C = a* M +  b*K;
         F=RHS1;

        %LHS
        LHS = ((2)./(dt*dt)).*M + K./2 + C./dt;

        %RHS
        RHS = F + (((2)./(dt*dt)).*M - K./2  + C./dt)*Sol.u + ((2)./(dt)).*M*Sol.du;

        %solving the solution  
        Sol.nu(unKnowns) = LHS(unKnowns,unKnowns)\RHS(unKnowns) ;

        %updating velocity 
        Sol.ndu = (2/dt).*(Sol.nu - Sol.u) - Sol.du ; 

        %Updating Variables 
        Sol.du = Sol.ndu ; 
        Sol.u = Sol.nu;
        Sol.ddu = Sol.nddu;

        %%%%% Post-processing %%%%%
        figure(3);
        mesh(X',Y',reshape(Sol.u,(nElemX+1),(nElemY+1)));
        set(gca,'TickLabelInterpreter','latex','FontSize',30);
        g1 = colorbar;
        colormap('jet');
        set(g1,'TickLabelInterpreter','latex','FontSize',30);
        xlabel('$x$','Interpreter','latex');
        ylabel('$y$','Interpreter','latex');
        zlabel('$Pa$','Interpreter','latex');
        title("Time marching ");
end
