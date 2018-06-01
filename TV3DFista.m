function [X_den,iter,fun_all]=TV3DFista(Xobs,lambda,l,u,pars)
%This function implements the FISTA method for TV-based denoising problems
%
% Based on the paper
% Amir Beck and Marc Teboulle, "Fast Gradient-Based Algorithms for Constrained
% Total Variation Image Denoising and Deblurring Problems"
% -----------------------------------------------------------------------
% Copyright (2008): Amir Beck and Marc Teboulle
% Modified by Nick Antipa, UC Berkeley, 2017
%
% FISTA is distributed under the terms of 
% the GNU General Public License 2.0.
% 
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ----------------------------------------------------------------------
% INPUT
% Xobs ..............................an observed noisy image.
% lambda ........................ parameter
% l ..................................... lower bound on the pixels' values
% u ..................................... upper bound on the pixels' values
% pars.................................parameters structure
% pars.MAXITER ..................... maximum number of iterations
%                                                      (Default=100)
% pars.epsilon ..................... tolerance for relative error used in
%                                                       the stopping criteria (Default=1e-4)
% 
% OUTPUT
% X_den ........................... The solution of the problem 
%                                            min{||X-Xobs||^2+2*lambda*TV(X
%                                            ) : l <= X_{i,j} <=u} 
% iter .............................  Number of iterations required to get
%                                            an optimal solution (up to a tolerance)
% fun_all ......................   An array containing all the function
%                                             values obtained during the
%                                             iterations


%Define the Projection onto the box

% Assigning parameres according to pars and/or default values
flag=exist('pars','var');
if (flag&&isfield(pars,'MAXITER'))
    MAXITER=pars.MAXITER;
else
    MAXITER=100;
end
if (flag&&isfield(pars,'epsilon'))
    epsilon=pars.epsilon;
else
    epsilon=1e-4;
end

project = @(x)min(max(x,l),u);    %Function that enforces bounding value

% Forward operator for dual problem. This is basically just finite
% difference each dimension. It is vector valued, so each pixel creates
% needs to output n values for an n-dimensional TV denoising problem.
Ltv3 = @(P1,P2,P3)cat(1,P1(1,:,:),diff(P1,1,1),-P1(end,:,:)) + ...
    cat(2,P2(:,1,:),diff(P2,1,2),-P2(:,end,:)) + ...
    cat(3,P3(:,:,1),diff(P3,1,3),-P3(:,:,end));

% Compute problem size
[m,n,k]=size(Xobs);

% Preallocate
R1 = zeros(m-1,n,k);
R2 = zeros(m,n-1,k);
R3 = zeros(m,n,k-1);

P1 = R1;
P2 = R2;
P3 = R3;

% Initialize
tkp1=1;
count=0;
i=0;

D=zeros(m,n,k);
fun_all=zeros(1,MAXITER);
st = 12;
while((i<MAXITER)&&(count<5))
    %%%%%%%%%
    % updating the iteration counter
    i=i+1;
    %%%%%%%%%
    % Storing the old value of the current solution
    Dold=D;
    %%%%%%%%%%
    %Computing the gradient of the objective function
    Pold1 = P1;
    Pold2 = P2;
    Pold3 = P3;
    
    tk=tkp1;  % Updtae nesterov coefficient
    
    D=project(Xobs-lambda*Ltv3(R1,R2,R3));  
    
    %Adjoint (this is the adjoint of Ltv3)
    Q1 = -diff(D,1,1);
    Q2 = -diff(D,1,2);
    Q3 = -diff(D,1,3);
    %%%%%%%%%%
    % Taking a step towards minus of the gradient
    P1 = R1 + Q1/(st*lambda);
    P2 = R2 + Q2/(st*lambda);
    P3 = R3 + Q3/(st*lambda);
    %%%%%%%%%%
    % Peforming the projection step

    A = cat(1,P1,zeros(1,n,k)).^2 + cat(2,P2,zeros(m,1,k)).^2 + cat(3,P3,zeros(m,n,1)).^2;
    A = sqrt(max(A,1));
    P1 = P1./A(1:end-1,:,:);
    P2 = P2./A(:,1:end-1,:);
    P3 = P3./A(:,:,1:end-1);

    %%%%%%%%%%
    %Updating R and t
    tkp1=(1+sqrt(1+4*tk^2))/2; 
    
    R1 = P1 + (tk-1)/tkp1*(P1-Pold1);
    R2 = P2 + (tk-1)/tkp1*(P2-Pold2);
    R3 = P3 + (tk-1)/tkp1*(P3-Pold3);
    re=norm(D(:)-Dold(:))/norm(D(:));
    if (re<epsilon)
        count=count+1;   %Only stop if residual has been blow epsilon more than once        
    else
        count=0;
    end
    
    %This block compute function values if desired. 
    %PC=project(Xobs-lambda*Ltv3(P1,P2,P3));
    %PC=project(C);
    %fval=-sum((C(:)-PC(:)).^2)+sum(C(:).^2);
    %fun_all(i) = fval;

end
%fun_all = fun_all(1:i);
X_den=D;
iter=i;
%fprintf('stopping TV at %i iters \n',iter);

