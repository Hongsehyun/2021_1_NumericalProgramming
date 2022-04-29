function Pcoef = myPolyfit_student(x, y,n)
% Calculates the coefficients of polynomial that best fits n data points.
% This is for explanation purpose. You should make it more efficient in
% terms of reducing iteration numbers
%
% Input variables:
%   x: A vector with the coordinates x of the data points.
%   y: A vector with the coordinates y of the data points.
% Output variables:
%   Pcoef=[....a4 a3 a2 a1 a0] 

%% Your code goes here

%  (n)th order polynomial:  (n+1) unknowns
% A: (m) by (n+1) //  A'A: (n+1) by (n+1)
% Polynomial order ...+(a4)x^4+(a3)x^3+...+(a1)x+a0.

    %   1-1. Check if length(X)~= length(Y) ? Exit: Continue 
    %   1-2. Check if   m<(n+1) ? Exit: Continue 
    mx = length(x);
    my = length(y);  
    % ADD YOUR CODE         
    % ADD YOUR CODE  
    
    
    %   2. Initialize S, b, z
    Pcoef=zeros(1,n+1);
    S=zeros(n+1,n+1);
    b=zeros(n+1,1);
    z=zeros(1,n+1);

    
    % 3.  SX(i)=sum(x.^(i)), for i=0 to 2n  #(2n+1)
    % ADD YOUR CODE         
    % ADD YOUR CODE        
    % ADD YOUR CODE        
    % ADD YOUR CODE     
  
    
    %   4. Create matrix A'A from SX
    % ADD YOUR CODE        
    % ADD YOUR CODE   
    
    
    %  5. Construct vector b=(A'y) as
    % ADD YOUR CODE         
    % ADD YOUR CODE        
    
    %  6. Calculate for optimal Z=inv(S)*(b)
    % where z=[a(0), a(1), a(2),.... a(n)] , (n+1)x1
    % ADD YOUR CODE         
    % ADD YOUR CODE    
    
    %  7-1.  Return z=[a0, a1, ....];
    %  7-2   If need to return coefficent as    Pcoef=[....a4 a3 a2 a1 a0]
    for i = 1:n+1
        Pcoef(i) = z(n+1-i+1);
    end


end