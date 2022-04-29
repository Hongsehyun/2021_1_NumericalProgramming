function z= linearFit_student(x, y)
% LinearRegration calculates the coefficients a1 and a0 of the linear
% equation y = a1*x + a0 that best fit n data points.
% Input variables:
% x    A vector with the coordinates x of the data points.
% y    A vector with the coordinates y of the data points.
% Output variable:
% a1   The coefficient a1.
% a0   The coefficient a0.

%% Your code goes here

    % Input data sets of (x, y):  # m data sets
 
    %   1-1. Check if the number of dataset is less than 2. (  m<2 )
    %   1-2. Check if length(X)~= length(Y) ? Exit: Continue 
    mx = length(x);
    my = length(y);  

    % ADD YOUR CODE         
    % ADD YOUR CODE         
    
    %   2. Initialize Sx, Sxx, Sy, Sxy
    Sx=0; Sxx=0; Sxy=0;Sy=0;  
    a1=0; a0=0; 
    
    %   3. Solve for Sx, Sxx, Sy, Sxy,  for k=1 to m
    % ADD YOUR CODE         
    % ADD YOUR CODE        
    
    %   4. Solve for a1, a2
    % ADD YOUR CODE         
    % ADD YOUR CODE        
    
    %   5. Return z=[a1, a0]
    
    z=[a1, a0];
    
end   % end of function