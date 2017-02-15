function AKF = POCAfunc(xinit,N,P,epsilon,a)

% DO NOT CHANGE THE FILE
% NOT randominzed
% a is the PAPR circle value
if nargin<4
    epsilon = .2*10^-14;           % do not change this value
end
if nargin<3
    P = 30;
end
        Q = P;                  %We just simulate the P = Q value. 

T = eye(P);
ind = [1:P]  ;        % is determined from the T matrix

%Golomb Siquence
golomb = zeros(N,1);
for n = 1:N
        golomb(n) = exp(1i*pi*(n-1)*n/N);
end
 x = xinit;
% x = rand(N,1); 
NumberOfIterations = 0;

 
 %Constructing Xhat
Xhat = zeros(N+P-1,P);
for i = 1:N+P-1
    for j = 1:P
        if ((i-j+1 <= N) && (i-j+1>0))
            Xhat(i,j) = x(i-j+1);
        end
%     end
end
end

XtildaNew = Xhat*T;
Xtilda = zeros(N+P-1,Q);
x = zeros(N,1);
 tic
%   while norm(Xtilda(:)-XtildaNew(:),'inf')>epsilon
for i = 1:10000

 
      
    
      %For debugging
%       NumberOfIterations = NumberOfIterations+1
      norm(Xtilda-XtildaNew,'inf')
  % [ISL] = sidelobe(x);
   
     Xtilda = XtildaNew;
     
     
     [U1, Sigma, U2] = svd(Xtilda','econ');                 %equation 16
     U = U2*U1';                                            %equation 17
     
     clear mu_k
     %% Step 2 
     for n = 1:N
         S = 0;
         for i = 1:Q
               mu_k(i) = (sqrt(N))*U(n-1+ind(i),i);
               
         end

             xxx =   (max(mu_k)+min(mu_k))/2;             %equation 20
             if (a ==1 || abs(xxx)>a)
                x(n) = a*xxx./abs(xxx);
             else
                 x(n) = xxx;
             end
           
                 
                 
     end
   
     % Constructing Xtilda from x
    Xhat = zeros(N+P-1,P);
    for i = 1:N+P-1
        for j = 1:P
            if ((i-j+1 <= N) && (i-j+1>0))
                Xhat(i,j) =x(i-j+1);
            end
        end
    end
    XtildaNew = Xhat*T;
    
    
   
  end
 toc
 
 %% Computing Autocorrelation function
 AKF = zeros(2*N-1,1);
 for k = 0:N-1
    S = 0;
    for  n = k+1:N
        S = S+x(n)*x(n-k)';
    end
        AKF(N+k) = S;
        AKF(N-k) = S';
 end     
AKF = AKF/AKF(N);




 