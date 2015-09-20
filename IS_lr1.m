% inputting variant number and determining selection size and number of
% properties

V = input('Input your variant\n');
if(V<7)
  S = 10*V;
  N = 5*V;
else
    if(V<10)
        N = 4*V;
        S = 10*V;
    else
        if(V<20)
            N = 3*V;
            S = 5*V;
        else
            N = 2*V;
            S = 3*V;
        end;
    end;
end;
x = zeros(S,N);
y = zeros(S,1);
K = 2;

%populating teaching selection <x,y>
c = clock;
rng(c(6)*100);
firstExpr = true;
for i=1:1:S
 for j=1:1:N
     if(mod(j,2)==0)
        x(i,j) = 0.01*j/V+0.3*i;    
     else
         if(firstExpr)
         x(i,j) = j * V - 0.1*i;
         firstExpr = false;   
         else
             x(i,j) = j*rand(1);
             firstExpr = true;
         end;
     end;    
 end;
 
 if((x(i,1)*x(i,1)+x(i,2)*x(i,2))<(V*V + 0.04*S*S))
    y(i) = 0;
 else
     y(i) = 1;
 end;
 
end;

%finding centers of each property of each class
tic
centers = zeros(K, N);
classSizes = zeros(K,1);
classSizes(2) = nnz(y);
classSizes(1) = S - classSizes(2);

for i=1:1:K
    for j=1:1:N
        for k=1:1:S
             if(y(k)==i-1)
                    centers(i,j) = centers(i,j)+(x(k,j)/classSizes(K));
             end;
        end;
    end;
end;
disp('Finding centers...');
toc

%setting the selection to be identified, randomly

x_s = zeros(S, N);

firstExpr = true;
for i=1:1:S
 for j=1:1:N
     if(mod(j,2)==0)
        x_s(i,j) = 0.01*j/V+0.3*i;    
     else
         if(firstExpr)
         x_s(i,j) = j * V - 0.1*i;
         firstExpr = false;   
         else
             x_s(i,j) = j*rand(1);
             firstExpr = true;
         end;
     end;    
 end;
 
 
end;
 
%determining range to the center of each class
tic
R = zeros(S,K);
for i=1:1:S
    for k=1:1:K
        for j=1:1:N
               R(i,k) = R(i,k) +(x_s(i,j)-centers(k,j))*(x_s(k,j)-centers(k,j));
        end; 
        R(i, k) = sqrt(R(i,k));
    end;
end;


%finding index of the minimum range, which is the identified class
y_s = zeros(S,1);
for i=1:1:S
    minimalIndices = find(R(i, 1:K)==min(R(i,1:K)),K);
    if(length(minimalIndices)>1)
        y_s(i) = max(classSizes)-1;
    else
    y_s(i) = minimalIndices-1;
    end;
end;
disp('Identifying classes...');
toc

%finding error and error probability
E = 0;
for i=1:1:S
 if((x_s(i,1)*x_s(i,1)+x_s(i,2)*x_s(i,2))<(V*V + 0.04*S*S))
    y_compared = 0;
 else
    y_compared = 1;
 end;
 if(y_compared~=y_s(i))
     E = E+1;
 end;
end;
P = E/S;

disp('Error : ');
disp(E);
disp('Probability of error : ');
disp(P);
disp('Probability of success : ');
disp(1-P);

% identifying instances with only one property which number is V
disp('Now identifying instances with only one property...');
tic
R1 = zeros(S,K);


for j=1:1:K
    for i=1:1:S
    R1(i,j) = sqrt((x_s(i,V)-centers(j, V))*(x_s(i,V)-centers(j, V)));
    end;
end;

y_s1 = zeros(S,1);
for i=1:1:S  
    minimalIndices = find(R1(i,1:K)==min(R1(i,1:K)),K);
    if(length(minimalIndices)>1)
        y_s1(i) = max(classSizes)-1;
    else
    y_s1(i) = minimalIndices-1;
    end;
end;
toc

E1 = 0;
for i=1:1:S
 if((x_s(i,1)*x_s(i,1)+x_s(i,2)*x_s(i,2))<(V*V + 0.04*S*S))
    y_compared = 0;
 else
    y_compared = 1;
 end;
 if(y_compared~=y_s1(i))
     E1 = E1+1;
 end;
end;
P1 = E1/S;

disp('Error : ');
disp(E1);
disp('Probability of error : ');
disp(P1);
disp('Probability of success : ');
disp(1-P1);