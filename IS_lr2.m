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

%finding clusters for not normalized selection
[distinct_clusters, s] = subclust(x, 0.5);
[fuzzy_clusters, u] = fcm(x, 3, [NaN, NaN, NaN, 0]);
disp('Displaying distinct clusters for not normalized selection :');
disp(distinct_clusters);
disp('Displaying fuzzy clusters for not normalized selection :');
disp(fuzzy_clusters);

%finding clusters for normalized selection
x_n = x/norm(x);
[distinct_clusters_n, s_n] = subclust(x_n, 0.5);
[fuzzy_clusters_n, u_n] = fcm(x_n, 3, [NaN, NaN, NaN, 0]);
disp('Displaying distinct clusters for normalized selection :');
disp(distinct_clusters_n);
disp('Displaying fuzzy clusters for normalized selection :');
disp(fuzzy_clusters_n);
