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
y2 = zeros(S,1);
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
 
 y2(i) = 2*x(i,1)+0.1*x(i,2);

end;

% coefficient of sing correlation
sign_cor = zeros(S,1);

x_mean = zeros(N,1);
y_mean = 0;

for i=1:1:N
    x_mean(i) = mean(x(1:S, i));
end;

y_mean = mean(y);

c_x_y = zeros(N,1);
c_x = zeros(N,1);
c_y = 0;

for i=1:1:N % calculating c_x_y , calculating c_x 
    for j=1:1:S
        if( (x(j,i)-x_mean(i))>=0 && (y(j)-y_mean)>=0)
            c_x_y(i)= c_x_y(i) +1;
        end; 
        if( (x(j,i)-x_mean(i))>=0 )
            c_x(i)= c_x(i) +1;
        end; 
    end;
    c_x_y(i)= c_x_y(i)/S;
    c_x(i) = c_x(i)/S;
end;

for i=1:1:S % calculating c_y 
    if( (y(i)-y_mean)>=0)
        c_y = c_y+1;
    end;
end;
c_y= c_y/S;

for i=1:1:N %calculating coefficient of sign correlation itself
    sign_cor(i) = (c_x_y(i)-c_x(i)*c_y)/sqrt(c_x(i)*c_y*(1-c_x(i))*(1-c_y));
end;

% Fechner correlation coefficient
fechner = zeros(N,1);
c_f = zeros(N,1);
d_f = zeros(N,1);

for i=1:1:N % calculating c_f , calculating d_f
    for j=1:1:S
        if( ((x(j,i)-x_mean(i))>=0 && (y(j)-y_mean)>=0) || ((x(j,i)-x_mean(i))<0 && (y(j)-y_mean)<0) )
            c_f(i)= c_f(i)+1;
        else
            d_f(i) = d_f(i)+1;
        end; 
    end;
end;

for i=1:1:N %calculating Fechner correlation coefficient itself
    fechner(i) = (c_f(i)-d_f(i))/(c_f(i)+d_f(i));
end;


% coefficient of sing correlation (FOR Y2)
sign_cor_2 = zeros(S,1);

x_mean = zeros(N,1);
y_mean = 0;

for i=1:1:N
    x_mean(i) = mean(x(1:S, i));
end;

y2_mean = mean(y2);

c_x_y = zeros(N,1);
c_x = zeros(N,1);
c_y = 0;

for i=1:1:N % calculating c_x_y , calculating c_x 
    for j=1:1:S
        if( (x(j,i)-x_mean(i))>=0 && (y2(j)-y2_mean)>=0)
            c_x_y(i)= c_x_y(i) +1;
        end; 
        if( (x(j,i)-x_mean(i))>=0 )
            c_x(i)= c_x(i) +1;
        end; 
    end;
    c_x_y(i)= c_x_y(i)/S;
    c_x(i) = c_x(i)/S;
end;

for i=1:1:S % calculating c_y 
    if( (y2(i)-y2_mean)>=0)
        c_y = c_y+1;
    end;
end;
c_y= c_y/S;

for i=1:1:N %calculating coefficient of sign correlation itself
    sign_cor_2(i) = (c_x_y(i)-c_x(i)*c_y)/sqrt(c_x(i)*c_y*(1-c_x(i))*(1-c_y));
end;

% Fechner correlation coefficient
fechner_2 = zeros(N,1);
c_f = zeros(N,1);
d_f = zeros(N,1);

for i=1:1:N % calculating c_f , calculating d_f
    for j=1:1:S
        if( ((x(j,i)-x_mean(i))>=0 && (y2(j)-y2_mean)>=0) || ((x(j,i)-x_mean(i))<0 && (y2(j)-y2_mean)<0) )
            c_f(i)= c_f(i)+1;
        else
            d_f(i) = d_f(i)+1;
        end; 
    end;
end;

for i=1:1:N %calculating Fechner correlation coefficient itself
    fechner_2(i) = (c_f(i)-d_f(i))/(c_f(i)+d_f(i));
end;