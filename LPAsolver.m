%% Linearized Proximal Algorithm with semismooth Newton iteration (LPA-SN)
%%
%% min h(F(x)) 
%% h : R^m -> R convex function   
%% F(x) : R^n -> R^m continuously differentiable
%%
%% LPAsolverSN
%%
%% X0 : unknown sensor position
%% P0 : anchor position
%% DD : sensor-sensor & sensor-anchor distance

function [X0, rmsd] = LPAsolver(P0,PP,DD,dim,num,R)
M = 1;
alpha = 2;
p=2;
dp = rand(dim*num,1); %(2*100,1)
v = 100;
%X0 = -0.5+rand(dim,num);
X0 = PP;
count = 0;
method=2;
%% Outer iteration
while count < 30
    count= count+1;
    F=zeros(size(X0,2)+size(P0,2),size(X0,2));
    directF=zeros((size(X0,2)+size(P0,2))*size(X0,2),1);
    flag=zeros(size(X0,2)+size(P0,2),size(X0,2));
    dF=zeros((size(X0,2)+size(P0,2))*size(X0,2),dim*size(X0,2)); 
    
%% Calculate F(x)
    for i =1:size(X0,2)     
        for j =1:size(X0,2)           
            F(j,i) = norm(X0(:,i)-X0(:,j))^2; %distance(X0(:,i),X0(:,j));         
            if DD(i,j)~=0 && i~=j
                if F(j,i)-DD(i,j)^2>=0
                    F(j,i)=F(j,i)-DD(i,j)^2;
                    flag(j,i)=1;
                else
                    F(j,i)=DD(i,j)^2-F(j,i);                   
                end
            else
                if i==j
                    F(j,i)=0;
                else
                    F(j,i)= R*R-F(j,i);
                end
            end
        end
        
        for j = 1:size(P0,2)
            F(j+length(X0),i) = distance(X0(:,i),P0(:,j));
            if DD(i,length(X0)+j)~=0     
                if F(j+length(X0),i)-DD(i,length(X0)+j)^2>=0
                F(j+length(X0),i)= F(j+length(X0),i)-DD(i,length(X0)+j)^2;
                flag(j+length(X0),i)=1;
                else
                    F(j+length(X0),i)= DD(i,length(X0)+j)^2-F(j+length(X0),i);
                end              
            else
                F(j+length(X0),i)= R*R-F(j+length(X0),i);
            end
        end
    end
    
%% Reform F(x)    
    for i=1:100
        directF(110*(i-1)+1:110*i)=F(:,i);
    end
    
%% Calculate F'(x)
    for i=1:11000
        t=floor((i-1)/110)+1;
        k=mod(i-1,110)+1;
        if k<=100
            if flag(k,t)==0
                dF(i,2*t-1)=-2*(X0(1,t)-X0(1,k));
                dF(i,2*t)=-2*(X0(2,t)-X0(2,k));
                dF(i,2*k-1)=2*(X0(1,t)-X0(1,k));
                dF(i,2*k)=2*(X0(2,t)-X0(2,k));
            else
                dF(i,2*t-1)=2*(X0(1,t)-X0(1,k));
                dF(i,2*t)=2*(X0(2,t)-X0(2,k));
                dF(i,2*k-1)=-2*(X0(1,t)-X0(1,k));
                dF(i,2*k)=-2*(X0(2,t)-X0(2,k));
            end
        else
            if flag(k,t)==0
                dF(i,2*t-1)=-2*(X0(1,t)-P0(1,k-100));
                dF(i,2*t)=-2*(X0(2,t)-P0(2,k-100));
                
            else
                dF(i,2*t-1)=2*(X0(1,t)-P0(1,k-100));
                dF(i,2*t)=2*(X0(2,t)-P0(2,k-100));
            end
        end
    end
    
    %% Calculate dp
    dFt = dF';
    H0=dFt*directF;
    %    tmp=dF*dp+directF;
    %    tmp(tmp<0)=0;
    iteration=1;
    
    if method==1
        while norm(dFt*max((dF*dp+directF),0)+dp/v)>=max(norm(dp)^3,1e-6) && iteration<=50
            
            %        if norm(H0)<=M*norm(dp)^alpha
            %            dp=(norm(dp)^(alpha-1))*dp;
            %        else
            dp=myFun(dFt,dF,directF,v,dp);
            iteration=iteration+1
        end
    else
        u=((1/2)*norm(max(directF,0))^p)^alpha;                                                                                     %
        while   iteration<=1             %norm(dFt*max((dF*dp+directF),0)+2*u*dp)>=max(norm(dp)^3,1e-6)&&
            
            %        if norm(H0)<=M*norm(dp)^alpha
            %            dp=(norm(dp)^(alpha-1))*dp;
            %        else
            dp=myFun(dFt,dF,directF,u,dp);
            iteration=iteration+1
            
        end
        
        
    end
    
    
    
    
    for i=1:size(X0,2)        
        Dp(:,i)=dp(2*i-1:2*i,1);
    end
    X0=X0+Dp;
    rmsd(count) = RMSD(X0,PP)
end
end



    
  
