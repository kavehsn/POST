
tic;

B=500; %Number of Monte Carlo Iterations
MCDist=zeros(B,1);

load('EmpiricalChapter1.mat')

%Split the sample
x_Est1=DivYield(1:round(0.1*(length(DivYield))));
x_Test1=DivYield(round(0.1*(length(DivYield)))+1:end);

y_Est=AnnExRet(1:round(0.1*(length(AnnExRet))));
y_Test=AnnExRet(round(0.1*(length(AnnExRet)))+1:end);

lagX_Est1=x_Est1(1:end-1);
lagX_Test1=x_Test1(1:end-1);
lagX_1=DivYield(1:end-1);

X_Est=[ones(length(lagX_Est1),1) lagX_Est1];
X_Test=[ones(length(lagX_Test1),1) lagX_Test1];
X=[ones(length(lagX_1),1) lagX_1];


% Estimation of beta by OLS */
betahat=((X_Est')*X_Est)\(X_Est')*y_Est(2:end);

%%Simulated Annealing Algorithm


%Set the parameters


C=0.25; %Controls how fast V adjusts
Beta_0=betahat; %Starting vector of Parameters
V=[0.05 0.05]; %Should cover the entire range of interest in Beta
vareps=0.001; %Convergence Criteria
r_t=0.75; %Temperature Reduction factor
T_0=50; %Initial Temperature
N_vareps=10; %No. of times tolerance is achieved before termination
N_s=20; %No. of times through function before V adjustment
N_t=20; %No. of times through N_s loop before T reduction


%Evalue g(b) at those parameters
g_b=Beta_0(1);

Beta_opt=Beta_0;
g_b_opt=g_b;


counter_vareps=0;

a=0;
while a==0
    
    counter_t=0;
    while counter_t<N_t
        
        
        counter_s=0;
        while counter_s<N_s
            
            counter_acceptance=[0;0];
            counter_rejection=[0;0];
            
            for i=1:numel(Beta_0)
                
                Beta_Trial=Beta_0;
                
                Beta_Trial(i)=Beta_0(i)+unifrnd(-1,1)*V(i);
                
                g_b_dash=Beta_Trial(1);
                
                p=normcdf(X_Test*(betahat-Beta_Trial));
                
                %The signs
                s_y=double(y_Test(2:end)-(X_Test*Beta_Trial)>=0);
                
                %Our Nonparametric Test
                a2=(1./p)-1;
                Test=log(1./a2).*s_y;
                Test=sum(Test);
                
                for k=1:B
                    
                    
                    mu_sim=0;
                    sigma_sim=1;
                    eps_sim=normrnd(mu_sim,sigma_sim,round(0.9*(length(DivYield))),1);
                    
                    
                    y_sim=(X_Test*Beta_Trial)+eps_sim(2:end);
                    
                    %The signs
                    s_ysim=double(y_sim-(X_Test*Beta_Trial)>=0);
                    p_sim=normcdf(X_Test*(betahat-Beta_Trial));
                    
                    %Our Nonparametric Test
                    a2_s=(1./p_sim)-1;
                    Test_s=log(1./a2_s).*s_ysim;
                    Test_s=sum(Test_s);
                    
                    MCDist(k)=Test_s;
                    
                end
                
                UB=quantile(MCDist,0.95);
                LB=quantile(MCDist,0.025);
                
                if Test<UB && Test>LB && isnan(UB)==0 && isnan(LB)==0
                    
                    if g_b_dash<=g_b
                        
                        p=exp( (g_b_dash - g_b) / T_0 );
                        
                        p_dash=rand;
                        
                        if p>p_dash
                            
                            Beta_0=Beta_Trial;
                            g_b=g_b_dash;
                            
                            counter_acceptance(i)=counter_acceptance(i)+1;
                        else
                            
                            counter_rejection(i)=counter_rejection(i)+1;
                            
                        end
                        
                    end
                    if g_b_dash>g_b
                        
                        Beta_0=Beta_Trial;
                        
                        g_b=g_b_dash;
                        counter_acceptance(i)=counter_acceptance(i)+1;
                        
                    end
                    
                    if g_b_dash>g_b_opt
                        
                        if g_b_dash-g_b_opt<vareps
                            
                            counter_vareps=counter_vareps+1;
                            
                        end
                        
                        Beta_0=Beta_Trial;
                        g_b=g_b_dash;
                        Beta_opt=Beta_Trial;
                        g_b_opt=g_b_dash;
                        
                    end
                end
                
                
            end
            counter_s=counter_s+1;
        end
        
        if (counter_acceptance+counter_rejection)~=0
           
            
            counter_avg=counter_acceptance./(counter_acceptance+counter_rejection);
            
            
            for j=1:numel(counter_avg)
                
                if counter_avg(j)>0.5
                    
                    V(j)=(1+C)*V(j);
                    
                end
                
            end
        end
        
        counter_t=counter_t+1;
        
    end
    
    if counter_vareps>N_vareps && abs(g_b-g_b_dash)<vareps
        
        disp(['The optimal value of Beta is: ',num2str(Beta_opt(1)),' and ',num2str(Beta_opt(2)),' and ',num2str(Beta_opt(3))]);
        disp(['The optimal value of function is: ',num2str(g_b_opt)]);
        disp(['The speed value is: ',num2str(V(1))]);
        
        Power=[Beta_opt;g_b_opt];
        
        break
        
    else
        
        Beta_0=Beta_opt;
        T_0=r_t*T_0;
        
    end
    
    disp(['The temperature is: ',num2str(T_0),' and the Iterations left are ',num2str(N_vareps-counter_vareps)])
    disp(['Function Difference is ',num2str(abs(g_b-g_b_dash))]);
    disp(['New Optimum Beta is= ',num2str(Beta_opt(1)),' and ',num2str(Beta_opt(2))]);
    
end

pp=toc;
disp(['Code took ',num2str(pp/60), ' minutes to run!'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if results are stored in MATLAB variable(s) e.g. 'OUT' these could be saved to ASCII text file(s), with
% base name e.g. output and labeled by the data set index. The file name is stored in the variable output_file as

output_file = strcat( 'DivYieldMaxBeta_0.txt' );

% and the value of OUT saved

save( output_file , 'Power' ,  '-ascii' );

%

%
