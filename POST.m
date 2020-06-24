tic;
T=500; %Set the number of observations
theta=0.99; %Autocorrelation coefficient
rho=-0.8; %Contemporaneous correlation coefficient
m_1=2000; %Number of iterations for simulating the distribution of the test statistic 
m_2=1000; %Number of iterations for simulating the power of the test
mp=20;
Power_POS=zeros(mp,1);
Power_T=zeros(mp,1);
Power_WT=zeros(mp,1);
Power_CD=zeros(mp,1);
beta=0;
counter_beta=1;

while counter_beta<=mp
    
    Reject=0;
    counterPower_T=0;
    counterPower_WT=0;
    counterPower_CD=0;
    
    for j=1:m_2
        
        %Generate the data within a predictive regression context
        
        y=zeros(T,1); 
        x=zeros(T,1);
        
        mu=0;
        sigma=1;
        w=normrnd(mu,sigma,T,1);
        u=zeros(T,1);
        eps=normrnd(0,1,T,1); %distribution of the residuals (change for Cauchy, Student's t, etc).
        POS_sim=zeros(m_1,1);
        
        x(1)=w(1)/sqrt(1-theta^2);
        
        for t=2:T
            y(t)=beta*x(t-1)+eps(t);
            u(t)=rho*eps(t)+sqrt(1-rho^2)*w(t);
            x(t)=theta*x(t-1)+u(t);
        end
        
        %Use the 10% split sample technique to prevent data snooping bias
        y1=y(2:round(0.1*T));
        x1=x(2:round(0.1*T));
        
        y2=y(round(0.1*T)+1:T);
        x2=x(round(0.1*T)+1:T);
        
        lagX=[ones(T-1,1) x(1:end-1)];
        
        % Estimation of beta by OLS */
        betahat=((lagX')*lagX)\(lagX')*y(2:end);
        betahatPOS=regress(y1(2:end),x1(1:end-1));
        
        %Find the sign function of y
        s=double(y2>=0);
        
        %Calculate the test-statistic
        Phi=normcdf(x2(1:end-1)*betahatPOS);
        POS=sum(log(Phi./(1-Phi)).*s(2:end));
        
        %Simulate the distribtution of the test-statistic
        for t=1:m_1
            
            s_sim=binornd(1,0.5,length(y2),1);
            POS_sim(t,1)=sum(log(Phi./(1-Phi)).*s_sim(2:end));
            
        end
        
        %find the 5% critical value
        CritVal=quantile(POS_sim,0.95);
        
        if POS>=CritVal
            
            Reject=Reject+1;
            
        end
        
        %Other Tests->
        
        %NonParametric Test of CD (1995)
        CD=y(2:end).*lagX(:,2);
        S_CD=double(CD>=0);
        Test_CD=sum(S_CD);
        
        
        if Test_CD>binoinv(0.95,T,0.5)
            
            counterPower_CD=counterPower_CD+1;
            
        end
        
        
        % Estimation of variance of betahat */
        uhat=y(2:end)-lagX*betahat;
        seghat=((uhat')*uhat)/(T-2);
        segbhat=seghat/(((lagX(:,2)')*lagX(:,2)));
        
        
        % Calculation of robust standard deviation of White */
        moy=(1/(T-1))*sum(lagX(:,2));
        SST=sum((lagX(:,2)-moy.*ones(T-1,1)).^2);
        var3=(((lagX(:,2)-moy.*ones(T-1,1)).^2)')*  (uhat.^2);
        Varob=var3/(SST^2);
        
        
        
        if betahat(2)/sqrt(segbhat)> 1.675
            counterPower_T=counterPower_T+1;
        end
        
        % Decision rule of testing H0: beta=0  using corrected white T-test*/
        
        if betahat(2)/sqrt(Varob)> 1.675
            counterPower_WT=counterPower_WT+1;
        end
        
        
    end
    
    
    beta=0.008*counter_beta;
    Power_POS(counter_beta)=Reject/m_2;
    Power_T(counter_beta)=counterPower_T/m_2;
    Power_WT(counter_beta)=counterPower_WT/m_2;
    Power_CD(counter_beta)=counterPower_CD/m_2;
    
    counter_beta=counter_beta+1;
    disp(counter_beta)
    
end
z=toc;

disp(['The code took ',num2str(z/60),' to run!'])