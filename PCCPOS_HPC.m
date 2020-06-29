tic;

% create a local cluster object
pc = parcluster('local')

% explicitly set the JobStorageLocation to the temp directory that was created in your sbatch script
pc.JobStorageLocation = strcat('slurm_task_id_/',getenv('SLURM_ARRAY_TASK_ID'))

% start the matlabpool with maximum available workers
% control how many workers by setting ntasks in your sbatch script
parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')))


Coef=0;
CurveCounter=1;
Power=zeros((0.5/0.05)+1,5);

while Coef<=0.5
    
    Iter=500;
    counterPower=0;
    counterPower_POS=0;
    counterPower_CD=0;
    counterPower_WT=0;
    counterPower_T=0;
    Round=0;
    
    for ll=1:Iter
        
        count = 0;
        err_count = 0;
        while count == err_count
            
            try
                
                T=50; %Set the number of observations
                Tr=0;
                B=500; %Number of Monte Carlo Iterations
                y=zeros(T,1); %Create a vector of dim Tx2 to store the simulated dependents
                x=zeros(T,1);
                Measures=zeros(T-1,2);
                MCDist=zeros(B,1);
                MCDist_Old=zeros(B,1);
                
                mu=0;
                sigma=1;
                w=normrnd(mu,sigma,T,1);
                u=zeros(T,1);
                eps=trnd(10,T,1);
                %True Data Generating Process
                
                Intercept=0;
                Slope=Coef;
                rho=0;
                theta=0.9;
                
                x(1)=w(1)/sqrt(1-theta^2);
                u(1)=rho*eps(1)+w(1)*sqrt(1-rho^2);
                
                for t=2:T
                    y(t)=Intercept+Slope*x(t-1)+eps(t);
                    u(t)=rho*eps(t)+w(t)*sqrt(1-rho^2);
                    x(t)=theta*x(t-1)+u(t);
                end
                
                %The signs
                s_y=double(y>=0);
                
                lagX=x(1:end-1);
                
                X=[ones(length(lagX),1) lagX];
                
                
                Beta=[0;Coef];
                
                if Coef==0
                Beta_CFC=robustfit(lagX,y(2:end));
                else
                Beta_CFC=Beta;
                end
                
                p=tcdf(X*Beta,1);
                
                
                CDF_P=binocdf(s_y(2:end),1,p);
                CDF_N=binocdf(s_y(2:end)-1,1,p);
                Pdf=CDF_P-CDF_N;
                
                 S_star=s_y(2:end)+(rand(numel(s_y(2:end)),1)-1);
                
                
                %ALL_1=[CDF_P(2:end) CDF_N(1:T-1) Pdf(2:end)];
                
                %for t=1:3
                
                %indices=find(ALL_1(:,t)==0);
                %ALL_1(indices,:)=[];
                %indices=find(ALL_1(:,t)==1);
                %ALL_1(indices,:)=[];
                
                %end
                
                %CDF_P=ALL_1(:,1);
                %CDF_N=ALL_1(:,2);
                %Pdf=ALL_1(:,3);
                
                C_pp_tt1=zeros(T-2,1);
                
                CDDF_P=CDF_P(1:end-1);
                CDDF_P1=CDF_P(2:end);
                
                
                Tau_pp_tt1=kendalltau([S_star(2:end),S_star(1:end-1)]);
                A_pp_tt1=copulaparam('Gaussian',Tau_pp_tt1(1,2));
                
                for i=1:T-1
                    if i+1<=T-1
                        if CDF_P(i)==1
                            C_pp_tt1(i)=CDF_P(i+1);
                        elseif CDF_P(i+1)==1
                            C_pp_tt1(i)=CDF_P(i);
                        elseif CDF_P(i)==0 || CDF_P(i+1)==0
                            C_pp_tt1(i)=0;
                        else
                            C_pp_tt1(i)=copulacdf('Gaussian',[CDF_P(i) CDF_P(i+1)],A_pp_tt1);
                        end
                    end
                end
                
                
                C_pn_tt1=zeros(T-2,1);
                
               
                CDDF_N1=CDF_N(2:end);
                
                Tau_pn_tt1=kendalltau([S_star(2:end)-1,S_star(1:end-1)]);
                A_pn_tt1=copulaparam('Gaussian',Tau_pn_tt1(1,2));
                
                
                for i=1:T-1
                    if i+1<=T-1
                        if CDF_P(i)==1
                            C_pn_tt1(i)=CDF_N(i+1);
                        elseif CDF_N(i+1)==1
                            C_pn_tt1(i)=CDF_P(i);
                        elseif CDF_P(i)==0 || CDF_N(i+1)==0
                            C_pn_tt1(i)=0;
                        else
                            C_pn_tt1(i)=copulacdf('Gaussian',[CDF_P(i) CDF_N(i+1)],A_pn_tt1);
                        end
                    end
                end
                
                C_np_tt1=zeros(T-2,1);
                
                CDDF_N=CDF_N(1:end-1);
                CDDF_P1=CDF_P(2:end);
                
                Tau_np_tt1=kendalltau([S_star(2:end),S_star(1:end-1)-1]);
                A_np_tt1=copulaparam('Gaussian',Tau_np_tt1(1,2));
                
                
                for i=1:T-1
                    if i+1<=T-1
                        if CDF_N(i)==1
                            C_np_tt1(i)=CDF_P(i+1);
                        elseif CDF_P(i+1)==1
                            C_np_tt1(i)=CDF_N(i);
                        elseif CDF_N(i)==0 || CDF_P(i+1)==0
                            C_np_tt1(i)=0;
                        else
                            C_np_tt1(i)=copulacdf('Gaussian',[CDF_N(i) CDF_P(i+1)],A_np_tt1);
                        end
                    end
                end
                
                C_nn_tt1=zeros(T-2,1);
                
                CDDF_N=CDF_N(1:end-1);
                CDDF_N1=CDF_N(2:end);
                
                Tau_nn_tt1=kendalltau([S_star(2:end)-1,S_star(1:end-1)-1]);
                A_nn_tt1=copulaparam('Gaussian',Tau_nn_tt1(1,2));
                
                for i=1:T-1
                    if i+1<=T-1
                        if CDF_N(i)==1
                            C_nn_tt1(i)=CDF_N(i+1);
                        elseif CDF_N(i+1)==1
                            C_nn_tt1(i)=CDF_N(i);
                        elseif CDF_N(i+1)==0 || CDF_N(i)==0
                            C_nn_tt1(i)=0;
                        else
                            C_nn_tt1(i)=copulacdf('Gaussian',[CDF_N(i) CDF_N(i+1)],A_nn_tt1);
                        end
                    end
                end
                
                
                
                
                %Step3(a)
                F_p_tgt1=(C_pp_tt1(1:end-1)-C_pn_tt1(1:end-1))./Pdf(2:end-1);
                F_n_tgt1=(C_np_tt1(1:end-1)-C_nn_tt1(1:end-1))./Pdf(2:end-1);
                f_tgt1=F_p_tgt1-F_n_tgt1;
                
                F_p_tgt1(isnan(F_p_tgt1))=0;
                F_n_tgt1(isnan(F_n_tgt1))=0;
                
                F_p_second=(C_pp_tt1(end-1)-C_pn_tt1(end-1))./Pdf(end);
                F_n_second=(C_np_tt1(end-1)-C_nn_tt1(end-1))./Pdf(end);
                f_second=F_p_second-F_n_second;
                
                %Step3(b)
                F_p_t2gt1=(C_pp_tt1(2:end)-C_np_tt1(2:end))./Pdf(2:end-1);
                F_n_t2gt1=(C_pn_tt1(2:end)-C_nn_tt1(2:end))./Pdf(2:end-1);
                f_t2gt1=F_p_t2gt1-F_n_t2gt1;
                
                F_p_t2gt1(isnan(F_p_t2gt1))=0;
                F_n_t2gt1(isnan(F_n_t2gt1))=0;
                
                
                %Step3(c)
                
                
                
                C_pp_tt2gt1=zeros(T-3,1);
                
                
                
                Tau_pp_tt2gt1=kendalltau([S_star(3:end),S_star(1:end-2)]);
                A_pp_tt2gt1=copulaparam('Gaussian',Tau_pp_tt2gt1(1,2));
                
                for i=1:T-3
                    if F_p_tgt1(i)==1
                        C_pp_tt2gt1(i)=F_p_t2gt1(i);
                    elseif F_p_t2gt1(i)==1
                        C_pp_tt2gt1(i)=F_p_tgt1(i);
                    elseif F_p_tgt1(i)==0 || F_p_t2gt1(i)==0
                        C_pp_tt2gt1(i)=0;
                    else
                        C_pp_tt2gt1(i)=copulacdf('Gaussian',[F_p_tgt1(i) F_p_t2gt1(i)],A_pp_tt2gt1);
                    end
                end
                
                
                C_pn_tt2gt1=zeros(T-3,1);
                
                Tau_pn_tt2gt1=kendalltau([S_star(3:end)-1,S_star(1:end-2)]);
                A_pn_tt2gt1=copulaparam('Gaussian',Tau_pn_tt2gt1(1,2));
                
                
                for i=1:T-3
                    if F_p_tgt1(i)==1
                        C_pn_tt2gt1(i)=F_n_t2gt1(i);
                    elseif F_n_t2gt1(i)==1
                        C_pn_tt2gt1(i)=F_p_tgt1(i);
                    elseif F_p_tgt1(i)==0 || F_n_t2gt1(i)==0
                        C_pn_tt2gt1(i)=0;
                    else
                        C_pn_tt2gt1(i)=copulacdf('Gaussian',[F_p_tgt1(i) F_n_t2gt1(i)],A_pn_tt2gt1);
                    end
                end
                
                C_np_tt2gt1=zeros(T-3,1);
                
                Tau_np_tt2gt1=kendalltau([S_star(3:end),S_star(1:end-2)-1]);
                A_np_tt2gt1=copulaparam('Gaussian',Tau_np_tt2gt1(1,2));
                
                
                for i=1:T-3
                    if F_n_tgt1(i)==1
                        C_np_tt2gt1(i)=F_p_t2gt1(i);
                    elseif F_p_t2gt1(i)==1
                        C_np_tt2gt1(i)=F_n_tgt1(i);
                    elseif F_n_tgt1(i)==0 || F_p_t2gt1(i)==0
                        C_np_tt2gt1(i)=0;
                    else
                        C_np_tt2gt1(i)=copulacdf('Gaussian',[F_n_tgt1(i) F_p_t2gt1(i)],A_np_tt2gt1);
                    end
                end
                
                C_nn_tt2gt1=zeros(T-3,1);
                
                Tau_nn_tt2gt1=kendalltau([S_star(3:end)-1,S_star(1:end-2)-1]);
                A_nn_tt2gt1=copulaparam('Gaussian',Tau_nn_tt2gt1(1,2));
                
                
                for i=1:T-3
                    if F_n_tgt1(i)==1
                        C_nn_tt2gt1(i)=F_n_t2gt1(i);
                    elseif F_n_t2gt1(i)==1
                        C_nn_tt2gt1(i)=F_n_tgt1(i);
                    elseif F_n_tgt1(i)==0 || F_n_t2gt1(i)==0
                        C_nn_tt2gt1(i)=0;
                    else
                        C_nn_tt2gt1(i)=copulacdf('Gaussian',[F_n_tgt1(i) F_n_t2gt1(i)],A_nn_tt2gt1);
                    end
                end
                
                C_pp_ttdm1gt1tdm1=zeros(T-3,T-4);
                C_pp_ttdm1gt1tdm1(1:T-3,1)=C_pp_tt2gt1;
                
                
                C_pn_ttdm1gt1tdm1=zeros(T-3,T-4);
                C_pn_ttdm1gt1tdm1(1:T-3,1)=C_pn_tt2gt1;
                
                
                C_np_ttdm1gt1tdm1=zeros(T-3,T-4);
                C_np_ttdm1gt1tdm1(1:T-3,1)=C_np_tt2gt1;
                
                
                C_nn_ttdm1gt1tdm1=zeros(T-3,T-4);
                C_nn_ttdm1gt1tdm1(1:T-3,1)=C_nn_tt2gt1;
                
                f_conda=zeros(T-3,T-4);
                f_conda(1:T-3,1)=f_t2gt1;
                
                f_condb=zeros(T-3,T-4);
                f_condb(1:T-3,1)=f_tgt1;
                
                clear C_pp_tt2gt1 C_pn_tt2gt1 C_np_tt2gt1 C_nn_tt2gt1 f_t2gt1
                
                counter=3;
                
                
                for d=1:T-4
                    
                    
                    %Step 4(a)
                    F_p_tgt1t2=(C_pp_ttdm1gt1tdm1(1:T-counter,d)-C_pn_ttdm1gt1tdm1(1:T-counter,d))./f_conda(1:T-counter,d);
                    F_n_tgt1t2=(C_np_ttdm1gt1tdm1(1:T-counter,d)-C_nn_ttdm1gt1tdm1(1:T-counter,d))./f_conda(1:T-counter,d);
                    
                    F_p_tgt1t2(isnan(F_p_tgt1t2))=0;
                    F_n_tgt1t2(isnan(F_n_tgt1t2))=0;
                    
                    %Step 4(b)
                    F_p_t3gt1t2=(C_pp_ttdm1gt1tdm1(2:T-counter,d)-C_np_ttdm1gt1tdm1(2:T-counter,d))./f_condb(2:T-counter,d);
                    F_n_t3gt1t2=(C_pn_ttdm1gt1tdm1(2:T-counter,d)-C_nn_ttdm1gt1tdm1(2:T-counter,d))./f_condb(2:T-counter,d);
                    
                    F_p_t3gt1t2(isnan(F_p_t3gt1t2))=0;
                    F_n_t3gt1t2(isnan(F_n_t3gt1t2))=0;
                    
                    f_conda(1:T-(counter+1),d+1)=F_p_t3gt1t2-F_n_t3gt1t2;
                    f_condb(1:T-counter,d+1)=F_p_tgt1t2-F_n_tgt1t2;
                    
                    if d<=Tr
                        
                    TauFinal_pp=kendalltau([S_star(counter+1:end),S_star(1:end-counter)]);
                    AFinal_pp=copulaparam('Gaussian',TauFinal_pp(1,2));

                        %Step 4(c)
                        for i=1:T-(counter+1)
                            if F_p_tgt1t2(i)==1
                                C_pp_ttdm1gt1tdm1(i,d+1)=F_p_t3gt1t2(i);
                            elseif F_p_t3gt1t2(i)==1
                                C_pp_ttdm1gt1tdm1(i,d+1)=F_p_tgt1t2(i);
                            elseif F_p_tgt1t2(i)==0 || F_p_t3gt1t2(i)==0
                                C_pp_ttdm1gt1tdm1(i,d+1)=0;
                            else
                                C_pp_ttdm1gt1tdm1(i,d+1)=copulacdf('Gaussian',[F_p_tgt1t2(i) F_p_t3gt1t2(i)],AFinal_pp);
                            end
                        end
                        
                    else
                        
                        for i=1:T-(counter+1)
                            if F_p_tgt1t2(i)==1
                                C_pp_ttdm1gt1tdm1(i,d+1)=F_p_t3gt1t2(i);
                            elseif F_p_t3gt1t2(i)==1
                                C_pp_ttdm1gt1tdm1(i,d+1)=F_p_tgt1t2(i);
                            elseif F_p_tgt1t2(i)==0 || F_p_t3gt1t2(i)==0
                                C_pp_ttdm1gt1tdm1(i,d+1)=0;
                            else
                                C_pp_ttdm1gt1tdm1(i,d+1)=F_p_tgt1t2(i)*F_p_t3gt1t2(i);
                            end
                        end
                        
                    end
                    
                    
                    
                    if d<=Tr

                    TauFinal_pn=kendalltau([S_star(counter+1:end)-1,S_star(1:end-counter)]);
                    AFinal_pn=copulaparam('Gaussian',TauFinal_pn(1,2));
                        
                        for i=1:T-(counter+1)
                            if F_p_tgt1t2(i)==1
                                C_pn_ttdm1gt1tdm1(i,d+1)=F_n_t3gt1t2(i);
                            elseif F_n_t3gt1t2(i)==1
                                C_pn_ttdm1gt1tdm1(i,d+1)=F_p_tgt1t2(i);
                            elseif F_p_tgt1t2(i)==0 || F_n_t3gt1t2(i)==0
                                C_pn_ttdm1gt1tdm1(i,d+1)=0;
                            else
                                C_pn_ttdm1gt1tdm1(i,d+1)=copulacdf('Gaussian',[F_p_tgt1t2(i) F_n_t3gt1t2(i)],AFinal_pn);
                            end
                        end
                    else
                        for i=1:T-(counter+1)
                            if F_p_tgt1t2(i)==1
                                C_pn_ttdm1gt1tdm1(i,d+1)=F_n_t3gt1t2(i);
                            elseif F_n_t3gt1t2(i)==1
                                C_pn_ttdm1gt1tdm1(i,d+1)=F_p_tgt1t2(i);
                            elseif F_p_tgt1t2(i)==0 || F_n_t3gt1t2(i)==0
                                C_pn_ttdm1gt1tdm1(i,d+1)=0;
                            else
                                C_pn_ttdm1gt1tdm1(i,d+1)=F_p_tgt1t2(i)*F_n_t3gt1t2(i);
                            end
                        end
                    end
                    
                    
                    if d<=Tr

                    TauFinal_np=kendalltau([S_star(counter+1:end),S_star(1:end-counter)-1]);
                    AFinal_np=copulaparam('Gaussian',TauFinal_np(1,2));
                        
                        for i=1:T-(counter+1)
                            if F_n_tgt1t2(i)==1
                                C_np_ttdm1gt1tdm1(i,d+1)=F_p_t3gt1t2(i);
                            elseif F_p_t3gt1t2(i)==1
                                C_np_ttdm1gt1tdm1(i,d+1)=F_n_tgt1t2(i);
                            elseif F_n_tgt1t2(i)==0 || F_p_t3gt1t2(i)==0
                                C_np_ttdm1gt1tdm1(i,d+1)=0;
                            else
                                C_np_ttdm1gt1tdm1(i,d+1)=copulacdf('Gaussian',[F_n_tgt1t2(i) F_p_t3gt1t2(i)],AFinal_np);
                            end
                        end
                        
                    else
                        
                        for i=1:T-(counter+1)
                            if F_n_tgt1t2(i)==1
                                C_np_ttdm1gt1tdm1(i,d+1)=F_p_t3gt1t2(i);
                            elseif F_p_t3gt1t2(i)==1
                                C_np_ttdm1gt1tdm1(i,d+1)=F_n_tgt1t2(i);
                            elseif F_n_tgt1t2(i)==0 || F_p_t3gt1t2(i)==0
                                C_np_ttdm1gt1tdm1(i,d+1)=0;
                            else
                                C_np_ttdm1gt1tdm1(i,d+1)=F_n_tgt1t2(i)*F_p_t3gt1t2(i);
                            end
                        end
                    end
                    
                    
                    if d<=Tr

                    TauFinal_nn=kendalltau([S_star(counter+1:end)-1,S_star(1:end-counter)-1]);
                    AFinal_nn=copulaparam('Gaussian',TauFinal_nn(1,2));
                        
                        for i=1:T-(counter+1)
                            if F_n_tgt1t2(i)==1
                                C_nn_ttdm1gt1tdm1(i,d+1)=F_n_t3gt1t2(i);
                            elseif F_n_t3gt1t2(i)==1
                                C_nn_ttdm1gt1tdm1(i,d+1)=F_n_tgt1t2(i);
                            elseif F_n_tgt1t2(i)==0 || F_n_t3gt1t2(i)==0
                                C_nn_ttdm1gt1tdm1(i,d+1)=0;
                            else
                                C_nn_ttdm1gt1tdm1(i,d+1)=copulacdf('Gaussian',[F_n_tgt1t2(i) F_n_t3gt1t2(i)],AFinal_nn);
                            end
                        end
                        
                    else
                        
                        for i=1:T-(counter+1)
                            if F_n_tgt1t2(i)==1
                                C_nn_ttdm1gt1tdm1(i,d+1)=F_n_t3gt1t2(i);
                            elseif F_n_t3gt1t2(i)==1
                                C_nn_ttdm1gt1tdm1(i,d+1)=F_n_tgt1t2(i);
                            elseif F_n_tgt1t2(i)==0 || F_n_t3gt1t2(i)==0
                                C_nn_ttdm1gt1tdm1(i,d+1)=0;
                            else
                                C_nn_ttdm1gt1tdm1(i,d+1)=F_n_tgt1t2(i)*F_n_t3gt1t2(i);
                            end
                        end
                    end
                    
                    
                    counter=counter+1;
                    
                end
                %Step 5
                F_p_1g2m=(C_pp_ttdm1gt1tdm1(1,end)-C_pn_ttdm1gt1tdm1(1,end))/f_conda(1,end);
                F_n_1g2m=(C_np_ttdm1gt1tdm1(1,end)-C_nn_ttdm1gt1tdm1(1,end))/f_conda(1,end);
                f_1g2m=F_p_1g2m-F_n_1g2m;
                
                %Step 6
                countpmf1=T-2;
                countpmf2=0;
                ppmmff=0;
                
                while countpmf2<=T-5
                    ppmmff=ppmmff+log(f_condb(T-countpmf1,end-countpmf2));
                    countpmf1=countpmf1-1;
                    countpmf2=countpmf2+1;
                end
                
                ppmmff=ppmmff+log(f_second)+log(Pdf(T-1));
                
                p_CFC=normcdf(X*Beta_CFC);
                POS_CFC=log(p_CFC./(1-p_CFC)).*s_y(2:end);
                Test_CFC=sum(POS_CFC);
                
            catch MyErr
                err_count = err_count + 1;
            end
            count = count + 1;
        end
        
        parfor k=1:B
            try
                mu_sim=0;
                sigma_sim=1;
                eps_sim=trnd(1,T,1);
                y_sim=eps_sim;
                %The signs
                s_ysim=double(y_sim>=0);
                
                p_sim=tcdf(X*Beta,1);
                
                
                CDF_Ps=binocdf(s_ysim(2:end),1,p_sim);
                CDF_Ns=binocdf(s_ysim(2:end)-1,1,p_sim);
                Pdfs=CDF_Ps-CDF_Ns;
                
                
                S_star_sim=s_ysim(2:end)+(rand(numel(s_ysim(2:end)),1)-1);
                
                
                sC_pp_tt1=zeros(T-2,1);
                sCDDF_P=CDF_Ps(1:end-1);
                sCDDF_P1=CDF_Ps(2:end);
                
                sTau_pp_tt1=kendalltau([S_star_sim(2:end),S_star_sim(1:end-1)]);
                sA_pp_tt1=copulaparam('Gaussian',sTau_pp_tt1(1,2));
                
                
                for i=1:T-1
                    if i+1<=T-1
                        if CDF_Ps(i)==1
                            sC_pp_tt1(i)=CDF_Ps(i+1);
                        elseif CDF_Ps(i+1)==1
                            sC_pp_tt1(i)=CDF_Ps(i);
                        elseif CDF_Ps(i)==0 || CDF_Ps(i+1)==0
                            sC_pp_tt1(i)=0;
                        else
                            sC_pp_tt1(i)=copulacdf('Gaussian',[CDF_Ps(i) CDF_Ps(i+1)],sA_pp_tt1);
                        end
                    end
                end
                
                
                sC_pn_tt1=zeros(T-2,1);
                
                sCDDF_P=CDF_Ps(1:end-1);
                sCDDF_N1=CDF_Ns(2:end);
                
                sTau_pn_tt1=kendalltau([S_star_sim(2:end)-1,S_star_sim(1:end-1)]);
                sA_pn_tt1=copulaparam('Gaussian',sTau_pn_tt1(1,2));
                
                
                for i=1:T-1
                    if i+1<=T-1
                        if CDF_Ps(i)==1
                            sC_pn_tt1(i)=CDF_Ns(i+1);
                        elseif CDF_Ns(i+1)==1
                            sC_pn_tt1(i)=CDF_Ps(i);
                        elseif CDF_Ps(i)==0 || CDF_Ns(i+1)==0
                            sC_pn_tt1(i)=0;
                        else
                            sC_pn_tt1(i)=copulacdf('Gaussian',[CDF_Ps(i) CDF_Ns(i+1)],sA_pn_tt1);
                        end
                    end
                end
                
                sC_np_tt1=zeros(T-2,1);
                
                sCDDF_N=CDF_Ns(1:end-1);
                sCDDF_P1=CDF_Ps(2:end);
                
                sTau_np_tt1=kendalltau([S_star_sim(2:end),S_star_sim(1:end-1)-1]);
                sA_np_tt1=copulaparam('Gaussian',sTau_np_tt1(1,2));
                
                for i=1:T-1
                    if i+1<=T-1
                        if CDF_Ns(i)==1
                            sC_np_tt1(i)=CDF_Ps(i+1);
                        elseif CDF_Ps(i+1)==1
                            sC_np_tt1(i)=CDF_Ns(i);
                        elseif CDF_Ns(i)==0 || CDF_Ps(i+1)==0
                            sC_np_tt1(i)=0;
                        else
                            sC_np_tt1(i)=copulacdf('Gaussian',[CDF_Ns(i) CDF_Ps(i+1)],sA_np_tt1);
                        end
                    end
                end
                
                sC_nn_tt1=zeros(T-2,1);
                
                sCDDF_N=CDF_Ns(1:end-1);
                sCDDF_N1=CDF_Ns(2:end);
                
                sTau_nn_tt1=kendalltau([S_star_sim(2:end)-1,S_star_sim(1:end-1)-1]);
                sA_nn_tt1=copulaparam('Gaussian',sTau_nn_tt1(1,2));
                
                
                for i=1:T-1
                    if i+1<=T-1
                        if CDF_Ns(i)==1
                            sC_nn_tt1(i)=CDF_Ns(i+1);
                        elseif CDF_Ns(i+1)==1
                            sC_nn_tt1(i)=CDF_Ns(i);
                        elseif CDF_Ns(i+1)==0 || CDF_Ns(i)==0
                            sC_nn_tt1(i)=0;
                        else
                            sC_nn_tt1(i)=copulacdf('Gaussian',[CDF_Ns(i) CDF_Ns(i+1)],sA_nn_tt1);
                        end
                    end
                end
                
                
                
                %Step3(a)
                sF_p_tgt1=(sC_pp_tt1(1:end-1)-sC_pn_tt1(1:end-1))./Pdfs(2:end-1);
                sF_n_tgt1=(sC_np_tt1(1:end-1)-sC_nn_tt1(1:end-1))./Pdfs(2:end-1);
                sf_tgt1=sF_p_tgt1-sF_n_tgt1;
                
                %Step3(a)
                sF_p_second=(sC_pp_tt1(end-1)-sC_pn_tt1(end-1))./Pdfs(end);
                sF_n_second=(sC_np_tt1(end-1)-sC_nn_tt1(end-1))./Pdfs(end);
                sf_second=sF_p_second-sF_n_second;
                
                %Step3(b)
                sF_p_t2gt1=(sC_pp_tt1(2:end)-sC_np_tt1(2:end))./Pdfs(2:end-1);
                sF_n_t2gt1=(sC_pn_tt1(2:end)-sC_nn_tt1(2:end))./Pdfs(2:end-1);
                sf_t2gt1=sF_p_t2gt1-sF_n_t2gt1;
                
                %Step3(c)
                
                sC_pp_tt2gt1=zeros(T-3,1);
                
                sTau_pp_tt2gt1=kendalltau([S_star_sim(3:end),S_star_sim(1:end-2)]);
                sA_pp_tt2gt1=copulaparam('Gaussian',sTau_pp_tt2gt1);
                
                for i=1:T-3
                    if sF_p_tgt1(i)==1
                        sC_pp_tt2gt1(i)=sF_p_t2gt1(i);
                    elseif sF_p_t2gt1(i)==1
                        sC_pp_tt2gt1(i)=sF_p_tgt1(i);
                    elseif sF_p_tgt1(i)==0 || sF_p_t2gt1(i)==0
                        sC_pp_tt2gt1(i)=0;
                    else
                        sC_pp_tt2gt1(i)=copulacdf('Gaussian',[sF_p_tgt1(i) sF_p_t2gt1(i)],sA_pp_tt2gt1);
                    end
                end
                
                
                sC_pn_tt2gt1=zeros(T-3,1);
                
                sTau_pn_tt2gt1=kendalltau([S_star_sim(3:end)-1,S_star_sim(1:end-2)]);
                sA_pn_tt2gt1=copulaparam('Gaussian',sTau_pn_tt2gt1);
                
                
                for i=1:T-3
                    if sF_p_tgt1(i)==1
                        sC_pn_tt2gt1(i)=sF_n_t2gt1(i);
                    elseif sF_n_t2gt1(i)==1
                        sC_pn_tt2gt1(i)=sF_p_tgt1(i);
                    elseif sF_p_tgt1(i)==0 || sF_n_t2gt1(i)==0
                        sC_pn_tt2gt1(i)=0;
                    else
                        sC_pn_tt2gt1(i)=copulacdf('Gaussian',[sF_p_tgt1(i) sF_n_t2gt1(i)],sA_pn_tt2gt1);
                    end
                end
                
                sC_np_tt2gt1=zeros(T-3,1);
                sTau_np_tt2gt1=kendalltau([S_star_sim(3:end),S_star_sim(1:end-2)-1]);
                sA_np_tt2gt1=copulaparam('Gaussian',sTau_np_tt2gt1);                
                
                for i=1:T-3
                    if sF_n_tgt1(i)==1
                        sC_np_tt2gt1(i)=sF_p_t2gt1(i);
                    elseif sF_p_t2gt1(i)==1
                        sC_np_tt2gt1(i)=sF_n_tgt1(i);
                    elseif sF_n_tgt1(i)==0 || sF_p_t2gt1(i)==0
                        sC_np_tt2gt1(i)=0;
                    else
                        sC_np_tt2gt1(i)=copulacdf('Gaussian',[sF_n_tgt1(i) sF_p_t2gt1(i)],sA_np_tt2gt1);
                    end
                end
                
                sC_nn_tt2gt1=zeros(T-3,1);
                sTau_nn_tt2gt1=kendalltau([S_star_sim(3:end)-1,S_star_sim(1:end-2)-1]);
                sA_nn_tt2gt1=copulaparam('Gaussian',sTau_nn_tt2gt1);                
                
                
                for i=1:T-3
                    if sF_n_tgt1(i)==1
                        sC_nn_tt2gt1(i)=sF_n_t2gt1(i);
                    elseif sF_n_t2gt1(i)==1
                        sC_nn_tt2gt1(i)=sF_n_tgt1(i);
                    elseif sF_n_tgt1(i)==0 || sF_n_t2gt1(i)==0
                        sC_nn_tt2gt1(i)=0;
                    else
                        sC_nn_tt2gt1(i)=copulacdf('Gaussian',[sF_n_tgt1(i) sF_n_t2gt1(i)],sA_nn_tt2gt1);
                    end
                end
                
                sC_pp_ttdm1gt1tdm1=zeros(T-3,T-4);
                sC_pp_ttdm1gt1tdm1(1:T-3,1)=sC_pp_tt2gt1;
                
                
                sC_pn_ttdm1gt1tdm1=zeros(T-3,T-4);
                sC_pn_ttdm1gt1tdm1(1:T-3,1)=sC_pn_tt2gt1;
                
                
                sC_np_ttdm1gt1tdm1=zeros(T-3,T-4);
                sC_np_ttdm1gt1tdm1(1:T-3,1)=sC_np_tt2gt1;
                
                
                sC_nn_ttdm1gt1tdm1=zeros(T-3,T-4);
                sC_nn_ttdm1gt1tdm1(1:T-3,1)=sC_nn_tt2gt1;
                
                sf_conda=zeros(T-3,T-4);
                sf_conda(1:T-3,1)=sf_t2gt1;
                
                sf_condb=zeros(T-3,T-4);
                sf_condb(1:T-3,1)=sf_tgt1;
                
                %      clear sC_pp_tt2gt1 sC_pn_tt2gt1 sC_np_tt2gt1 sC_nn_tt2gt1 sf_t2gt1
                
                counter=3;
                
                for d=1:T-4
                    
                    
                    %Step 4(a)
                    sF_p_tgt1t2=(sC_pp_ttdm1gt1tdm1(1:T-counter,d)-sC_pn_ttdm1gt1tdm1(1:T-counter,d))./sf_conda(1:T-counter,d);
                    sF_n_tgt1t2=(sC_np_ttdm1gt1tdm1(1:T-counter,d)-sC_nn_ttdm1gt1tdm1(1:T-counter,d))./sf_conda(1:T-counter,d);
                    
                    sF_p_tgt1t2(isnan(sF_p_tgt1t2))=0;
                    sF_n_tgt1t2(isnan(sF_n_tgt1t2))=0;
                    
                    %Step 4(b)
                    sF_p_t3gt1t2=(sC_pp_ttdm1gt1tdm1(2:T-counter,d)-sC_np_ttdm1gt1tdm1(2:T-counter,d))./sf_condb(2:T-counter,d);
                    sF_n_t3gt1t2=(sC_pn_ttdm1gt1tdm1(2:T-counter,d)-sC_nn_ttdm1gt1tdm1(2:T-counter,d))./sf_condb(2:T-counter,d);
                    
                    
                    sF_p_t3gt1t2(isnan(sF_p_t3gt1t2))=0;
                    sF_n_t3gt1t2(isnan(sF_n_t3gt1t2))=0;
                    
                    
                    sf_conda(1:T-(counter+1),d+1)=sF_p_t3gt1t2-sF_n_t3gt1t2;
                    sf_condb(1:T-counter,d+1)=sF_p_tgt1t2-sF_n_tgt1t2;
                    
                    %Step 4(c)
                    
                    if d<=Tr
                    sTauFinal_pp=kendalltau([S_star_sim(counter+1:end),S_star_sim(1:end-counter)]);
                    sAFinal_pp=copulaparam('Gaussian',sTauFinal_pp(1,2));
                        
                        for i=1:T-(counter+1)
                            if sF_p_tgt1t2(i)==1
                                sC_pp_ttdm1gt1tdm1(i,d+1)=sF_p_t3gt1t2(i);
                            elseif sF_p_t3gt1t2(i)==1
                                sC_pp_ttdm1gt1tdm1(i,d+1)=sF_p_tgt1t2(i);
                            elseif sF_p_tgt1t2(i)==0 || sF_p_t3gt1t2(i)==0
                                sC_pp_ttdm1gt1tdm1(i,d+1)=0;
                            else
                                sC_pp_ttdm1gt1tdm1(i,d+1)=copulacdf('Gaussian',[sF_p_tgt1t2(i) sF_p_t3gt1t2(i)],sAFinal_pp);
                            end
                        end
                        
                    else
                        
                        for i=1:T-(counter+1)
                            if sF_p_tgt1t2(i)==1
                                sC_pp_ttdm1gt1tdm1(i,d+1)=sF_p_t3gt1t2(i);
                            elseif sF_p_t3gt1t2(i)==1
                                sC_pp_ttdm1gt1tdm1(i,d+1)=sF_p_tgt1t2(i);
                            elseif sF_p_tgt1t2(i)==0 || sF_p_t3gt1t2(i)==0
                                sC_pp_ttdm1gt1tdm1(i,d+1)=0;
                            else
                                sC_pp_ttdm1gt1tdm1(i,d+1)=sF_p_tgt1t2(i)*sF_p_t3gt1t2(i);
                            end
                        end
                    end
                    
                    
                    if d<=Tr
                        
                    sTauFinal_pn=kendalltau([S_star_sim(counter+1:end)-1,S_star_sim(1:end-counter)]);
                    sAFinal_pn=copulaparam('Gaussian',sTauFinal_pn(1,2));
                        
                        for i=1:T-(counter+1)
                            if sF_p_tgt1t2(i)==1
                                sC_pn_ttdm1gt1tdm1(i,d+1)=sF_n_t3gt1t2(i);
                            elseif sF_n_t3gt1t2(i)==1
                                sC_pn_ttdm1gt1tdm1(i,d+1)=sF_p_tgt1t2(i);
                            elseif sF_p_tgt1t2(i)==0 || sF_n_t3gt1t2(i)==0
                                sC_pn_ttdm1gt1tdm1(i,d+1)=0;
                            else
                                sC_pn_ttdm1gt1tdm1(i,d+1)=copulacdf('Gaussian',[sF_p_tgt1t2(i) sF_n_t3gt1t2(i)],sAFinal_pn);
                            end
                        end
                        
                    else
                        
                        for i=1:T-(counter+1)
                            if sF_p_tgt1t2(i)==1
                                sC_pn_ttdm1gt1tdm1(i,d+1)=sF_n_t3gt1t2(i);
                            elseif sF_n_t3gt1t2(i)==1
                                sC_pn_ttdm1gt1tdm1(i,d+1)=sF_p_tgt1t2(i);
                            elseif sF_p_tgt1t2(i)==0 || sF_n_t3gt1t2(i)==0
                                sC_pn_ttdm1gt1tdm1(i,d+1)=0;
                            else
                                sC_pn_ttdm1gt1tdm1(i,d+1)=sF_p_tgt1t2(i)*sF_n_t3gt1t2(i);
                            end
                        end
                    end
                    
                    
                    if d<=Tr
                        
                    sTauFinal_np=kendalltau([S_star_sim(counter+1:end),S_star_sim(1:end-counter)-1]);
                    sAFinal_np=copulaparam('Gaussian',sTauFinal_np(1,2));
                        
                        for i=1:T-(counter+1)
                            if sF_n_tgt1t2(i)==1
                                sC_np_ttdm1gt1tdm1(i,d+1)=sF_p_t3gt1t2(i);
                            elseif sF_p_t3gt1t2(i)==1
                                sC_np_ttdm1gt1tdm1(i,d+1)=sF_n_tgt1t2(i);
                            elseif sF_n_tgt1t2(i)==0 || sF_p_t3gt1t2(i)==0
                                sC_np_ttdm1gt1tdm1(i,d+1)=0;
                            else
                                sC_np_ttdm1gt1tdm1(i,d+1)=copulacdf('Gaussian',[sF_n_tgt1t2(i) sF_p_t3gt1t2(i)],sAFinal_np);
                            end
                        end
                        
                    else
                        
                        for i=1:T-(counter+1)
                            if sF_n_tgt1t2(i)==1
                                sC_np_ttdm1gt1tdm1(i,d+1)=sF_p_t3gt1t2(i);
                            elseif sF_p_t3gt1t2(i)==1
                                sC_np_ttdm1gt1tdm1(i,d+1)=sF_n_tgt1t2(i);
                            elseif sF_n_tgt1t2(i)==0 || sF_p_t3gt1t2(i)==0
                                sC_np_ttdm1gt1tdm1(i,d+1)=0;
                            else
                                sC_np_ttdm1gt1tdm1(i,d+1)=sF_n_tgt1t2(i)*sF_p_t3gt1t2(i);
                            end
                        end
                    end
                    
                    if d<=Tr
                        
                    sTauFinal_nn=kendalltau([S_star_sim(counter+1:end)-1,S_star_sim(1:end-counter)-1]);
                    sAFinal_nn=copulaparam('Gaussian',sTauFinal_nn(1,2));
                        
                        for i=1:T-(counter+1)
                            if sF_n_tgt1t2(i)==1
                                sC_nn_ttdm1gt1tdm1(i,d+1)=sF_n_t3gt1t2(i);
                            elseif sF_n_t3gt1t2(i)==1
                                sC_nn_ttdm1gt1tdm1(i,d+1)=sF_n_tgt1t2(i);
                            elseif sF_n_tgt1t2(i)==0 || sF_n_t3gt1t2(i)==0
                                sC_nn_ttdm1gt1tdm1(i,d+1)=0;
                            else
                                sC_nn_ttdm1gt1tdm1(i,d+1)=copulacdf('Gaussian',[sF_n_tgt1t2(i) sF_n_t3gt1t2(i)],sAFinal_nn);
                            end
                        end
                        
                    else
                        
                        for i=1:T-(counter+1)
                            if sF_n_tgt1t2(i)==1
                                sC_nn_ttdm1gt1tdm1(i,d+1)=sF_n_t3gt1t2(i);
                            elseif sF_n_t3gt1t2(i)==1
                                sC_nn_ttdm1gt1tdm1(i,d+1)=sF_n_tgt1t2(i);
                            elseif sF_n_tgt1t2(i)==0 || sF_n_t3gt1t2(i)==0
                                sC_nn_ttdm1gt1tdm1(i,d+1)=0;
                            else
                                sC_nn_ttdm1gt1tdm1(i,d+1)=sF_n_tgt1t2(i)*sF_n_t3gt1t2(i);
                            end
                        end
                    end
                    
                    
                    counter=counter+1;
                    
                end
                %Step 5
                sF_p_1g2m=(sC_pp_ttdm1gt1tdm1(1,end)-sC_pn_ttdm1gt1tdm1(1,end))/sf_conda(1,end);
                sF_n_1g2m=(sC_np_ttdm1gt1tdm1(1,end)-sC_nn_ttdm1gt1tdm1(1,end))/sf_conda(1,end);
                sf_1g2m=sF_p_1g2m-sF_n_1g2m;
                
                %Step 6
                countpmf1=T-2;
                countpmf2=0;
                sppmmff=0;
                
                while countpmf2<=T-5
                    sppmmff=sppmmff+log(sf_condb(T-countpmf1,end-countpmf2));
                    countpmf1=countpmf1-1;
                    countpmf2=countpmf2+1;
                end
                
                sppmmff=sppmmff+log(sf_second)+log(Pdfs(T-1));
                
                
                p_CFCsim=normcdf(X*Beta_CFC);
                POS_CFCsim=log(p_CFCsim./(1-p_CFCsim)).*s_ysim(2:end);
                Test_CFCsim=sum(POS_CFCsim);
                
                MCDist(k)=sppmmff;
                MCDist_Old(k)=Test_CFCsim;
                
            end
        end
        
        
        UB=quantile(MCDist,0.95);
        %LB=quantile(MCDist,0.025);
        
        if ppmmff>UB %|| ppmmff<LB
            
            counterPower=counterPower+1;
            
        end
        
        UB_CFC=quantile(MCDist_Old,0.95);
        %LB=quantile(MCDist,0.025);
        
        if Test_CFC>UB_CFC %|| ppmmff<LB
            
            counterPower_POS=counterPower_POS+1;
            
        end
        
        %Other Tests->
        %NonParametric Test of CD
        CD=y(2:end).*lagX;
        S_CD=double(CD>=0);
        Test_CD=sum(S_CD);
        
        if Test_CD>30
            
            counterPower_CD=counterPower_CD+1;
            
        end
        
        
        % Estimation of beta by OLS */
        betahat=((lagX')*lagX)\(lagX')*y(2:end);
        
                % Estimation of variance of betahat */
        uhat=y(2:end)-betahat*lagX;
        seghat=((uhat')*uhat)/(T-2);
        segbhat=seghat/(((lagX')*lagX));
        
        
        % Calculation of robust standard deviation of White */
        moy=(1/(T-1))*sum(lagX);
        SST=sum((lagX-moy.*ones(T-1,1)).^2);
        var3=(((lagX-moy.*ones(T-1,1)).^2)')*  (uhat.^2);
        Varob=var3/(SST^2);
 
        
        
        if betahat/sqrt(segbhat)> 1.675
            counterPower_T=counterPower_T+1;
        end
        
        % Decision rule of testing H0: beta=0  using corrected white T-test*/
        
        if betahat/sqrt(Varob)> 1.675
            counterPower_WT=counterPower_WT+1;
        end
        
        
        
        clearvars -except counterPower counterPower_POS counterPower_CD counterPower_WT counterPower_T Iter Tr Round CurveCounter Coef Power
        
        Round=Round+1;
        
        disp(['Round ',num2str(Round)])
        
        
        Power(CurveCounter,1)=(counterPower/Iter)*100;
        Power(CurveCounter,2)=(counterPower_CD/Iter)*100;
        Power(CurveCounter,3)=(counterPower_T/Iter)*100;
        Power(CurveCounter,4)=(counterPower_WT/Iter)*100;
        Power(CurveCounter,5)=(counterPower_POS/Iter)*100;

        
        
    end
    
    
    CurveCounter=CurveCounter+1;
    Percento=(((Coef/0.05)+1)/((0.5/0.05)+1))*100;
    disp([num2str(Percento), '% is complete! Next coefficient is ' num2str(Coef+0.05)])
    Coef=Coef+0.05;

    
end

p=toc;
disp(['Code took ',num2str(p/60), ' minutes to run!'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if results are stored in MATLAB variable(s) e.g. 'OUT' these could be saved to ASCII text file(s), with
% base name e.g. output and labeled by the data set index. The file name is stored in the variable output_file as

output_file = strcat( 'Powert1.txt' );

% and the value of OUT saved

save( output_file , 'Power' ,  '-ascii' );

%
