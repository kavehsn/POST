
% this is a wrapper for a .m file to use with SLURM on Hamilton
% read the index of the data set to be processed from standard input and store in the variable data_set_index

data_set_index = input('');

% print it as a reminder (can be commented with %)

fprintf('\nIndex of data set is:  %d \n', data_set_index);

% compose the name of the inputfile(s) from a base name ,e.g. 'input' 
% and the index

input_file = strcat('GridSearch_', int2str(data_set_index));

% and print it

fprintf( 'Input file for data set is: %s \n' , input_file );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here the code or call of a function should be included to import the data set from the file input_file 

tic
%Extract the Null Matrix

SS=load('GridSearch.mat',input_file);

C = fieldnames(SS);

for k = 1:numel(C)
     GridSearch=SS.(C{k});
end

dims=size(GridSearch);

NullRetain=zeros(dims(1),dims(2));

load ('ProjectionTechnique.mat')

B=1000; %Number of Monte Carlo Iterations
MCDist=zeros(B,1);


%Split the sample
x_Est1=x_1(1:round(0.1*(length(x_1))));
x_Test1=x_1(round(0.1*(length(x_1)))+1:end);

x_Est2=x_2(1:round(0.1*(length(x_2))));
x_Test2=x_2(round(0.1*(length(x_2)))+1:end);

y_Est=y(1:round(0.1*(length(y))));
y_Test=y(round(0.1*(length(y)))+1:end);

lagX_Est1=x_Est1(1:end-1);
lagX_Test1=x_Test1(1:end-1);
lagX_1=x_1(1:end-1);

lagX_Est2=x_Est2(1:end-1);
lagX_Test2=x_Test2(1:end-1);
lagX_2=x_2(1:end-1);

X_Est=[lagX_Est1 lagX_Est2];
X_Test=[ones(length(lagX_Test1),1) lagX_Test1 lagX_Test2];
X=[ones(length(lagX_1),1) lagX_1 lagX_2];


% Estimation of beta by OLS */
betahat=robustfit(X_Est,y_Est(2:end));

counter=1;


for i=1:length(GridSearch)
    
    Beta_Trial=GridSearch(i,:);
    Beta_Trial=Beta_Trial';

                
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
		    sigma_alt_sim=100000;
		    sgn_sim=binornd(1,0.05,0.9*T,1);
                    eps_sim=sgn_sim.*normrnd(mu_sim,sigma_alt_sim,0.9*T,1)+(1-sgn_sim).*normrnd(mu_sim,sigma_sim,0.9*T,1);                    
                    
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
                
                UB=quantile(MCDist,0.975);
                LB=quantile(MCDist,0.025);
                
                if Test<UB && Test>LB
                    
                    NullRetain(counter,:)=Beta_Trial;
                    counter=counter+1;
                    
                end
                
                
                if mod(i,100)==0
                    
                    disp(['Elapased ',num2str(i),'! Percentage= ',num2str(i/length(GridSearch))])
                
                end
end
pp=toc;

NullRetain=NullRetain(1:counter,:);

disp(['Code took ',num2str(pp/60),' minutes to complete!'])
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if results are stored in MATLAB variable(s) e.g. 'OUT' these could be saved to ASCII text file(s), with
% base name e.g. output and labeled by the data set index. The file name is stored in the variable output_file as

output_file = strcat( 'retainNull' , int2str( data_set_index ) , '.txt' );

% and the value of OUT saved

save( output_file , 'NullRetain' ,  '-ascii' );
                
    