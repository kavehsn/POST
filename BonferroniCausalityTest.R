

library(MASS)
library(pracma)
library(CompQuadForm)
library(robustbase)


# Start the clock!
ptm <- proc.time()
    #Assign the values
counter_SampleSize<-50

RatioR<-matrix(data=NA,nrow=2,ncol=1)
RatioA<-matrix(data=NA,nrow=2,ncol=1)
RatioI<-matrix(data=NA,nrow=2,ncol=1)
Coef<-matrix(data=NA,nrow=2,ncol=1)


while(counter_SampleSize<=100){

    T<-counter_SampleSize
    B<-500
    Sims<-250
    Reject<-0
    Accept<-0
    Inconclusive<-0


    Alpha<-0.05
    Alpha1<-0.008
    Alpha2<-Alpha-Alpha1
    Alpha3<-Alpha+Alpha1


    for(SimRun in 1:Sims){


        Test_Sim<-matrix(data=NA,nrow=B,ncol=1)



        Intercept_TrueDGP<-matrix(
            c(0,0),
            nrow=2,
            ncol=1,
            byrow=TRUE)
        Beta_Gen_TrueDGP<-matrix(
            c(0,0,0,0),
            nrow=2,
            ncol=2,
            byrow=TRUE)

        mu_e<-matrix(
            c(0,0),
            nrow=2,
            ncol=1,
            byrow=TRUE)

        sigma<-matrix(
            c(1,0,0,1),
            nrow=2,
            ncol=2,
            byrow=TRUE)

        eps<-mvrnorm(T,mu_e,sigma,empirical=FALSE)

    #Generate the data

        z <- matrix(data=NA,nrow=T,ncol=2)
        z[1,]<-eps[1,]
        for(i in 2:T){
            z[i,]=Intercept_TrueDGP+Beta_Gen_TrueDGP%*%z[i-1,]+eps[i,]
        }

    #Separate the arguments

        y<-matrix(
            z[,1],
            nrow=T,
            ncol=1,
            byrow=TRUE)

        x<-matrix(
            z[,2],
            nrow=T,
            ncol=1,
            byrow=TRUE)

    #Transform X
    #x<-(x-median(x))    


    #LS Estimator for LR Test

        lagY<-c(NA, y[seq_along(y) -1])

        lagX<-c(NA, x[seq_along(x) -1])

    ##OLSBeta<-ar.ols(y, aic=FALSE,order.max=1)
    ##Beta_LS<-c(OLSBeta$x.intercept,OLSBeta$ar)

    fit_Alt<-ltsReg(y ~ lagY + lagX)
    Coefs_Alt<-coefficients(fit_Alt)

    #Constant_UR<-Coefs_Alt[1]
    #CoefY_UR<-Coefs_Alt[2]
    #CoefX_UR<-Coefs_Alt[3]

    #y_Trans_Alt<-y-Constant_UR_Alt-(CoefY_UR_Alt*lagY)



    #fit_AltMain<-ltsReg(y_Trans_Alt ~ lagY + lagX)
    #Coefs_AltMain<-coefficients(fit_AltMain)

    #Constant_UR<-Coefs_AltMain[1]
    #CoefY_UR<-Coefs_AltMain[2]
    #CoefX_UR<-Coefs_AltMain[3]

    #LS Estimator for Imhof calculation

        I<-diag(T-1)
        zeros<-matrix(data=0,nrow=T-1,ncol=1)

        D_0<-cbind(zeros,I)
        D_T<-cbind(I,zeros)

        X_1<-matrix(data=1,nrow=T-1,ncol=1)
        x_2<-matrix(data=0, nrow=T-1,ncol=1)
        counter<-0

        for(s in 1:T-1){

            x_2[s,1]<-counter+1
            counter<-counter+1
        }

        X_2<-cbind(X_1,x_2)

        P_2<-X_2%*%solve(t(X_2)%*%X_2)%*%t(X_2)

        Alpha_Est<-t(y)%*%t(D_0)%*%(diag(T-1)-P_2)%*%D_T%*%y%*%solve((t(y)%*%t(D_T)%*%(diag(T-1)-P_2)%*%D_T%*%y))
        Beta_Est<-t(x)%*%t(D_0)%*%(diag(T-1)-P_2)%*%D_T%*%x%*%solve((t(x)%*%t(D_T)%*%(diag(T-1)-P_2)%*%D_T%*%x))


        y_hat<-y-(Alpha_Est*lagY)
        x_hat<-x-(Beta_Est*lagX)


        y_hat <- y_hat[!is.na(y_hat)]
        x_hat <- x_hat[!is.na(x_hat)]

        Int_Est<-mean(y_hat)
        Int_Est_x<-mean(x_hat)


    #Calculate the matrix components

        b<-1/(1-Alpha_Est^2)^0.5

        counter1<-0

        R_a<-diag(T)

        for(ss in 1:T){

            if (ss==1){

                for (sss in 1:T){

                    R_a[sss,ss]<-b*(Alpha_Est^counter1)
                    counter1<-counter1+1;

                }

                counter1<-1}
                else if(ss<T){

                    counter1<-1

                    for(sss in (ss+1):length(y)){

                        R_a[sss,ss]<-Alpha_Est^counter1
                        counter1<-counter1+1

                    }


                }

            }
    #Calculate the residuals

            U<-solve(R_a)%*%y    

            rrr<-seq(-1,1,0.01)


            for(j in 1:length(rrr)){

    #Calculate the weight matrix
                W<-t(R_a)%*%((t(D_0)%*%(diag(T-1)-P_2)%*%D_T)/2+(t(D_T)%*%(diag(T-1)-P_2)%*%D_0)/2-rrr[j]*t(D_T)%*%(diag(T-1)-P_2)%*%D_T)%*%R_a


    #Calculate characteristic roots
                CP<-charpoly(W, info=TRUE)
                lambda<-roots(CP$cp)

    #Imhof Algorithm
                PV<-imhof(0, lambda, h = rep(1, length(lambda)),
                    delta = rep(0, length(lambda)),
                    epsabs = 10^(-6), epsrel = 10^(-6), limit = 10000)
                PVal<-PV$Qq
                ActPVal<-1-PVal

                if(round(ActPVal, digits=3)==(1-(Alpha1/2))|round(ActPVal, digits=3)==(1-(Alpha1/2))-10^(-3)|round(ActPVal, digits=3)==(1-(Alpha1/2))+10^(-3)){
                    UpperBound<-rrr[j]
                }else if(round(ActPVal, digits=3)==(Alpha1/2)|round(ActPVal, digits=3)==(Alpha1/2)-10^(-3)|round(ActPVal, digits=3)==(Alpha1/2)+10^(-3)){
                    LowerBound<-rrr[j]
                }

            }

    #Run iteratios for each value in the interval

            CI<-seq(LowerBound,UpperBound,0.01)

            Intercept <- matrix(data=NA,nrow=length(CI),ncol=1)

            Test_Actual<-matrix(data=NA,nrow=length(CI),ncol=1)

        #SCDFU<-matrix(data=NA,nrow=T,ncol=length(CI))

            counterS<-1

            for(ii in 1:length(CI)){

    #Transformation

                y_hatt<-y-(CI[ii]*lagY)

                y_hatt <- y_hatt[!is.na(y_hatt)]

                Intercept[ii]<-mean(y_hatt)

                y_Trans<-y-Intercept[ii]-(CI[ii]*lagY)


                s<-as.double(y_Trans>=0)


                SCDFU[,counterS]<-pnorm((Intercept[ii]-Coefs_Alt[1])+(CI[ii]-Coefs_Alt[2])*lagY-(Coefs_Alt[3]*lagX))
                #SCDFU<-pnorm((CoefX_UR*lagX))


                counterS<-counterS+1


                BetaGamma<-s*(log(SCDFU[,ii]/(1-SCDFU[,ii])))
                #BetaGamma<-s*(log(SCDFU/(1-SCDFU)))
                BeforeMean<- BetaGamma[2:length(BetaGamma)]
                BeforeMean <- BeforeMean[!is.na(BeforeMean)]
                Test_Actual[ii]<-sum(BeforeMean) 

            }



            Infim<-min(Test_Actual)
            Suprem<-max(Test_Actual)


            #InfLoc<-which.min(Test_Actual)
            #SupLoc<-which.max(Test_Actual)



        #Simulate distribution under the null hypothesis (Known)

            for(t in 1:B){


               y_Trans_s<-rnorm(T,0,1)



               ss<-as.double(y_Trans_s>=0)


           #SCDFU_S<-sample(SCDFU,T,replace=TRUE)
           #SCDFU_S<-pnorm((Intercept[MedLoc]-Intercept_TrueDGP[1,1])+(CI[MedLoc]-Beta_Gen_TrueDGP[1,1])*lagY+(Beta_Gen_TrueDGP[1,2]*lagX))
               SCDFU_S<-pnorm((Int_Est-Coefs_Alt[1])+(Alpha_Est-Coefs_Alt[2])*lagY-(Coefs_Alt[3]*lagX))

               BetaGamma_s<-ss*(log(SCDFU_S/(1-SCDFU_S)))
               BeforeMean_s<- BetaGamma_s[2:length(BetaGamma_s)]
               BeforeMean_s<- BeforeMean_s[!is.na(BeforeMean_s)]
               Test_Sim[t]<-sum(BeforeMean_s)


           }


           TestLib<-quantile(Test_Sim,(Alpha3))
           TestCons<-quantile(Test_Sim,(Alpha2))

           if(Suprem<TestCons){
            Reject<-Reject+1
        } else if(Infim>=TestLib){
            Accept<-Accept+1
        }else{
            Inconclusive<-Inconclusive+1
        }


    }

    Coef[(counter_SampleSize/50)]<-counter_SampleSize
    RatioR[(counter_SampleSize/50)]<-Reject/Sims
    RatioA[(counter_SampleSize/50)]<-Accept/Sims
    RatioI[(counter_SampleSize/50)]<-Inconclusive/Sims 

    counter_SampleSize<-counter_SampleSize+50

}


RatioAll<-cbind(Coef,RatioR,RatioA,RatioI)
df<-data.frame(RatioAll)
Col_Headings<-c('Coef','RatioR','RatioA','RatioI')
names(df)<-Col_Headings


write.table(df,"DGP1YX.txt",sep="\t",row.names=FALSE)

    # Stop the clock
proc.time() - ptm
