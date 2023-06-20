
#' Monte Carlo estimation of mean and HCR based on SAE GAMLSS
#' @description A function to estimate the mean or/and the HCR     based on unit-level small area estimation based on generalized additive models for location, scale and shape
#'
#' @param sample A dataset of sampled units. With dependent variable (named y), covariates and small area (named sa)
#' @param nonsample  A dataset of non-sampled units. With covariates and small area (named sa)
#' @param D Number of Areas
#' @param Ni 1xD vector containing the number of units (in population) for each area
#' @param ni 1xD vector containing the number of sampled-units for each area
#' @param f1 A formula object for fitting a model to the mu parameter, e.g.f1=y~x+random(sa)
#' @param f2 A formula object for fitting a model to the sigma parameter, e.g.f2=y~x+random(sa)
#' @param f3 A formula object for fitting a model to the nu parameter, e.g.f3=y~x+random(sa)
#' @param f4 A formula object for fitting a model to the tau parameter, e.g.f4=y~x+random(sa)
#' @param fdis The distribution family of the GAMLSS object
#' @param nRS Number of loop to do with the RS() algorithm
#' @param nCG Number of loop to do with the CG() algorithm
#' @param R 	Number of Monte-Carlo repetition. Default is 50
#' @param Dis Type of distribution in form of rDis where Dis is one of the distribution allowed by GAMLSS
#' @param np Number of parameters of the distribution. i.e. for the normal distribution np=2, for the GB2 distribution np=4
#' @param param The parameter to estimate, "mean" or "HCR". Default is both.
#' @param seed The seed
#' @param tau.fix A value to be fixed to 1 to obtain the Singh-Maddala distribution. Dist must be GB2 and np=4
#' @param nu.fix 	A value to be fixed to 1 to obtain the Dagum distribution. Dist must be GB2 and np=4
#' @param z 	The Poverty line. Default is equal to 0.6 of the dependent variable
#' @return a list of two list. The first list, named est has:
#' @return ME 1xD vector with the estimate of the mean and
#' @return HCR 1xD vector with the estimate of the HCR.
#' @return The second, named input_var has:
#' @return fit a GAMLSS object with all the components of the regression
#' @return f1 a formula object used for fitting a model to the mu parameter
#' @return f2 a formula object used for fitting a model to the sigma parameter
#' @return f3 a  formula object used for fitting a model to the nu parameter
#' @return f4 a formula object used for fitting a model to the tau parameter
#' @return np the number of the parameters of the choose distribution
#' @return nRS the number of loop used within the RS() algorithm. Default is 150.
#' @return nCG the number of loop used within the CG() algorithm. Default is 150.
#' @return R the number of loop used to obtain the monte-carlo estimation
#' @return D the number of area
#' @return Ni 1xD vector containing the number of units (in population) for each area
#' @return ni 1xD vector containing the number of sampled-units for each area
#' @return originaldata The dataset of sampled-units
#' @return fids The distribution used for the estimation
#' @return z The poverty line
#' @export
#' @import gamlss
#' @import dplyr
#' @import splitstackshape
#' @examples #Generate data
#' @examples #
#' @examples dep.y=data_gen(Ni=rep(10,4), D=4, M=2 ,ty ="no", k=4, b1=100,
#' @examples                x1=rnorm(40, 0,1), b2 = NULL, x2 = NULL, b3 = NULL,
#' @examples                b4 = NULL, x4 = NULL, xh = NULL, Dis=NO,
#' @examples                l=c(identity), sigma = 6, sigmah = NULL,
#' @examples                sigmae = 22, costh = NULL, seed = 1234)
#'
#' @examples data=dep.y[[1]]
#' @examples #
#' @examples #sample data with a sample fraction of 0.1
#' @examples #
#' @examples library(splitstackshape)
#' @examples #sample data
#' @examples #
#' @examples sample=stratified(data, "sa", size=0.1)
#' @examples #nonsample data
#' @examples #
#' @examples nonsample=nonsample(data=data, sample=sample, id=id)
#' @examples #estimate
#' @examples est_saegamlss(sample=sample, nonsample=nonsample,
#' @examples               D=4, Ni=rep(10,4), ni=rep(1,4),
#' @examples               f1=y~x1+random(sa), f2=NULL, f3=NULL,
#' @examples               f4=NULL, fdis=NO, nRS=150, nCG=150, R=200,
#' @examples               Dis=rNO, np=2, param=NULL,
#' @examples               seed=1234, tau.fix=NULL, nu.fix=NULL)
#' @references Mori, L., & Ferrante, M. R. (2023). Small area estimation under unit-level generalized additive models for location, scale and shape. arXiv e-prints, arXiv-2302.
#' @references Graf, M., Marin, J. M., & Molina, I. (2019). A generalized mixed model for skewed distributions applied to small area estimation. Test, 28(2), 565–597.
#' @author Lorenzo Mori and Maria Rosaria Ferrante
#' @note With "object"$input_var$fit is possible to use all the classical function used by gamlss
#' @note The Small Area have to be denoted with a number from 1 to D

est_saegamlss=function(sample, nonsample, D, Ni, ni, f1, f2=NULL, f3=NULL, f4=NULL, fdis, nRS=NULL, nCG=NULL, R=NULL,
                      Dis, np, param=NULL, seed=NULL, tau.fix=NULL, nu.fix=NULL, z=NULL){

  sa <- y <- NULL
  mixed <- NULL
  dpvar=dplyr::bind_rows(sample, nonsample)$y

  if (is.null(seed)==TRUE){seed=123}
  if (is.null(param)==TRUE){param=c("both")}

  if (is.null(R)==TRUE){R=50}
  if (is.null(nRS)==TRUE){R=150}
  if (is.null(nCG)==TRUE){R=150}

  if (is.null(f2)==TRUE){f2=y ~1}
  if (is.null(f3)==TRUE){f3=y ~1}
  if (is.null(f4)==TRUE){f4=y ~1}
  if (is.null(z)==TRUE){z=0.6*median(dpvar, na.rm=TRUE)}

  set.seed=seed;

  if(np==1){
    gam1=gamlss::gamlss(f1,
               trace=F,family=substitute(fdis), data=sample, method=mixed(substitute(nRS), substitute(nRG)))

    yns=data.frame("ys"=as.vector(gamlss::predictAll(gam1,data=sample, newdata= nonsample)$mu),
                   "sa"=nonsample$sa)
    ME=rep(0,D)
    HCR=rep(0,D)
    for (x in 1:R) {

      for (i  in 1:D) {
        yns1=subset(yns, sa==i)
        samp=as.vector(Dis(Ni[i]-ni[i], mu=yns1$ys))
        sp=dplyr::pull(subset(sample, sa==i), f1[[2]])
        samp1=c(samp,sp )
     if(param=="both" | param=="mean") {ME[i]=ME[i]+mean(samp1)}
     if(param=="both" | param=="HCR")  {HCR[i]=HCR[i]+povinc(samp1, z=z)}

      }
      message('Processing estimation loop ', x, ' of ', R)

    }
    ME=ME/x
    HCR=HCR/x
  } else if (np==2){

    gam1=gamlss::gamlss(f1, sigma.fo=f2, nu.fo=f3, tau.fo=f4,
               trace=F,family=substitute(fdis), data=sample, method=mixed(substitute(nRS), substitute(nRG)))

    yns=data.frame("ys"=as.vector(gamlss::predictAll(gam1,data=sample, newdata= nonsample)$mu),
                   "yss"=as.vector(gamlss::predictAll(gam1,data=sample, newdata=nonsample)$sigma),
                   "sa"=nonsample$sa)
    ME=rep(0,D)
    HCR=rep(0,D)
    for (x in 1:R) {

      for (i  in 1:D) {
        yns1=subset(yns, sa==i)
        samp=as.vector(Dis(Ni[i]-ni[i], mu=yns1$ys, sigma=yns1$yss))
        sp=dplyr::pull(subset(sample, sa==i), f1[[2]])
        samp1=c(samp,sp )
        if(param=="both" | param=="mean") {ME[i]=ME[i]+mean(samp1)}
        if(param=="both" | param=="HCR")  {HCR[i]=HCR[i]+povinc(samp1, z=z)}

      }
      message('Processing estimation loop ', x, ' of ', R)

    }
    ME=ME/x
    HCR=HCR/x


  } else if (np==3){

    gam1=gamlss::gamlss(f1, sigma.fo=f2, nu.fo=f3, tau.fo=f4,
               trace=F,family=substitute(fdis), data=sample, method=mixed(substitute(nRS), substitute(nRG)))

    yns=data.frame("ys"=as.vector(gamlss::predictAll(gam1,data=sample, newdata= nonsample)$mu),
                   "yss"=as.vector(gamlss::predictAll(gam1,data=sample, newdata=nonsample)$sigma),
                   "ysn"=as.vector(gamlss::predictAll(gam1,data=sample, newdata=nonsample)$nu),
                   "sa"=nonsample$sa)
    ME=rep(0,D)
    HCR=rep(0,D)
    for (x in 1:R) {

      for (i  in 1:D) {
        yns1=subset(yns, sa==i)
        samp=as.vector(Dis(Ni[i]-ni[i], mu=yns1$ys, sigma=yns1$yss, nu=yns1$ysn))
        sp=dplyr::pull(subset(sample, sa==i), f1[[2]])
        samp1=c(samp,sp )
        if(param=="both" | param=="mean") {ME[i]=ME[i]+mean(samp1)}
        if(param=="both" | param=="HCR")  {HCR[i]=HCR[i]+povinc(samp1, z=z)}

      }
      message('Processing estimation loop ', x, ' of ', R)

    }
    ME=ME/x
    HCR=HCR/x


  } else {
    gam1=gamlss::gamlss(f1, sigma.fo=f2, nu.fo=f3, tau.fo=f4, tau.fixed=substitute(tau.fixed), nu.fixed=substitute(nu.fixed),
               trace=F,family=substitute(fdis), data=sample, method=mixed(substitute(nRS), substitute(nRG)))

    yns=data.frame("ys"=as.vector(gamlss::predictAll(gam1,data=sample, newdata= nonsample)$mu),
                   "yss"=as.vector(gamlss::predictAll(gam1,data=sample, newdata=nonsample)$sigma),
                   "ysn"=as.vector(gamlss::predictAll(gam1,data=sample, newdata=nonsample)$nu),
                   "yst"=as.vector(gamlss::predictAll(gam1,data=sample, newdata=nonsample)$tau),
                   "sa"=nonsample$sa)
    ME=rep(0,D)
    HCR=rep(0,D)
    for (x in 1:R) {

      for (i  in 1:D) {
        yns1=subset(yns, sa==i)
        samp=as.vector(Dis(Ni[i]-ni[i], mu=yns1$ys, sigma=yns1$yss, nu=yns1$ysn, tau=yns1$yst))
        sp=dplyr::pull(subset(sample, sa==i), f1[[2]])
        samp1=c(samp,sp )
        if(param=="both" | param=="mean") {ME[i]=ME[i]+mean(samp1)}
        if(param=="both" | param=="HCR")  {HCR[i]=HCR[i]+povinc(samp1, z=z)}
      }
      message('Processing estimation loop ', x, ' of ', R)
    }
    ME=ME/x
    HCR=HCR/x
  }

  if (param =="both"){estim=list("ME"=ME,"HCR"= HCR)
  } else if (param =="mean"){estim=list("ME"=ME)
  } else {estim=list("HCR"= HCR)
  }

  input_var=list("fit"=gam1, "f1"=f1,"f2"=f2, "f3"=f3, "f4"=f4,"fdis"=fdis,
                 "np"=np, "nRS"=nRS, "nCG"=nCG, "R"=R, "D"=D, "Ni"=Ni, "ni"=ni,
                 "nu.fix"=nu.fix, "tau.fix"=tau.fix, "param"=param, "origindata"=sample, "z"=z)
  return(list("est"=estim, "input_var"=input_var))
}

