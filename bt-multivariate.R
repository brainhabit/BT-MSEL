#####################################################################################################################
# This script is run to test different multivariate twin modelling approaches
# using the BATSS sample for the paper: 
# The latent structure of emerging cognitive abilities: an infant twin study.
# Bussu G., Taylor M., Tammimies K., Ronald A., Falck-Ytter T.
#
# coded by Giorgia Bussu based on script from M. Taylor, June 2022
# 
####################################
# The script needs as input the selected dataset in the wide format (twin_data.csv from the bt-prepdata script).
#
# The script gives as output results from the different multivariate models:
# saturated model, correlated factors solution, independent pathway solution, and common pathway solution.
#
#####################################################################################################################

rm(list=ls())

setwd('your folder')
list.files()

### load OpenMx and helper functions (available upon request):
require(OpenMx)
source('your folder/miFunctions.R')

mxOption(key="Number of Threads",
         value=parallel::detectCores())

### load twin data:
data <- read.csv(file='twin_data.csv',header=T,sep=',')
names(data);dim(data)

Vars <- c('GMstd','FMstd','VRstd','RLstd','ELstd')
nv <- 5 # number of phenotypes
ntv <- nv*2 # number of measures/variables
(selVars <- paste(Vars,c(rep(1,nv),rep(2,nv)),sep=''))


mz <- subset(data,zygosity=='MZ',c(selVars))
dz <- subset(data,zygosity=='DZ',c(selVars))

############## Fully Saturated model ##################################################################################

svCov <- c(1,rep(.5,9),1,rep(.5,8),1,rep(.5,7),1,rep(.5,6),1,rep(.5,5),1,rep(.5,4),1,rep(.5,3),1,rep(.5,2),1,.5,1)


# Create Matrices for Covariates and linear Regression Coefficients
expMeanMZ <- mxMatrix(type='Full',nrow=1,ncol=ntv,free=T,values=0.01,
                      labels=labFull('m_mz',1,ntv),name='ExpMeanMZ')
expMeanDZ <- mxMatrix(type='Full',nrow=1,ncol=ntv,free=T,values=0.01,
                      labels=labFull('m_dz',1,ntv),name='ExpMeanDZ')


expCovMZ <- mxMatrix(type='Symm',nrow=ntv,ncol=ntv,free=T,values=svCov,
labels=labSymm('mz_cov',ntv),name='ExpCovMZ')

expCovDZ <- mxMatrix(type='Symm',nrow=ntv,ncol=ntv,free=T,values=svCov,
labels=labSymm('dz_cov',ntv),name='ExpCovDZ')

matI <- mxMatrix(type='Iden',nrow=ntv,ncol=ntv,name='I')

expCorMZ <- mxAlgebra(solve(sqrt(I*ExpCovMZ))%&%ExpCovMZ,name='ExpCorMZ')
expCorDZ <- mxAlgebra(solve(sqrt(I*ExpCovDZ))%&%ExpCovDZ,name='ExpCorDZ')

dataMZ <- mxData(mz,type='raw')
dataDZ <- mxData(dz,type='raw')

# objectives:
objMZ <- mxExpectationNormal(covariance='ExpCovMZ',means='ExpMeanMZ',dimnames=selVars)
objDZ <- mxExpectationNormal(covariance='ExpCovDZ',means='ExpMeanDZ',dimnames=selVars)

funcML <- mxFitFunctionML()

# Create Model Objects for Multiple Groups

# specify the submodels for each zygosity group:
modelMZ <- mxModel('MZ',expMeanMZ,expCovMZ,matI,expCorMZ,dataMZ,objMZ,funcML)
modelDZ <- mxModel('DZ',expMeanDZ,expCovDZ,matI,expCorDZ,dataDZ,objDZ,funcML)

# combine submodels into one object to make sure they are evaluated together:
multi <- mxFitFunctionMultigroup(c('MZ','DZ'))
ci <- mxCI(c('MZ.ExpCorMZ','DZ.ExpCorDZ'))

# combine all model objects:
SatModel <- mxModel('Sat',modelMZ,modelDZ,multi,ci)

mxOption(model= SatModel, key="Number of Threads", value= (omxDetectCores() - 1))

SatFit <- mxTryHard(SatModel,intervals=T)
sumSatFit<-summary(SatFit)

####################################### Constrained Model (estimate twin correlation) ############################################################

CorModel <- mxModel(SatFit,name='Cor')

# equal means across twins and zygosity:
CorModel <- omxSetParameters(CorModel,labels=c(labFull('m_mz',1,ntv),
labFull('m_dz',1,ntv)),
free=T,values=.001,newlabels=labFull('m',1,nv))

# equal variances:
CorModel <- omxSetParameters(CorModel,labels=c(labDiag('mz_cov',ntv),
labDiag('dz_cov',ntv)),
free=T,values=1,newlabels=labDiag('v',nv))

# equal phenotypic correlation in twin 1 and twin 2:
CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_2_1','mz_cov_7_6',
'dz_cov_2_1','dz_cov_7_6'),
free=T,values=.5,newlabel='rph_GMFM')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_3_1','mz_cov_8_6',
'dz_cov_3_1','dz_cov_8_6'),
free=T,values=.5,newlabel='rph_GMVR')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_4_1','mz_cov_9_6',
'dz_cov_4_1','dz_cov_9_6'),
free=T,values=.5,newlabel='rph_RLGM')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_5_1','mz_cov_10_6',
                                               'dz_cov_5_1','dz_cov_10_6'),
                             free=T,values=.5,newlabel='rph_GMEL')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_3_2','mz_cov_8_7',
                                               'dz_cov_3_2','dz_cov_8_7'),
                             free=T,values=.5,newlabel='rph_FMVR')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_4_2','mz_cov_9_7',
                                               'dz_cov_4_2','dz_cov_9_7'),
                             free=T,values=.5,newlabel='rph_RLFM')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_5_2','mz_cov_10_7',
                                               'dz_cov_5_2','dz_cov_10_7'),
                             free=T,values=.5,newlabel='rph_ELFM')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_4_3','mz_cov_9_8',
                                               'dz_cov_4_3','dz_cov_9_8'),
                             free=T,values=.5,newlabel='rph_RLVR')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_5_3','mz_cov_10_8',
                                               'dz_cov_5_3','dz_cov_10_8'),
                             free=T,values=.5,newlabel='rph_ELVR')
CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_5_4','mz_cov_10_9',
                                               'dz_cov_5_4','dz_cov_10_9'),
                             free=T,values=.5,newlabel='rph_ELRL')

# fix the CTCT to be the same between twins:
CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_6_2','mz_cov_7_1'),
free=T,values=.5,newlabel='ctct_mz_GMFM')
CorModel <- omxSetParameters(CorModel,labels=c('dz_cov_6_2','dz_cov_7_1'),
free=T,values=.5,newlabel='ctct_dz_GMFM')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_6_3','mz_cov_8_1'),
free=T,values=.5,newlabel='ctct_mz_GMVR')
CorModel <- omxSetParameters(CorModel,labels=c('dz_cov_6_3','dz_cov_8_1'),
free=T,values=.5,newlabel='ctct_dz_GMVR')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_6_4','mz_cov_9_1'),
                             free=T,values=.5,newlabel='ctct_mz_GMRL')
CorModel <- omxSetParameters(CorModel,labels=c('dz_cov_6_4','dz_cov_9_1'),
                             free=T,values=.5,newlabel='ctct_dz_GMRL')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_6_5','mz_cov_10_1'),
                             free=T,values=.5,newlabel='ctct_mz_GMEL')
CorModel <- omxSetParameters(CorModel,labels=c('dz_cov_6_5','dz_cov_10_1'),
                             free=T,values=.5,newlabel='ctct_dz_GMEL')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_7_4','mz_cov_9_2'),
free=T,values=.5,newlabel='ctct_mz_FMRL')
CorModel <- omxSetParameters(CorModel,labels=c('dz_cov_7_4','dz_cov_9_2'),
free=T,values=.5,newlabel='ctct_dz_FMRL')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_7_3','mz_cov_8_2'),
                             free=T,values=.5,newlabel='ctct_mz_FMVR')
CorModel <- omxSetParameters(CorModel,labels=c('dz_cov_7_3','dz_cov_8_2'),
                             free=T,values=.5,newlabel='ctct_dz_FMVR')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_7_5','mz_cov_10_2'),
                             free=T,values=.5,newlabel='ctct_mz_FMEL')
CorModel <- omxSetParameters(CorModel,labels=c('dz_cov_7_5','dz_cov_10_2'),
                             free=T,values=.5,newlabel='ctct_dz_FMEL')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_9_3','mz_cov_8_4'),
                             free=T,values=.5,newlabel='ctct_mz_RLVR')
CorModel <- omxSetParameters(CorModel,labels=c('dz_cov_9_3','dz_cov_8_4'),
                             free=T,values=.5,newlabel='ctct_dz_RLVR')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_10_3','mz_cov_8_5'),
                             free=T,values=.5,newlabel='ctct_mz_ELVR')
CorModel <- omxSetParameters(CorModel,labels=c('dz_cov_10_3','dz_cov_8_5'),
                             free=T,values=.5,newlabel='ctct_dz_ELVR')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_9_5','mz_cov_10_4'),
                             free=T,values=.5,newlabel='ctct_mz_RLEL')
CorModel <- omxSetParameters(CorModel,labels=c('dz_cov_9_5','dz_cov_10_4'),
                             free=T,values=.5,newlabel='ctct_dz_RLEL')

# fit the model
CorFit <- mxTryHard(CorModel,intervals=T)
sumCorFit<-summary(CorFit)

# compare with fully saturated to check assumptions
mxCompare(SatFit,CorFit)

########################################## Correlated factors solution ####################################################

### full model:

# path coefficients:
coefpath<-c(1,rep(.5,4),1,rep(.5,3),1,rep(.5,2),1,.5,1)

pathA <- mxMatrix(type='Lower',nrow=nv,ncol=nv,free=T,values=coefpath,
                  labels=labLower('a',nv),name='a')
pathC <- mxMatrix(type='Lower',nrow=nv,ncol=nv,free=T,values=coefpath,
                  labels=labLower('c',nv),name='c')  
pathE <- mxMatrix(type='Lower',nrow=nv,ncol=nv,free=T,values=coefpath,
                  labels=labLower('e',nv),name='e')  

# calculate the variance components:

varA <- mxAlgebra(a%*%t(a),name='A')
varC <- mxAlgebra(c%*%t(c),name='C')
varE <- mxAlgebra(e%*%t(e),name='E')

# calculate the total variance and the standard deviation:

varP <- mxAlgebra(A+C+E,name='V')

matI <- mxMatrix(type='Iden',nrow=nv,ncol=nv,name='I')
isd <- mxAlgebra(solve(sqrt(I*V)),name='iSD')

# calculate phenotypic and etiological correlations:

corP <- mxAlgebra(solve(sqrt(I*V))%&%V,name='rPH')
corA <- mxAlgebra(solve(sqrt(I*A))%&%A,name='rA')
corC <- mxAlgebra(solve(sqrt(I*C))%&%C,name='rC')
corE <- mxAlgebra(solve(sqrt(I*E))%&%E,name='rE')

# means:

expMean <- mxMatrix(type='Full',nrow=1,ncol=ntv,free=T,values=0.001,
                  labels=labFull('m',1,nv),name='ExpMean')


# variance/covariance:

expCovMZ <- mxAlgebra(rbind(cbind(V,A+C),
                            cbind(A+C,V)),name='ExpCovMZ')
expCovDZ <- mxAlgebra(rbind(cbind(V,0.5%x%A+C),
                            cbind(0.5%x%A+C,V)),name='ExpCovDZ')

# convert the variance components to proportions:

estVC <- mxAlgebra(cbind(A/V,C/V,E/V),name='EstVC')

# observed data:

dataMZ <- mxData(mz,type='raw')
dataDZ <- mxData(dz,type='raw')

# objectives:

objMZ <- mxExpectationNormal(covariance='ExpCovMZ',means='ExpMean',dimnames=selVars)
objDZ <- mxExpectationNormal(covariance='ExpCovDZ',means='ExpMean',dimnames=selVars)

# method of estimation:

funcML <- mxFitFunctionML()

# specify submodels for each zygosity group:

pars <- list(pathA,pathC,pathE,varA,varC,varE,varP,matI,isd,corP,corA,corC,corE,expMean)

modelMZ <- mxModel('MZ',pars,estVC,expCovMZ,dataMZ,objMZ,funcML)
modelDZ <- mxModel('DZ',pars,estVC,expCovDZ,dataDZ,objDZ,funcML)

# specify that we need to evaluate all submodels together:

multi <- mxFitFunctionMultigroup(c('MZ','DZ'))

# confidence intervals:

ci <- mxCI(c('MZ.EstVC','MZ.rPH','MZ.rA','MZ.rC','MZ.rE'))

# combine all objects:

CorFac<- mxModel('Cor',modelMZ,modelDZ,multi,ci)

# fit the model:

FacFit <- mxTryHard(CorFac,intervals=T)
sumCorFac<-summary(FacFit)

# compare fit to saturated model:

mxCompare(SatFit,FacFit)

######################################## Independent pathways Solution #######################################################################################################################################33

nf        <- 1       # number of factors

# Create Matrices for Covariates and linear Regression Coefficients
expMean <- mxMatrix(type='Full',nrow=1,ncol=ntv,free=T,values=0.001,
                  labels=labFull('m',1,nv),name='ExpMean')


# Matrices ac, cc, and ec to store a, d, and e path coefficients for common factors
pathAc    <- mxMatrix( type="Full", nrow=nv, ncol=nf, free=TRUE, values=.5, labels=labFull("ac",nv,nf), name="ac" )
pathCc    <- mxMatrix( type="Full", nrow=nv, ncol=nf, free=TRUE, values=.5, labels=labFull("cc",nv,nf), name="cc" )
pathEc    <- mxMatrix( type="Full", nrow=nv, ncol=nf, free=TRUE, values=.5, labels=labFull("ec",nv,nf), name="ec" )

# Matrices as, ds, and es to store a, d, and e path coefficients for specific factors
pathAs    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=1, labels=labDiag("as",nv), name="as" )
pathCs    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=1, labels=labDiag("cs",nv), name="cs" )
pathEs    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=1, labels=labDiag("es",nv), name="es" )

# Matrices A, D, and E compute variance components
varA      <- mxAlgebra( expression=ac %*% t(ac) + as %*% t(as), name="A" )
varC      <- mxAlgebra( expression=cc %*% t(cc) + cs %*% t(cs), name="C" )
varE      <- mxAlgebra( expression=ec %*% t(ec) + es %*% t(es), name="E" )

varP      <- mxAlgebra( expression= A+C+E, name="V" )


### expected variance/covariance:

expCovMZ <- mxAlgebra(rbind(cbind(V,A+C),
                            cbind(A+C,V)),name='ExpCovMZ')
expCovDZ <- mxAlgebra(rbind(cbind(V,0.5%x%A+C),
                            cbind(0.5%x%A+C,V)),name='ExpCovDZ')

# raw data
dataMZ <- mxData(mz,type='raw')
dataDZ <- mxData(dz,type='raw')

# objectives:
objMZ <- mxExpectationNormal(covariance='ExpCovMZ',means='ExpMean',dimnames=selVars)
objDZ <- mxExpectationNormal(covariance='ExpCovDZ',means='ExpMean',dimnames=selVars)

funcML <- mxFitFunctionML()

# calculate proportions of variance explained by common and specific 
# variance components:

estVC <- mxAlgebra(rbind(cbind(A/V,C/V,E/V),
                         cbind((ac%*%t(ac))/V,(cc%*%t(cc))/V,(ec%*%t(ec))/V),
                         cbind((as%*%t(as))/V,(cs%*%t(cs))/V,(es%*%t(es))/V)),
                   name='EstVC')

# Create Model Objects for Multiple Groups
pars      <- list(pathAc, pathCc, pathEc, pathAs, pathCs, pathEs, varA, varC, varE, varP, expMean)


modelMZ   <- mxModel( name="MZ", pars, expCovMZ, dataMZ, objMZ, funcML, estVC )
modelDZ   <- mxModel( name="DZ", pars, expCovDZ, dataDZ, objDZ, funcML, estVC )
multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )

### specify the confidence intervals:

ci <- mxCI(c('MZ.EstVC'))


# Build & Run Model 
modelIP   <- mxModel( "mulIPc", modelMZ, modelDZ, multi, ci )
fitIP     <- mxTryHard( modelIP, intervals=F)
sumIP     <- summary( fitIP )

mxCompare( FacFit, fitIP )


########################################### Common Pathway Solution ######################################################################################################################

nl        <- 1

# Create Matrices for Covariates and linear Regression Coefficients
expMean <- mxMatrix(type='Full',nrow=1,ncol=ntv,free=T,values=0.001,
                  labels=labFull('m',1,nv),name='ExpMean')

# Matrices ac, cc, and ec to store a, c, and e path coefficients for latent phenotype(s)
pathAl    <- mxMatrix( type="Lower", nrow=nl, ncol=nl, free=TRUE, values=.5, labels=labLower("al",nl), name="al" )
pathCl    <- mxMatrix( type="Lower", nrow=nl, ncol=nl, free=TRUE, values=.5, labels=labLower("cl",nl), name="cl" )
pathEl    <- mxMatrix( type="Lower", nrow=nl, ncol=nl, free=TRUE, values=.5, labels=labLower("el",nl), name="el" )


# Matrix and Algebra for constraint on variance of latent phenotype
varLP   <- mxAlgebra( expression=al %*% t(al) + cl %*% t(cl) + el %*% t(el), name="VarLP" )

unit      <- mxMatrix( type="Unit", nrow=nl, ncol=1, name="Unit")
varLP1    <- mxConstraint(diag2vec(VarLP)==Unit, name="varLP1")

# Matrix f for factor loadings on latent phenotype
pathFl    <- mxMatrix( type="Full", nrow=nv, ncol=nl, free=TRUE, values=1, labels=labFull("fl",nv,nl), name="fl" )

# Matrices as, cs, and es to store a, d, and e path coefficients for specific factors
pathAs    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=1, labels=labDiag("as",nv), name="as" )
pathCs    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=1, labels=labDiag("cs",nv), name="cs" )
pathEs    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=1, labels=labDiag("es",nv), name="es" )


# Matrices A, C, and E compute variance components
covA      <- mxAlgebra( expression=fl %&% (al %*% t(al)) + as %*% t(as), name="A" )
covC      <- mxAlgebra( expression=fl %&% (cl %*% t(cl)) + cs %*% t(cs), name="C" )
covE      <- mxAlgebra( expression=fl %&% (el %*% t(el)) + es %*% t(es), name="E" )

covP      <- mxAlgebra( expression= A+C+E, name="V" )

### expected variance/covariance:

expCovMZ <- mxAlgebra(rbind(cbind(V,A+C),
                            cbind(A+C,V)),name='ExpCovMZ')
expCovDZ <- mxAlgebra(rbind(cbind(V,0.5%x%A+C),
                            cbind(0.5%x%A+C,V)),name='ExpCovDZ')

# raw data
dataMZ <- mxData(mz,type='raw')
dataDZ <- mxData(dz,type='raw')

# objectives:
objMZ <- mxExpectationNormal(covariance='ExpCovMZ',means='ExpMean',dimnames=selVars)
objDZ <- mxExpectationNormal(covariance='ExpCovDZ',means='ExpMean',dimnames=selVars)

funcML <- mxFitFunctionML()

# calculate proportions:
estVC <- mxAlgebra(cbind(A/V,C/V,E/V),name='EstVC')

estVCl <- mxAlgebra(cbind((al%*%t(al))/VarLP,(cl%*%t(cl))/VarLP,(el%*%t(el))/VarLP),
                    name='EstVCl')
estVCc <- mxAlgebra(cbind((fl%&%(al%*%t(al)))/V,(fl%&%(cl%*%t(cl)))/V,
                          (fl%&%(el%*%t(el)))/V),name='EstVCc')
estVCs <- mxAlgebra(cbind((as%*%t(as))/V,(cs%*%t(cs))/V,(es%*%t(es))/V),name='EstVCs')


# Create Model Objects for Multiple Groups
pars      <- list(expMean, pathAl, pathCl, pathEl, varLP, unit, pathFl, pathAs, pathCs, pathEs, covA, covC, covE, covP)

modelMZ   <- mxModel( name="MZ", pars, expCovMZ, dataMZ, objMZ, funcML,estVC,estVCl,estVCc,estVCs )
modelDZ   <- mxModel( name="DZ", pars, expCovDZ, dataDZ, objDZ, funcML,estVC,estVCl,estVCc,estVCs )
multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )

### specify the confidence intervals:

ci <- mxCI(c('MZ.EstVC','MZ.EstVCl','MZ.EstVCc','MZ.EstVCs'))


# Build & Run Model 
modelCP   <- mxModel( "mulCPc", modelMZ, modelDZ, multi, ci )
fitCP    <- mxTryHard(modelCP, intervals=T)

sumCP     <- summary( fitCP )
mxCompare( FacFit,c(fitIP, fitCP ))

#####################################################################################################################
####### Save results
#####################################################################################################################

save.image("~/your folder/multivariate_BT_MSEL_modelling.RData")
