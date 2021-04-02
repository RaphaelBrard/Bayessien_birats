## Brard Raphaël
## Durand-Perdriel Flavie
## M'rabet Youssef
## Monceau Gaël

library(MASS)

# Initialisation
N = 30
t = 5
y = matrix(c(151, 145, 147, 155, 135, 159, 141, 159, 177, 134, 
              160, 143, 154, 171, 163, 160, 142, 156, 157, 152, 154, 139, 146, 
              157, 132, 160, 169, 157, 137, 153, 199, 199, 214, 200, 188, 210, 
              189, 201, 236, 182, 208, 188, 200, 221, 216, 207, 187, 203, 212, 
              203, 205, 190, 191, 211, 185, 207, 216, 205, 180, 200, 246, 249, 
              263, 237, 230, 252, 231, 248, 285, 220, 261, 220, 244, 270, 242, 
              248, 234, 243, 259, 246, 253, 225, 229, 250, 237, 257, 261, 248, 
              219, 244, 283, 293, 312, 272, 280, 298, 275, 297, 350, 260, 313, 
              273, 289, 326, 281, 288, 280, 283, 307, 286, 298, 267, 272, 285, 
              286, 303, 295, 289, 258, 286, 320, 354, 328, 297, 323, 331, 305, 
              338, 376, 296, 352, 314, 325, 358, 312, 324, 316, 317, 336, 321, 
              334, 302, 302, 323, 331, 345, 333, 316, 291, 324), nrow= 30 ,ncol= 5)
x = c(8.0, 15.0, 22.0, 29.0, 36.0)
xbar = 22
Omega = matrix(c(0.005, 0, 0, 5), nrow = 2) 
mu_beta = c(0,0)
sigmasq_y = 1
beta = matrix(c(100,100,100,100,100,
                    100,100,100,100,100,
                    100,100,100,100,100,
                    100,100,100,100,100,
                    100,100,100,100,100,
                    100,100,100,100,100,
                    6,6,6,6,6,
                    6,6,6,6,6,
                    6,6,6,6,6,
                    6,6,6,6,6,
                    6,6,6,6,6,
                    6,6,6,6,6),
                  nrow = 30, ncol=2)
Sigma_beta = matrix(c(1,0,0,1),
                       nrow = 2)


priors1 = c(0.001)
V = matrix(c(10**5,0,0,10**5), nrow =2)
R = matrix(c(200,0,0,0.2), nrow =2)
update_rho = dim(y)[1]+2





# Création de la fonction Birats :

birats = function(nchain,data,priors){
  init = c(1,0,0)
  chain = matrix(NA, nchain+1, 3)# columns tau , mu_beta1, mu_beta2
  colnames(chain) <- c("tau","mu_beta1","mu_beta2")
  chain[1,] <- init
  
  for(iter in 1:nchain){
    
    ### Mise a jour tau 
    update_alpha = priors[1]+ (dim(y)[1]*dim(y)[2])/2
    update_beta = 0
    for (i in 1:dim(y)[1]){
      for(j in 1:dim(y)[2]){
        update_beta =+ (y[i,j]-beta[i,1]-beta[i,2]*x[j])**2
      }
    }
    update_beta =0.5*update_beta + priors[1]
    
    tau = rgamma(1,update_alpha,update_beta)
    print(tau)
    # Mise a jour mu_beta
    update_mean = 0
    for (i in 1:dim(y)[1]){
      update_mean =+ Omega %*% as.matrix(beta[i,])
    }
    update_variance = solve(solve(V) + dim(y)[1]*Omega)
    
    mu_beta = as.matrix(mvrnorm(1, update_mean,update_variance))
    mu_beta1 = mu_beta[1]
    mu_beta2 = mu_beta[2]
   
    
    # Mise a jour de omega 
    update_R = matrix(0,nrow = 2 , ncol = 2)
    for(i in 1:dim(y)[1]){
      update_R = update_R +   (as.matrix(beta[i,])-mu_beta)%*%(t(as.matrix(beta[i,]))-t(mu_beta))
    }
     update_R = solve(solve(R)+update_R) 
     Omega = matrix(rWishart(1,update_rho,update_R),nrow = 2,ncol=2)
    
    
    #Mise a jour beta
     for (i in 1:dim(y)[1]){
       update1 = c(NA,NA)
       mat = matrix(NA, nrow = 2,ncol = 2)
       for (j in 1:t){
         update1 =+ y[i,j]*c(1,x[j])
         mat =+ matrix(c(1,1,1,x[j]**2),nrow = 2,ncol = 2)
       }
       update_mean1 = Omega %*% mu_beta - tau*t*update1
       
       update_variance1 = solve(Omega + tau*mat)
       beta[i,] = mvrnorm(1,update_mean1,update_variance1)
     }
    chain[iter+1,] <- c(tau, mu_beta1, mu_beta2)
  }
  
  return(chain)
  
}




chain = birats(10^4, y, priors = priors1)
chain = chain[-(1:1000),]



# Création des graphiques 

par(mfrow=c(1,3))
para = c("tau","mu_beta1","mu_beta2")
for (i in 1:3){
  plot(chain[,i], type="l", xlab="Iterations", ylab="", main=para[i], ylim = c(-50,50))
}









