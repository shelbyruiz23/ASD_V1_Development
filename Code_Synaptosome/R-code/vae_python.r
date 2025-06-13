
##################################################################
#
# VAE
#
##################################################################
# import Python package (call after setting up python environment with reticulate)
VAEIT <- import("VAEIT") 

#' run VAE
#'
#' @param Xna An p-by-n input matrix containing missing entries.
#' @param batches_cate A n-by-b1 matrix containing categorical/binary covariates.
#' @param batches_cont A n-by-b2 matrix containing continuous covariates.
#' @param X_init A p-by-n initialization matrix.
#' @param num_epoch A p-by-n initialization matrix.
#' @param print_every_epoch A p-by-n initialization matrix.
#' @param X (optional) The ground truth matrix for printing evaluating during training.
vae <- function(Xna, batches_cate, batches_cont, X_init, 
                num_epoch=20L, print_every_epoch=8L, X=NULL,
               beta_kl=1., beta_unobs=0.5, p_feat=0.2){
    data <- t(Xna)
    data[is.na(data)] <- 0.
    mask <- - as.matrix(is.na(t(Xna)))
    X_init <- t(X_init)
    config = list(
        # A network stucture of 
        # x     :              dim_input -> 64 -> 16 -> z 4 -> 16 -> 64 -> dim_input
        #                                 |                  |
        # masks : dim_input -> dim_embed ->                 ->
        'dim_input'=dim(data)[2],
        'dimensions'=c(64, 16), # hidden layers
        'dim_embed'=128L, 
        'dim_latent'=4L,

        # some hyperparameters
        'beta_kl'=beta_kl,
        'beta_unobs'=beta_unobs,
        'beta_reverse'=0.2,

        # prob of random maskings
        "p_feat"=p_feat

    )
    
    model <- VAEIT$VAEIT(config, data, mask, batches_cate, batches_cont)
    
    if(!is.null(X)){X <- t(X)}
    model$train(
        X_init,
        num_epoch=as.integer(num_epoch),# the number of iterations, generally ~10 would be 
                         # good. If this is too large, it may overfit the data.
        num_repeat=500L, # how many batches for each epoch
        batch_size=dim(data)[1],  # full batch of 48 samples
        learning_rate=3e-4,
        # if X is provided, then evaluate the model every 4 epochs
        save_every_epoch=as.integer(print_every_epoch),
        X=X,
        verbose=FALSE
    )
    mu <- matrix(as.vector(apply(t(X_init), 1, mean)), nrow=1)
               ##############################################
    sig <- 1   # this change was made in the 08-22-2022 version, this used to be: sig <- sd(X_init)
               ##############################################
    X_imp <- t(model$get_recon(mu, sig))

    X_blend <- Xna
    ina <- is.na(Xna)
    X_blend[ina] <- X_imp[ina]    
    
    list(X_blend=X_blend, X_imp=X_imp)
}


