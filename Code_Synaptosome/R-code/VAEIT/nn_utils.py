import numpy as np
import tensorflow as tf
import tensorflow_probability as tfp
tfd = tfp.distributions

from tensorflow.keras.layers import Layer, Dense, BatchNormalization, LeakyReLU, Lambda
use_bias = True



###########################################################################
#
# Sampling layers in the latent space
#
###########################################################################
class cdf_layer(Layer):
    '''
    The Normal cdf layer with custom gradients.
    '''
    def __init__(self):
        '''
        '''
        super(cdf_layer, self).__init__()
        
    @tf.function
    def call(self, x):
        return self.func(x)
        
    @tf.custom_gradient
    def func(self, x):
        '''Return cdf(x) and pdf(x).

        Parameters
        ----------
        x : tf.Tensor
            The input tensor.
        
        Returns
        ----------
        f : tf.Tensor
            cdf(x).
        grad : tf.Tensor
            pdf(x).
        '''   
        dist = tfp.distributions.Normal(
            loc = tf.constant(0.0, tf.keras.backend.floatx()), 
            scale = tf.constant(1.0, tf.keras.backend.floatx()), 
            allow_nan_stats=False)
        f = dist.cdf(x)
        def grad(dy):
            gradient = dist.prob(x)
            return dy * gradient
        return f, grad
    

class Sampling(Layer):
    """Sampling latent variable \(z\) from \(N(\\mu_z, \\log \\sigma_z^2\)).    
    Used in Encoder.
    """
    def __init__(self, seed=0, **kwargs):
        super(Sampling, self).__init__(**kwargs)
        self.seed = seed

    @tf.function
    def call(self, z_mean, z_log_var):
        '''Return cdf(x) and pdf(x).

        Parameters
        ----------
        z_mean : tf.Tensor
            \([B, L, d]\) The mean of \(z\).
        z_log_var : tf.Tensor
            \([B, L, d]\) The log-variance of \(z\).

        Returns
        ----------
        z : tf.Tensor
            \([B, L, d]\) The sampled \(z\).
        '''   
   #     seed = tfp.util.SeedStream(self.seed, salt="random_normal")
   #     epsilon = tf.random.normal(shape = tf.shape(z_mean), seed=seed(), dtype=tf.keras.backend.floatx())
        epsilon = tf.random.normal(shape = tf.shape(z_mean), dtype=tf.keras.backend.floatx())
        z = z_mean + tf.exp(0.5 * z_log_var) * epsilon
        z = tf.clip_by_value(z, -1e6, 1e6)
        return z



###########################################################################
#
# Encoder
# 
###########################################################################
class Encoder(Layer):
    '''
    Encoder, model \(p(Z_i|Y_i,X_i)\).
    '''
    def __init__(self, dimensions, dim_latent, name='encoder', **kwargs):
        '''
        Parameters
        ----------
        dimensions : np.array
            The dimensions of hidden layers of the encoder.
        dim_latent : int
            The latent dimension of the encoder.
        dim_block_inputs : list of int
            (num_block,) The dimension of each input block, where the last block 
            is assumed to be the batch effects.
        dim_block_latents : list of int
            (num_block,) The dimension of output of first layer for each block.
        block_names : list of str, optional
            (num_block,) The name of first layer for each block.  
        name : str, optional
            The name of the layer.
        **kwargs : 
            Extra keyword arguments.
        ''' 
        super(Encoder, self).__init__(name = name, **kwargs)
        self.input_layer = Dense(dimensions[0], use_bias=False, activation = LeakyReLU(), name='input')
        self.dense_layers = [Dense(dim, activation = LeakyReLU(),
                                          name = 'encoder_%i'%(i+1)) \
                             for (i, dim) in enumerate(dimensions[1:])]
        self.batch_norm_layers = [BatchNormalization(center=True) \
                                    for _ in range(len(dimensions[1:]))]
        self.batch_norm_layers.append(BatchNormalization(center=True))
        self.latent_mean = Dense(dim_latent, name = 'latent_mean')
        self.latent_log_var = Dense(dim_latent, name = 'latent_log_var')
        self.sampling = Sampling()
    
    @tf.function
    def call(self, x, embed, batches, L=1, training=True):
        '''Encode the inputs and get the latent variables.

        Parameters
        ----------
        x : tf.Tensor
            \([B, L, d]\) The input.
        L : int, optional
            The number of MC samples.
        training : boolean, optional
            Whether in the training or inference mode.
        
        Returns
        ----------
        z_mean : tf.Tensor
            \([B, L, d]\) The mean of \(z\).
        z_log_var : tf.Tensor
            \([B, L, d]\) The log-variance of \(z\).
        z : tf.Tensor
            \([B, L, d]\) The sampled \(z\).
        '''
        
        z = tf.concat([x, embed, batches], axis=-1)
        tmp = self.input_layer(z, training=training)
        _z = tmp
        for dense, bn in zip(self.dense_layers, self.batch_norm_layers):
            _z = dense(_z, training=training)
            _z = bn(_z, training=training)
        z_mean = self.batch_norm_layers[-1](
            self.latent_mean(_z, training=training), training=training)
        z_log_var = self.latent_log_var(_z)
        _z_mean = tf.tile(tf.expand_dims(z_mean, 1), (1,L,1))
        _z_log_var = tf.tile(tf.expand_dims(z_log_var, 1), (1,L,1))
        z = self.sampling(_z_mean, _z_log_var)
        return z_mean, z_log_var, z, tmp



###########################################################################
#
# Decoder
# 
###########################################################################
class Decoder(Layer):
    '''
    Decoder, model \(p(Y_i|Z_i,X_i)\).
    '''
    def __init__(self, dimensions, dim_output, name = 'decoder', **kwargs):
        '''
        Parameters
        ----------
        dimensions : np.array
            The dimensions of hidden layers of the encoder.
        dim_block_outputs : list of int
            (B,) The dimension of each output block.
        dist_block_outputs : list of str
            (B,) `'NB'`, `'ZINB'`, `'Bernoulli'` or `'Gaussian'`.
        dim_block_latents : list of int
            (B,) The dimension of output of last layer for each block.
        block_names : list of str, optional
            (B,) The name of last layer for each block.        
        name : str, optional
            The name of the layer.
        '''
        super(Decoder, self).__init__(name = name, **kwargs)       
        self.dense_layers = [Dense(dim, activation = LeakyReLU(),
                                          name = 'decoder_%i'%(i+1)) \
                             for (i,dim) in enumerate(dimensions[:-1])]
        self.output_layer = Dense(dimensions[-1], activation = LeakyReLU(), name='output')
        self.mean_layer = Dense(dim_output, name='mean')
        self.batch_norm_layers = [BatchNormalization(center=True) \
                                    for _ in range(len((dimensions)))]
        self.disp = Dense(dim_output, use_bias=False, activation = tf.nn.softplus, name="disp")  
       
    @tf.function  
    def call(self, x, embed, masks, batches, z, tmp, training=True, return_prob=True):
        '''Decode the latent variables and get the reconstructions.

        Parameters
        ----------
        z : tf.Tensor
            \([B, L, d]\) the sampled \(z\).
        training : boolean, optional
            whether in the training or inference mode.

        Returns
        ----------
        log_probs : tf.Tensor
            \([B, block]\) The log probability.
        '''
        L = tf.shape(z)[1]
        _z = tf.concat([z, tf.tile(tf.expand_dims(batches, 1), (1,L,1))], axis=-1)

        for dense, bn in zip(self.dense_layers, self.batch_norm_layers):
            _z = dense(_z)
            _z = bn(_z, training=training)

#         _z = self.batch_norm_layers[-1](self.output_layer(_z) + tf.expand_dims(tmp,1), training=training)
        _z = self.output_layer(_z) + tf.expand_dims(tmp,1)
        x_hat = self.mean_layer(_z)
        
        disp = self.disp(tf.expand_dims(batches ,1))
        generative_dist = tfd.Independent(tfd.Masked(
            tfd.Normal(
                loc = x_hat, scale = disp, name='Gaussian_rv'
            ), tf.expand_dims(masks, 1)), reinterpreted_batch_ndims=1)
            
        if return_prob:
            return generative_dist.log_prob(tf.expand_dims(x,1))
        else:
            return generative_dist.mean()




