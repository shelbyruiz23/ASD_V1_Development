import numpy as np
import tensorflow as tf
import tensorflow_probability as tfp
tfd = tfp.distributions

from VAEIT.utils import ModalMaskGenerator
from VAEIT.nn_utils import Encoder, Decoder
from tensorflow.keras.layers import Layer, Dense, BatchNormalization, Lambda
from tensorflow.keras.utils import Progbar

            
class VariationalAutoEncoder(tf.keras.Model):
    """
    Combines the encoder, decoder and LatentSpace into an end-to-end model for training and inference.
    """
    def __init__(self, config, name = 'autoencoder', **kwargs):
        '''
        Parameters
        ----------
        dimensions : np.array
            The dimensions of hidden layers of the encoder.
        dim_latent : int
            The latent dimension of the encoder.
        
        dim_block : list of int
            (num_block,) The dimension of each input block.        
        dist_block : list of str
            (num_block,) `'NB'`, `'ZINB'`, `'Bernoulli'` or `'Gaussian'`.
        dim_block_enc : list of int
            (num_block,) The dimension of output of first layer of the encoder for each block.
        dim_block_dec : list of int
            (num_block,) The dimension of output of last layer of the decoder for each block.        
        block_names : list of str, optional
            (num_block,) The name of first layer for each block.
        name : str, optional
            The name of the layer.
        **kwargs : 
            Extra keyword arguments.
        '''
        super(VariationalAutoEncoder, self).__init__(name = name, **kwargs)
        self.config = config
        self.embed_layer = Dense(self.config.dim_embed, 
                                 activation = tf.nn.tanh, use_bias=False, name = 'embed') \
            if self.config.dim_embed > 0 else Lambda(lambda x,training: tf.identity(x))
        self.encoder = Encoder(self.config.dimensions, self.config.dim_latent)
        self.decoder = Decoder(self.config.dimensions[::-1], self.config.dim_input)
        
        self.mask_generator = ModalMaskGenerator(
            None, config.p_feat, None)
        
        
    def generate_mask(self, inputs, mask=None, p=None):
        return self.mask_generator(inputs, mask, p)
        
          
    def call(self, x, masks, batches, L=1, training=True):
        '''Feed forward through encoder, LatentSpace layer and decoder.

        Parameters
        ----------
        x_normalized : np.array
            \([B, G]\) The preprocessed data.
        c_score : np.array
            \([B, s]\) The covariates \(X_i\), only used when `has_cov=True`.
        x : np.array, optional
            \([B, G]\) The original count data \(Y_i\), only used when data_type is not `'Gaussian'`.
        scale_factor : np.array, optional
            \([B, ]\) The scale factors, only used when data_type is not `'Gaussian'`.
        pre_train : boolean, optional
            Whether in the pre-training phare or not.
        L : int, optional
            The number of MC samples.
        alpha : float, optional
            The penalty parameter for covariates adjustment.

        Returns
        ----------
        losses : float
            the loss.
        '''

        z_mean_obs, z_log_var_obs, z_obs, log_probs_obs = self._get_reconstruction_loss(
            x, masks!=-1., masks!=-1., batches, L, training=training)
        z_mean_unobs_1, z_log_var_unobs_1, z_unobs_1, log_probs_unobs = self._get_reconstruction_loss(
            x, masks==0., masks==1., batches, L, training=training)
        z_mean_unobs_2, z_log_var_unobs_2, z_unobs_2, log_probs_ = self._get_reconstruction_loss(
            x, masks==1., masks==0., batches, L, training=training)
        log_probs_unobs = (1-self.config.beta_reverse) * log_probs_unobs + self.config.beta_reverse*log_probs_

#         if alpha>0.0:
#             zero_in = tf.concat([tf.zeros([z.shape[0],1,z.shape[2]], dtype=tf.keras.backend.floatx()), 
#                                 tf.tile(tf.expand_dims(c_score,1), (1,1,1))], -1)
#             reconstruction_zero_loss = self._get_reconstruction_loss(x, zero_in, scale_factor, 1)
#             reconstruction_z_loss = (1-alpha)*reconstruction_z_loss + alpha*reconstruction_zero_loss
        
        self.add_loss(
            - (1-self.config.beta_unobs) * tf.reduce_sum(log_probs_obs)
        )
        self.add_loss(
            - self.config.beta_unobs * tf.reduce_sum(log_probs_unobs)
        )
        
        kl_1 = (1-self.config.beta_reverse) * self._get_kl_normal(
            z_mean_unobs_1, z_log_var_unobs_1, z_mean_obs, z_log_var_obs)
        self.add_loss(self.config.beta_kl * kl_1)
        kl_2 = self.config.beta_reverse * self._get_kl_normal(
            z_mean_unobs_2, z_log_var_unobs_2, z_mean_obs, z_log_var_obs)
        self.add_loss(self.config.beta_kl * kl_2)

        return self.losses
    
    @tf.function
    def _get_reconstruction_loss(self, x, bool_mask_in, bool_mask_out, batches, L, training=True):
        '''
        Parameters
        ----------
        bool_mask_in : tf.Tensor of type tf.bool
            False indicates missing.
        bool_mask_out : tf.Tensor of type tf.bool
            Compute likelihood for entries with value True.
        '''
        _masks = tf.where(bool_mask_in, 0., 1.)
        _x = tf.multiply(x, 1.-_masks)
        embed = self.embed_layer(_masks, training=training)
        z_mean, z_log_var, z, tmp = self.encoder(_x, embed, batches, L, training=training)       
        log_probs = tf.reduce_mean(
            self.decoder(x, embed, bool_mask_out, batches, z, tmp, training=training), axis=0)
        return z_mean, z_log_var, z, log_probs
    
    
    @tf.function
    def _get_kl_normal(self, mu_0, log_var_0, mu_1, log_var_1):
        kl = 0.5 * (
            tf.exp(tf.clip_by_value(log_var_0-log_var_1, -6., 6.)) + 
            (mu_1 - mu_0)**2 / tf.exp(tf.clip_by_value(log_var_1, -6., 6.)) - 1.
             + log_var_1 - log_var_0)
        return tf.reduce_mean(tf.reduce_sum(kl, axis=-1))
    
    @tf.function
    def _get_kl_loss(self, z, z_log_var, training=True):
        log_p_z = self.latent_space(z, training=training)

        E_qzx = - tf.reduce_mean(
            0.5 * self.config.dim_latent *
            (tf.math.log(tf.constant(2 * np.pi, tf.keras.backend.floatx())) + 1.0) +
            0.5 * tf.reduce_sum(z_log_var, axis=-1)
            )
        return log_p_z, E_qzx
    
    
    def get_recon(self, dataset_test, mu, sig, L=50):
        mu = tf.convert_to_tensor(mu, dtype=tf.float32)
        sig = tf.convert_to_tensor(sig, dtype=tf.float32)
        x_hat = []
        for x,m,b in dataset_test:
            _m = tf.where(m!=-1., 0., 1.)
            x = (x-mu)/sig
            x = tf.where(m==-1., 0., x)
            embed = self.embed_layer(_m, training=True)
            _, _, z, tmp = self.encoder(x, embed, b, L, True)
            _x_hat = tf.reduce_mean(
                self.decoder(x, embed, m!=-1., b, z, tmp, training=True, return_prob=False), axis=1)
            x_hat.append(_x_hat.numpy())
        x_hat = np.concatenate(x_hat)        

        return x_hat
    
    
    def get_z(self, dataset_test):
        '''Get \(q(Z_i|Y_i,X_i)\).

        Parameters
        ----------
        dataset_test : tf.Dataset
            Dataset containing (x, batches).

        Returns
        ----------
        z_mean : np.array
            \([B, d]\) The latent mean.
        '''        
        z_mean = []
        for x,m,b in dataset_test:
            m = tf.where(m!=-1., 0., 1.)
            embed = self.embed_layer(m, training=False)
            _z_mean, _, _, _ = self.encoder(x, embed, b, 1, training=True)         
            z_mean.append(_z_mean.numpy())
        z_mean = np.concatenate(z_mean)        

        return z_mean

