# -*- coding: utf-8 -*-
from typing import Optional

from VAEIT.utils import Early_Stopping
import numpy as np
import tensorflow as tf
from tensorflow.keras.utils import Progbar
import tensorflow_addons as tfa

from time import time

def clear_session():
    '''Clear Tensorflow sessions.
    '''
    tf.keras.backend.clear_session()
    return None


def pre_train(X_init, masks, dataset_train, dataset_valid, dataset_full, vae, 
              save_model, checkpoint_dir, save_every_epoch,
              learning_rate: float, L: int, alpha: float,
              num_epoch: int, num_step_per_epoch: int, 
              es_patience: int, es_tolerance: int, es_relative: bool,
              verbose: bool = True, X=None, sample_mask=True):
    '''Pretraining.

    Parameters
    ----------
    dataset_train : tf.Dataset
        The Tensorflow Dataset object.
    dataset_valid : tf.Dataset
        The Tensorflow Dataset object.
    vae : VariationalAutoEncoder
        The model.
    learning_rate : float
        The initial learning rate for the Adam optimizer.
    L : int
        The number of MC samples.
    alpha : float, optional
        The value of alpha in [0,1] to encourage covariate adjustment. Not used if there is no covariates.
    num_epoch : int
        The maximum number of epoches.
    num_step_per_epoch : int
        The number of step per epoch, it will be inferred from number of cells and batch size if it is None.            
    es_patience : int
        The maximum number of epoches if there is no improvement.
    es_tolerance : float
        The minimum change of loss to be considered as an improvement.
    es_relative : bool, optional
        Whether monitor the relative change of loss or not.        
    es_warmup : int, optional
        The number of warmup epoches.

    Returns
    ----------
    vae : VariationalAutoEncoder
        The pretrained model.
    '''
    mu = tf.convert_to_tensor(
        np.nanmean(X_init, axis=0, keepdims=True).astype(np.float32), dtype=tf.float32)
#     mu = np.zeros(X_init.shape, dtype=np.float32)
#     for i in range(4):
#         mu[12*i:12*(i+1)+1,:] = np.mean(X_init[12*i:12*(i+1)+1,:], axis=0, keepdims=True)
#     mu = tf.convert_to_tensor(mu, dtype=tf.float32)
    
#     sig = tf.convert_to_tensor(np.nanstd(X_init), dtype=tf.float32)
#     sig = tf.convert_to_tensor(np.nanstd(X_init, axis=0, keepdims=True), dtype=tf.float32)
    sig = tf.convert_to_tensor(1., dtype=tf.float32)
    masks = np.array(masks, dtype=np.float32)
    if X is not None:
        X = np.array(X, dtype=np.float32)
        
    optimizer = tfa.optimizers.AdamW(learning_rate=learning_rate, weight_decay=1e-4)
#     optimizer = tf.keras.optimizers.Adam(learning_rate=learning_rate)
    checkpoint = tf.train.Checkpoint(optimizer=optimizer, net=vae, step=tf.Variable(0),)
    manager = tf.train.CheckpointManager(checkpoint, checkpoint_dir, 
        max_to_keep=None if dataset_valid is None else es_patience+2)
    checkpoint.restore(manager.latest_checkpoint)
    if manager.latest_checkpoint:
        print("Restored from {}".format(manager.latest_checkpoint))
    else:
        print("Initializing from scratch.")
    
    loss_train = tf.metrics.Mean('train_loss', dtype=tf.float32)
    loss_train_list = {}
    loss_train_list['obs'] = [tf.metrics.Mean('train_loss_%d'%i, dtype=tf.float32) for i in range(1)]
    loss_train_list['unobs'] = [tf.metrics.Mean('train_loss_%d'%i, dtype=tf.float32) for i in range(1)]
    loss_train_list['kl'] = [tf.metrics.Mean('train_loss_kl_%d'%i, dtype=tf.float32) for i in range(1)]
    
    if dataset_valid is not None:
        loss_val = tf.metrics.Mean('val_loss', dtype=tf.float32)
        loss_val_list = {}
        loss_val_list['obs'] = [tf.metrics.Mean('val_loss_%d'%i, dtype=tf.float32) for i in range(1)]
        loss_val_list['unobs'] = [tf.metrics.Mean('val_loss_%d'%i, dtype=tf.float32) for i in range(1)]
        loss_val_list['kl'] = [tf.metrics.Mean('val_loss_kl_%d'%i, dtype=tf.float32) for i in range(1)]
        
    early_stopping = Early_Stopping(patience=es_patience, tolerance=es_tolerance, relative=es_relative)

    if not verbose:
        progbar = Progbar(num_epoch)
    
    start_time = time()
    epoch_start = int(checkpoint.step) + 1
    for epoch in range(epoch_start,num_epoch+epoch_start):
        checkpoint.step.assign_add(1)
        if verbose:
            progbar = Progbar(num_step_per_epoch)
            print('Train - Start of epoch %d' % (epoch,))
        else:
            if epoch%2==0 or epoch+1==num_epoch+epoch_start:
                progbar.update(epoch-epoch_start+1)

        # Iterate over the batches of the dataset.
        for step, (x, m, b) in enumerate(dataset_train):
            m = vae.generate_mask(x, m)
            x = (x - mu)/sig
            x = tf.where(m==-1., 0., x)
            
            if sample_mask:
                m = tf.where(
                    tf.convert_to_tensor(
                        (masks[np.random.choice(masks.shape[0],m.shape[0]),:]==-1.) & (m!=-1.), dtype=tf.bool),
                    1., m)

            with tf.GradientTape() as tape:
                losses = vae(x, m, b, L=L)#/tf.cast(tf.shape(x)[1], tf.float32)
                # Compute reconstruction loss
                loss = tf.reduce_sum(losses)
            grads = tape.gradient(loss, vae.trainable_weights,
                        unconnected_gradients=tf.UnconnectedGradients.ZERO)
            optimizer.apply_gradients(zip(grads, vae.trainable_weights))
            
            for i, l in enumerate(loss_train_list['obs']):
                l(losses[i])
            for i, l in enumerate(loss_train_list['unobs']):
                l(losses[i+1])
            loss_train_list['kl'][0](losses[-2]+losses[-1])
            loss_train(loss)

            if verbose:
                if (step+1)%10==0 or step+1==num_step_per_epoch:
                    progbar.update(step+1, [('Reconstructed Loss', float(loss))])
        
        if dataset_valid is not None:
            for step, (x, m, b) in enumerate(dataset_valid):
                m = vae.generate_mask(x, m, p=0.)
                x = (x - mu)/sig
                x = tf.where(m==-1., 0., x)
                losses = vae(x, m, b, L=L, training=True)

                for i, l in enumerate(loss_val_list['obs']):
                    l(losses[i])
                for i, l in enumerate(loss_val_list['unobs']):
                    l(losses[i+1])
                loss_val_list['kl'][0](losses[-2]+losses[-1])
                loss_val(losses[0])#tf.reduce_sum(losses))

        if verbose:
            if dataset_valid is not None:
                print('Epoch {}, Train ELBO: {:5.02f}, Val ELBO: {:5.02f}, Time elapsed: {} minutes'.\
                    format(epoch, float(loss_train.result()), 
                       float(loss_val.result()), round((time() - start_time) / 60, 2)))
                print(', '.join('{:>5.02f}'.format(l.result()) for key in loss_train_list for l in loss_train_list[key]))
                print(', '.join('{:>5.02f}'.format(l.result()) for key in loss_val_list for l in loss_val_list[key]))
            else:
                print('Epoch {}, Train ELBO: {:5.02f}, Time elapsed: {} minutes'.\
                    format(epoch, float(loss_train.result()), round((time() - start_time) / 60, 2)))
                print(', '.join('{:>5.02f}'.format(l.result()) for key in loss_train_list for l in loss_train_list[key]))
        
        if dataset_valid is not None:
            if save_model:
                save_path = manager.save()
                print("Saved checkpoint for epoch {}: {}".format(epoch, save_path))
            if early_stopping(float(loss_val.result())):
                print('Early stopping.')
                break
                
            loss_val.reset_states()
            [l.reset_states() for key in loss_val_list for l in loss_val_list[key]]
        else:
            if int(checkpoint.step) % save_every_epoch == 0:
                # print(checkpoint.step)
                if save_model:
                    save_path = manager.save()
                    print("Saved checkpoint for epoch {}: {}".format(epoch, save_path))
                
                x_hat = vae.get_recon(dataset_full, mu, sig, L=100)
                x_hat = x_hat * sig + mu
                
#                 if recomp:
#                     x_hat = np.where(masks==-1.,x_hat, X)
#                     mu = tf.reduce_mean(x_hat, axis=0, keepdims=True)
#                     # sig = tf.math.reduce_std(x_hat, axis=0, keepdims=True)
#                     x_hat = vae.get_recon(dataset_train, mu, sig, L=100)
#                     x_hat = x_hat * sig + mu
                if X is not None:
                    print(np.mean((x_hat - X)[masks==-1]**2))
        
        
        loss_train.reset_states()
        [l.reset_states() for key in loss_train_list for l in loss_train_list[key]]

    print('Train Done.')
    return vae
