import warnings
from typing import Optional, Union
from types import SimpleNamespace

import VAEIT.model as model 
import VAEIT.train as train
import tensorflow as tf

from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.model_selection import train_test_split
import numpy as np

import warnings
warnings.filterwarnings('ignore')

class VAEIT():
    """
    Variational Inference for Trajectory by AutoEncoder.
    """
    def __init__(self, config: SimpleNamespace, data, masks, batches, batches_cont,
        ):
        '''
        Get input data for model. Data need to be first processed using scancy and stored as an AnnData object
         The 'UMI' or 'non-UMI' model need the original count matrix, so the count matrix need to be saved in
         adata.layers in order to use these models.


        Parameters
        ----------
        config : SimpleNamespace
            Dict of config.

        Returns
        -------
        None.

        '''
        self.dict_method_scname = {
            'PCA' : 'X_pca',
            'UMAP' : 'X_umap',
            'TSNE' : 'X_tsne',
            'diffmap' : 'X_diffmap',
            'draw_graph' : 'X_draw_graph_fa'
        }

#         if model_type != 'Gaussian':
#             if adata_layer_counts is None:
#                 raise ValueError("need to provide the name in adata.layers that stores the raw count data")
#             if 'highly_variable' not in adata.var:
#                 raise ValueError("need to first select highly variable genes using scanpy")
#         if copy_adata:
#             self.adata = adata.copy()
#         else:
#             self.adata = adata

#         self._adata = sc.AnnData(X = self.adata.X, var = self.adata.var)
#         self._adata.obs = self.adata.obs
#         self._adata.uns = self.adata.uns

        self.data = np.array(data, dtype=np.float32)
        batches = np.array(batches)
        self.cat_enc = OneHotEncoder().fit(batches)
        self.batches = self.cat_enc.transform(batches).toarray()
        if batches_cont is not None:
            self.cont_enc = StandardScaler()            
            batches_cont = self.cont_enc.fit_transform(batches_cont)
            batches_cont = np.nan_to_num(batches_cont)
            batches_cont = np.array(batches_cont, dtype=np.float32)            
            self.batches = np.c_[self.batches, batches_cont]
            
        self.masks = np.array(masks, dtype=np.float32)
        
        if isinstance(config, dict):
            from types import SimpleNamespace
            config = SimpleNamespace(**config)
        
        self.config = config
        self.reset()
        

    def train(self, X_init, 
              valid = False, stratify = False, test_size = 0.1, random_state: int = 0,
              learning_rate: float = 1e-3, batch_size: int = 256, batch_size_inference: int = 512, 
              L: int = 1, alpha: float = 0.10,
              num_epoch: int = 200, num_step_per_epoch: Optional[int] = None, num_repeat=1., shuffle=False,
              early_stopping_patience: int = 10, early_stopping_tolerance: float = 1e-4, 
              early_stopping_relative: bool = True, verbose: bool = False,
              save_model = False, checkpoint_dir: Optional[str] = None, delete_existing: Optional[str] = True, save_every_epoch = 50,
              sample_mask=True, X=None):
        '''Pretrain the model with specified learning rate.

        Parameters
        ----------
        test_size : float or int, optional
            The proportion or size of the test set.
        random_state : int, optional
            The random state for data splitting.
        learning_rate : float, optional
            The initial learning rate for the Adam optimizer.
        batch_size : int, optional 
            The batch size for pre-training.  Default is 256. Set to 32 if number of cells is small (less than 1000)
        L : int, optional 
            The number of MC samples.
        alpha : float, optional
            The value of alpha in [0,1] to encourage covariate adjustment. Not used if there is no covariates.
        num_epoch : int, optional 
            The maximum number of epochs.
        num_step_per_epoch : int, optional 
            The number of step per epoch, it will be inferred from number of cells and batch size if it is None.            
        early_stopping_patience : int, optional 
            The maximum number of epochs if there is no improvement.
        early_stopping_tolerance : float, optional 
            The minimum change of loss to be considered as an improvement.
        early_stopping_relative : bool, optional
            Whether monitor the relative change of loss as stopping criteria or not.
        path_to_weights : str, optional 
            The path of weight file to be saved; not saving weight if None.
        '''
        # if valid:
        #     if stratify is False:
        #         stratify = None    

        #     id_train, id_valid = train_test_split(
        #                             np.arange(self.data.shape[0]),
        #                             test_size=test_size,
        #                             stratify=stratify,
        #                             random_state=random_state)

        #     self.dataset_train = tf.data.Dataset.from_tensor_slices((
        #         self.data[id_train].astype(tf.keras.backend.floatx()), 
        #         self.masks[id_train].astype(tf.keras.backend.floatx()),
        #         self.batches[id_train].astype(tf.keras.backend.floatx())
        #     )).repeat(num_repeat)
        #     self.dataset_valid = tf.data.Dataset.from_tensor_slices((
        #             self.data[id_valid].astype(tf.keras.backend.floatx()), 
        #             self.masks[id_valid].astype(tf.keras.backend.floatx()),
        #             self.batches[id_valid].astype(tf.keras.backend.floatx())
        #     )).batch(batch_size_inference).prefetch(tf.data.experimental.AUTOTUNE)
        # else:            
        id_train = np.arange(self.data.shape[0])

        self.dataset_train = tf.data.Dataset.from_tensor_slices((
            self.data.astype(tf.keras.backend.floatx()),
            self.masks.astype(tf.keras.backend.floatx()),
            self.batches.astype(tf.keras.backend.floatx())
        )).repeat(num_repeat)
        
                
        if shuffle:
            self.dataset_train = self.dataset_train.shuffle(
                buffer_size = len(id_train)*num_repeat, seed=0, reshuffle_each_iteration=True)
        self.dataset_train = self.dataset_train.batch(batch_size).prefetch(tf.data.experimental.AUTOTUNE)
        if valid:
            self.dataset_valid = tf.data.Dataset.from_tensor_slices((
                self.data.astype(tf.keras.backend.floatx()),
                self.masks.astype(tf.keras.backend.floatx()),
                self.batches.astype(tf.keras.backend.floatx())
            )).batch(batch_size).prefetch(tf.data.experimental.AUTOTUNE)
        else:
            self.dataset_valid = None
            
        self.dataset_full = tf.data.Dataset.from_tensor_slices((
            self.data.astype(tf.keras.backend.floatx()),
            self.masks.astype(tf.keras.backend.floatx()),
            self.batches.astype(tf.keras.backend.floatx())
            )).batch(batch_size).prefetch(tf.data.experimental.AUTOTUNE)
            
        if num_step_per_epoch is None:
            num_step_per_epoch = (len(id_train)*num_repeat)//batch_size+1
            
            
        checkpoint_dir = 'checkpoint/pretrain/' if checkpoint_dir is None else checkpoint_dir
        if save_model:
            if delete_existing and tf.io.gfile.exists(checkpoint_dir):
                print("Deleting old log directory at {}".format(checkpoint_dir))
                tf.io.gfile.rmtree(checkpoint_dir)
            tf.io.gfile.makedirs(checkpoint_dir)
        
        self.vae = train.pre_train(
            X_init, self.masks,
            self.dataset_train,
            self.dataset_valid,
            self.dataset_full,
            self.vae,
            save_model,
            checkpoint_dir,
            save_every_epoch,
            learning_rate,                        
            L, alpha,
            num_epoch,
            num_step_per_epoch,
            early_stopping_patience,
            early_stopping_tolerance,
            early_stopping_relative,
            verbose, X, sample_mask)
    
    def get_recon(self, mu, sig=1., L=100):
        if not hasattr(self, 'dataset_full'):
            self.dataset_full = tf.data.Dataset.from_tensor_slices((
                self.data.astype(tf.keras.backend.floatx()),
                self.masks.astype(tf.keras.backend.floatx()),
                self.batches.astype(tf.keras.backend.floatx())
                )).batch(batch_size).prefetch(tf.data.experimental.AUTOTUNE)
        x_hat = self.vae.get_recon(self.dataset_full, mu, sig, L=L)
        x_hat = x_hat * sig + mu
        return x_hat
    
    
    def reset(self):
        train.clear_session()
        if hasattr(self, 'vae'):
            del self.vae
            import gc
            gc.collect()
        self.vae = model.VariationalAutoEncoder(self.config)
    
    def save_model(self, path_to_weights):
        checkpoint = tf.train.Checkpoint(net=self.vae)
        manager = tf.train.CheckpointManager(
            checkpoint, path_to_weights, max_to_keep=None
        )
        save_path = manager.save()        
        print("Saved checkpoint: {}".format(save_path))
        
    def load_model(self, path_to_weights):
        checkpoint = tf.train.Checkpoint(net=self.vae)
        status = checkpoint.restore(path_to_weights)
        print("Loaded checkpoint: {}".format(status))

            
    def get_latent_z(self, batch_size_inference=512):
        ''' get the posterier mean of current latent space z (encoder output)

        Returns
        ----------
        z : np.array
            \([N,d]\) The latent means.
        ''' 
        if not hasattr(self, 'dataset_full'):
            self.dataset_full = tf.data.Dataset.from_tensor_slices((
                    self.data.astype(tf.keras.backend.floatx()),
                    self.masks.astype(tf.keras.backend.floatx()),
                    self.batches.astype(tf.keras.backend.floatx()),
                )).batch(batch_size_inference).prefetch(tf.data.experimental.AUTOTUNE)

        return self.vae.get_z(self.dataset_full)
            
    
    def visualize_latent(self, method: str = "UMAP", 
                         color = None, **kwargs):
        '''
        visualize the current latent space z using the scanpy visualization tools

        Parameters
        ----------
        method : str, optional
            Visualization method to use. The default is "draw_graph" (the FA plot). Possible choices include "PCA", "UMAP", 
            "diffmap", "TSNE" and "draw_graph"
        color : TYPE, optional
            Keys for annotations of observations/cells or variables/genes, e.g., 'ann1' or ['ann1', 'ann2'].
            The default is None. Same as scanpy.
        **kwargs :  
            Extra key-value arguments that can be passed to scanpy plotting functions (scanpy.pl.XX).   

        Returns
        -------
        None.

        '''
          
        if method not in ['PCA', 'UMAP', 'TSNE', 'diffmap', 'draw_graph']:
            raise ValueError("visualization method should be one of 'PCA', 'UMAP', 'TSNE', 'diffmap' and 'draw_graph'")
        
        temp = list(self.adata.obsm.keys())
        if method == 'PCA' and not 'X_pca' in temp:
            print("Calculate PCs ...")
            sc.tl.pca(self.adata)
        elif method == 'UMAP' and not 'X_umap' in temp:  
            print("Calculate UMAP ...")
            sc.tl.umap(self.adata)
        elif method == 'TSNE' and not 'X_tsne' in temp:
            print("Calculate TSNE ...")
            sc.tl.tsne(self.adata)
        elif method == 'diffmap' and not 'X_diffmap' in temp:
            print("Calculate diffusion map ...")
            sc.tl.diffmap(self.adata)
        elif method == 'draw_graph' and not 'X_draw_graph_fa' in temp:
            print("Calculate FA ...")
            sc.tl.draw_graph(self.adata)
            

#         self._adata.obsp = self.adata.obsp
#        self._adata.uns = self.adata.uns
#         self._adata.obsm = self.adata.obsm
    
        if method == 'PCA':
            axes = sc.pl.pca(self.adata, color = color, **kwargs)
        elif method == 'UMAP':            
            axes = sc.pl.umap(self.adata, color = color, **kwargs)
        elif method == 'TSNE':
            axes = sc.pl.tsne(self.adata, color = color, **kwargs)
        elif method == 'diffmap':
            axes = sc.pl.diffmap(self.adata, color = color, **kwargs)
        elif method == 'draw_graph':
            axes = sc.pl.draw_graph(self.adata, color = color, **kwargs)
            
        return axes


 

    