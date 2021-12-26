import scanpy as sc
import tensorflow.keras.layers as layers
import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
import sys
from tensorflow import keras
from sklearn.model_selection import train_test_split

#config
SIZE = 1024
config  = np.loadtxt(sys.argv[1],delimiter=',',dtype='str')
anndata_path,marker_path,threshold_quality,rae_epochs,pae_epochs,gen_epochs,save_path = config[:,1]

#data load
adata = sc.read(anndata_path)
df = np.loadtxt(marker_path,delimiter=',',dtype='str')
sc.pp.filter_genes(adata, min_counts=1)
adata.uns["rna_marker"] =df[:,0]
adata.uns["protein_marker"] =df[:,1]

protein = adata[:,adata.uns["protein_marker"]]
sc.pp.normalize_total(protein,target_sum=1e4)
sc.pp.scale(protein)

rna = adata[:,adata.var.feature_types.isin(["Gene Expression"])]
sc.pp.normalize_total(rna,target_sum=1e4)
sc.pp.scale(rna)

def qulity_check(rna,protein):
    qt = [[],[]] 
    def cossim(i,j):
        if np.linalg.norm(i) == 0 or np.linalg.norm(j) == 0:
            return 0
        else:
            return np.dot(i,j)/(np.linalg.norm(i)*np.linalg.norm(j))
    def eucdis(i,j):
            return np.linalg.norm(i-j)
    for i in range(rna.shape[0]):   
        qt[0].append(cossim(rna[i],protein[i]))
        qt[1].append(eucdis(rna[i],protein[i]))
    return np.array(qt)

adata.uns["marker_quality_before"] = qulity_check(rna[:,rna.uns["rna_marker"]].X,protein.X)
threshold = np.percentile(adata.uns["marker_quality_before"][0], int(threshold_quality), interpolation='nearest')
train_rna = rna[adata.uns["marker_quality_before"][0]>threshold,rna.uns["rna_marker"]].X
train_protein = protein[adata.uns["marker_quality_before"][0]>threshold].X
train_rna,val_rna,train_protein,val_protein = train_test_split(train_rna,train_protein)
 
def model_train_plot(model,path):
    _, loss_ax = plt.subplots()
    acc_ax = loss_ax.twinx()
    loss_ax.plot(model.history['loss'], 'y', label='train loss')
    loss_ax.plot(model.history['val_loss'], 'r', label='val loss')
    acc_ax.plot(model.history['mse'], 'b', label='train acc')
    acc_ax.plot(model.history['val_mse'], 'g', label='val acc')
    loss_ax.set_xlabel('epoch')
    loss_ax.set_ylabel('loss')
    acc_ax.set_ylabel('accuray')
    loss_ax.legend(loc='upper left')
    acc_ax.legend(loc='lower left')
    plt.savefig(path)
    plt.close()

#RNA encoding 1
size=train_rna.shape[1]
r_input = keras.Input(shape=size)
rx = layers.Dense(units=SIZE)(r_input)
rx = layers.Dropout(rate=0.3)(rx)
rx = layers.Dense(units=SIZE//2)(rx)
rx = layers.Dropout(rate=0.3)(rx)
r_encode = layers.Dense(units=SIZE//4)(rx)
rx = layers.Dense(units=SIZE//2)(rx)
rx = layers.Dense(units=SIZE)(rx)
r_out = layers.Dense(units=size)(rx)

rae = keras.Model(r_input,r_out)
r_encoder = keras.Model(r_input,r_encode)
rae.compile(optimizer=tf.optimizers.Adam(), loss='mse',metrics=['mse'])
r_encoder.compile(optimizer=tf.optimizers.Adam(), loss='mse')
rhist = rae.fit(train_rna,train_rna,validation_data=(val_rna,val_rna), epochs=int(rae_epochs), batch_size=256)
model_train_plot(rhist,"rae_train.png")

#protein encoding
size=train_protein.shape[1]
p_input = keras.Input(shape=size)
px = layers.Dense(units=SIZE)(p_input)
px = layers.Dropout(rate=0.3)(px)
px = layers.Dense(units=SIZE//2)(px)
px = layers.Dropout(rate=0.3)(px)
p_encode = layers.Dense(units=SIZE//4)(px)
px = layers.Dense(units=SIZE//2)(px)
px = layers.Dense(units=SIZE)(px)
p_out = layers.Dense(units=size)(px)

pae = keras.Model(p_input,p_out)
p_encoder = keras.Model(p_input,p_encode)
pae.compile(optimizer=tf.optimizers.Adam(), loss='mse',metrics=['mse'])
p_encoder.compile(optimizer=tf.optimizers.Adam(), loss='mse')
phist = pae.fit(train_protein,train_protein,validation_data=(val_protein,val_protein), epochs=int(pae_epochs), batch_size=256)
model_train_plot(phist,"pae_train.png")

#generator
in1=r_encoder(r_input)
in2=p_encoder(p_input)
size=train_protein.shape[1]
con = keras.layers.concatenate([in1,in2])
gx = keras.layers.Dense(units=256)(con)
gx = keras.layers.Dropout(0.3)(gx)
gx = keras.layers.Dense(units=1024)(gx)
gx = keras.layers.Dropout(0.3)(gx)
g_out = keras.layers.Dense(units=size)(gx)
gen = keras.Model([r_input,p_input],g_out)

r_encoder.trainable = False
p_encoder.trainable = False
gen.compile(optimizer=tf.optimizers.Adam(), loss='mse',metrics=['mse'])
ghist = gen.fit([train_rna,train_protein],train_protein,validation_data=([val_rna,val_protein],val_protein),epochs=int(gen_epochs), batch_size=256)
model_train_plot(ghist,"gen_train.png")

adata.uns["rna"] = rna.X
adata.uns["original_protein"] = protein.X
adata.uns["produced_protein"] = gen.predict([rna[:,rna.uns["rna_marker"]].X,protein.X])
adata.uns["marker_quality_after"] = qulity_check(rna[:,rna.uns["rna_marker"]].X,adata.uns["produced_protein"])
adata.write(save_path)
gen.save("gen.h5")

####ploting 
plt.plot(np.sort(adata.uns["marker_quality_before"][0]),label='before quality')
plt.plot(np.sort(adata.uns["marker_quality_after"][0]),label='after quality')
plt.legend()
plt.savefig("cos_sim.png")
plt.close()

plt.plot(np.sort(adata.uns["marker_quality_before"][1]),label='before quality')
plt.plot(np.sort(adata.uns["marker_quality_after"][1]),label='after quality')
plt.legend()
plt.savefig("euclid_distance.png")
plt.close()
