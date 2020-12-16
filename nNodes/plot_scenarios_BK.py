#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 19:52:32 2020

@author: guido
"""
import matplotlib.pyplot as plt
import numpy as np

def read_data(file):
 
    sent_1 = []
    coll_1 = []
    rec_1 = []
    sent_3 = []
    coll_3 = []
    rec_3 = []
    with open(file) as file:
        data = file.read().splitlines()
    
    title_1 = data[2]
    x = data[5].split(",")
    for i in range(7,16):
        aux = data[i].split(",")
        sent_1.append(float(aux[0]))
        coll_1.append(float(aux[1]))
        rec_1.append(float(aux[4]))
    
    title_3 = data[17]
    for i in range(22,31):
        aux = data[i].split(",")
        sent_3.append(float(aux[0]))
        coll_3.append(float(aux[1]))
        rec_3.append(float(aux[4]))

    return data,x,title_1,sent_1,coll_1,rec_1,title_3,sent_3,coll_3,rec_3
    
    

data_file = "ET_nNodes.txt"
data,x,title_1,sent_1,coll_1,rec_1,title_3,sent_3,coll_3,rec_3 = read_data(data_file)
sent_1 = np.array(sent_1)
coll_1 = np.array(coll_1)
rec_1 = np.array(rec_1)
sent_3 = np.array(sent_3)
coll_3 = np.array(coll_3)
rec_3 = np.array(rec_3)

rec_1_norm = np.divide(rec_1, sent_1, out=np.zeros_like(rec_1), where=sent_1!=0)
coll_1_norm = np.divide(coll_1, sent_1, out=np.zeros_like(coll_1), where=sent_1!=0)

rec_3_norm = np.divide(rec_3, sent_3, out=np.zeros_like(rec_3), where=sent_3!=0)
coll_3_norm = np.divide(coll_3, sent_3, out=np.zeros_like(coll_3), where=sent_3!=0)

plt.figure(figsize=(6, 4), dpi= 80, facecolor='w', edgecolor='k')
x = [int(i) for i in x]
plt.title(title_1.split("_")[0])
plt.xlabel("N° of nodes")
plt.ylabel("count")
plt.grid()

plt.plot(x,rec_1_norm,'c',label="nrReceived/sent - 1 ch",marker="^")
plt.plot(x,coll_1_norm,'c',label="nrCollisions/sent - 1 ch",marker="v")


plt.plot(x,rec_3_norm,'b',label="nrReceived/sent - 3 ch",marker="^")
plt.plot(x,coll_3_norm,'b',label="nrCollisions/sent - 3 ch",marker="v")
plt.xlim(1400,6000)
plt.show()
plt.legend()
plt.savefig(title_1.split("_")[0]+".png")
