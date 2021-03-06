#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 11:08:09 2020

@author: alvarezguido
GITHUB: https://github.com/alvarezguido
"""

"""
SYNOPSIS
----
----
-----

"""

import simpy
import random
import numpy as np
import math
#import sys
#import re
import matplotlib.pyplot as plt
#import os
#import operator
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
#import PIL
import random
import re


####WE START BY USING SF=12 ADN BW=125 AND CR=1, FOR ALL NODES AND ALL TRANSMISIONS######
####WE ALSO CONSIDER SIMPLE CHECK, WHERE TWO PACKETS COLLIDE WHEN THEY ARRIVE AT: SAME TIME, SAME FREQUENCY AND SAME SF####
nrNodes = 50 ##NUMBER OF NODES TO BE SIMULATED (IN ORDER FROM CSV FILE)
#multi_nodes = [1400,1000,500,250,100,50,25,10,5]
RANDOM_SEED = 6
random.seed(RANDOM_SEED) #RANDOM SEED IS FOR GENERATE ALWAYS THE SAME RANDOM NUMBERS (ie SAME RESULTS OF SIMULATION)

###PLOTS ##
plots_nodes = 0 ## FOR PLOT SIMULATION
plots_bar = 0 ##FOR PLOT BARS RESULTS OF SIMULATION

full_collision = False

###GLOBAL PARAMS ####
bsId = 1 ##ID OF BASE STATION (NOT USED)
channel = [0,1,2] ##NOT USED BY NOW

avgSendTime = 3  ## NOT USED! --> A NODE SENDS A PACKET EVERY X SECS
packetlen = 20   ##NODES SEND PACKETS OF JUST 20 Bytes
total_data = 60 ##TOTAL DATA ON BUFFER, FOR EACH NODE (IT'S THE BUFFER O DATA BEFORE START SENDING)

beacon_time = 120 ###SAT SENDS BEACON EVERY CERTAIN TIME
back_off = beacon_time * 0.95 ###BACK OFF TIME FOR SEND A PACKET
back_off = 110 ##NEW BACK-OFF TIME

packetsAtBS = [] ##USED FOR CHEK IF THERE ARE ALREADY PACKETS ON THE SATELLITE
c = 299792.458 ###SPEED LIGHT [km/s]
Ptx = 14
G_device = 0; ##ANTENNA GAIN FOR AN END-DEVICE
G_sat = 12;   ##ANTENNA GAIN FOR SATELLITE
nodes = [] ###EACH NODE WILL BE APPENDED TO THIS VARIABLE
freq =868e6 ##USED FOR PATH LOSS CALCULATION
frequency = [868100000, 868300000, 868500000] ##FROM LORAWAN REGIONAL PARAMETERS EU863-870 / EU868
#frequency = [868100000,868100000,868100000]
maxBSReceives = 1 ##MAX NUMBER OF PACKETS THAT BS (ie SATELLITE) CAN RECEIVE AT SAME TIME

nrLost = 0 ### TOTAL OF LOST PACKETS DUE Lpl
nrCollisions = 0 ##TOTAL OF COLLIDED PACKETS
nrProcessed = 0 ##TOTAL OF PROCESSED PACKETS
nrReceived = 0 ###TOTAL OF RECEIVED PACKETS

env = simpy.Environment()




##ARRAY WITH MEASURED VALUES FOR SENSIBILITY, NEW VALUES
##THE FOLLOWING VALUES CORRESPOND TO:
#   - FIRST ELEMENT: IT'S THE SF (NOT USABLE)
#   - SECOND ELEMENT: SENSIBILITY FOR 125KHZ BW
#   - THIRD ELEMENT: SENSIBILITY FOR 250KHZ BW
#   - FOURTH ELEMENT: SENSIBILITY FOR 500KHZ BW
# NOTICE THAT SENSIBILITY DECREASE ALONG BW INCREASES, ALSO WITH LOWER SF
# THIS VALUES RESPONDS TO:
# wf = -174 + 10 log(BW) +NF +SNRf
sf7 = np.array([7,-123,-120,-117.0])
sf8 = np.array([8,-126,-123,-120.0])
sf9 = np.array([9,-129,-126,-123.0])
sf10 = np.array([10,-132,-129,-126.0])
sf11 = np.array([11,-134.53,-131.52,-128.51])
sf12 = np.array([12,-137,-134,-131.0])

sensi = np.array([sf7,sf8,sf9,sf10,sf11,sf12])

## READ PARAMS FROM DIRECTORY ##
path = "./narrow_scenario/"

### -137dB IS THE MINIMUN TOLERABLE SENSIBILITY, FOR SF=12 AND BW=125KHz ###

leo_pos=np.loadtxt( path + "LEO-XYZ-Pos.csv",skiprows=1,delimiter=',',usecols=(1,2,3))
## WHERE:
    ## leo_pos[i,j]:
        ## i --> the step time in sat pass
        ## j --> 0 for x-position, 1 for y-position, 2 for z-position

sites_pos = np.loadtxt( path + "SITES-XYZ-Pos.csv",skiprows=1,delimiter=',',usecols=(1,2,3))
## WHERE:
    ## sites_pos[i,j]:
        ## i --> the node i
        ## j --> 0 for x-position, 1 for y-position, 2 for z-position


dist_sat = np.zeros((sites_pos.shape[0],3,leo_pos.shape[0]))
t = 0
for i in range(leo_pos.shape[0]):
    t+=1
    dist_sat [:,:,i] = leo_pos[i,:] - sites_pos
## WHERE:
    ## dist_sat[i,j,k]:
        ## i --> the node i
        ## j --> 0 for x-position, 1 for y-position, 2 for z-position
        ## k --> the step time in sat pass
    
#### FOR COMPUTE DISTANCE MAGNITUDE (ABS) FROM END-DEVICE TO SAT PASSING BY ####
distance = np.zeros((sites_pos.shape[0],leo_pos.shape[0]))
distance[:,:] = (dist_sat[:,0,:]**2 + dist_sat[:,1,:]**2 + dist_sat[:,2,:]**2)**(1/2)
## WHERE:
    ## distance[i,j]:
        ## i --> the node i
        ## j --> the step time in sat pass

##MATRIX FOR LINK BUDGET Lpl ###
Lpl = np.zeros((sites_pos.shape[0],leo_pos.shape[0])) 
Lpl = 20*np.log10(distance*1000) + 20*np.log10(freq) - 147.55 #DISTANCE MUST BE IN METERS
## WHERE:
    ## Lpl[i,j]:
        ## i --> the node i
        ## j --> the step time in sat pass 

##MATRIX FOR LINK BUDGET, USING Prx ###
Prx = np.zeros((sites_pos.shape[0],leo_pos.shape[0])) 
Prx = Ptx + G_sat + G_device -20*np.log10(distance*1000) - 20*np.log10(freq) + 147.55 #DISTANCE IS CONVERTED TO METERS
## WHERE:
    ## Prx[i,j]:
        ## i --> the node i
        ## j --> the step time in sat pass 

def checkcollision(packet):
    col = 0 # flag needed since there might be several collisions for packet
    processing = 0
    #print ("MAX RECEIVE IS: ", maxBSReceives)
    for i in range(0,len(packetsAtBS)):
        if packetsAtBS[i].packet.processed == 1:
            processing = processing + 1
    if (processing > maxBSReceives):
        print ("------{:3.5f} || Too much packets on Base Sattion.. Packet will be lost!", len(packetsAtBS))
        packet.processed = 0
    else:
        packet.processed = 1

    if packetsAtBS:
        print ("{:3.5f} || >> FOUND overlap... node {} (sf:{} bw:{} freq:{}) others: {}".format(env.now,packet.nodeid, packet.sf, packet.bw,packet.freq,len(packetsAtBS)))
        for other in packetsAtBS:
            if other.nodeid != packet.nodeid:
               print ("{:3.5f} || >> node {} overlapped with node {} (sf:{} bw:{} freq:{}). Let's check Freq...".format(env.now,packet.nodeid, other.nodeid, other.packet.sf, other.packet.bw,other.packet.freq))
               # simple collision
               #if frequencyCollision(packet, other.packet) and sfCollision(packet, other.packet):
               if frequencyCollision(packet, other.packet) and sfCollision(packet, other.packet):
# =============================================================================
#                     if timingCollision(packet, other.packet):
#                        # check who collides in the power domain
#                        c = powerCollision(packet, other.packet)
#                        # mark all the collided packets
#                        # either this one, the other one, or both
#                        for p in c:
#                            p.collided = 1
#                            if p == packet:
#                                col = 1
# =============================================================================
                
                    packet.collided = 1
                    other.packet.collided = 1  # other also got lost, if it wasn't lost already
                    col = 1
                                   
        return col
    return 0


###frequencyCollision, CONDITIONS###

##|f1-f2| <= 120 kHz if f1 or f2 has bw 500
##|f1-f2| <= 60 kHz if f1 or f2 has bw 250
##|f1-f2| <= 30 kHz if f1 or f2 has bw 125
def frequencyCollision(p1,p2):
    if (abs(p1.freq-p2.freq)<=120 and (p1.bw==500 or p2.freq==500)):
        print ("{:3.5f} || >> freq coll on node {} and node {}.. Let's check SF...".format(env.now,p1.nodeid, p2.nodeid))
        return True
    elif (abs(p1.freq-p2.freq)<=60 and (p1.bw==250 or p2.freq==250)):
        print ("{:3.5f} || >> freq coll on node {} and node {}.. Let's check SF...".format(env.now,p1.nodeid, p2.nodeid))
        return True
    else:
        if (abs(p1.freq-p2.freq)<=30):
            print( "{:3.5f} || >> Freq coll on node {} and node {}.. Let's check SF...".format(env.now,p1.nodeid, p2.nodeid))
            return True
        #else:
    print ("{:3.5f} || >> No frequency collision..".format(env.now))
    return False

#FOLLOWING FUNCTION NOT USED
def channelCollision(p1,p2):
    if (p1.ch == p2.ch):
        print ("{:3.5f} || >> channel coll for ch {} on node {} and ch {} on node {}.. Let's check SF...".format(env.now,p1.ch,p1.nodeid,p2.ch,p2.nodeid))
        return True
    else:
        print ("{:3.5f} || >> No channel collision..".format(env.now))
        return False

def sfCollision(p1, p2):
    if p1.sf == p2.sf:
        print ("{:3.5f} || >> COLLISION! SF coll on node {} and node {} (ie same SF)...".format(env.now,p1.nodeid, p2.nodeid))
        # p2 may have been lost too, will be marked by other checks
        return True
    print ("{:3.5f} || >> No SF Collision!".format(env.now))
    return False

def timingCollision(p1, p2):
    # assuming p1 is the freshly arrived packet and this is the last check
    # we've already determined that p1 is a weak packet, so the only
    # way we can win is by being late enough (only the first n - 5 preamble symbols overlap)

    # assuming 8 preamble symbols
    Npream = 8

    # we can lose at most (Npream - 5) * Tsym of our preamble
    Tpreamb = 2**p1.sf/(1.0*p1.bw) * (Npream - 5)

    # check whether p2 ends in p1's critical section
    p2_end = p2.addTime + p2.rectime
    p1_cs = env.now + (Tpreamb/1000.0)  # to sec
    ##print ("{} || >> collision timing node {} ({},{},{}) node {} ({},{})".format(env.now,p1.nodeid, env.now - env.now, p1_cs - env.now, p1.rectime,p2.nodeid, p2.addTime - env.now, p2_end - env.now))
    if p1_cs < p2_end:
        # p1 collided with p2 and lost
        print ("{:3.5f} || not late enough.. Timing collision...".format(env.now))
        return True
    print ("{:3.5f} || Saved by the preamble.. No timing collision!".format(env.now))
    return False

def powerCollision(p1, p2):
    powerThreshold = 6 # dB
    print ("{:3.5f} || power: node {} {:3.2f} dBm, node {} {:3.2f}; diff is {}dBm".format(env.now,p1.nodeid,p1.rssi[math.ceil(env.now)],p2.nodeid, p2.rssi[math.ceil(env.now)], round(p1.rssi[math.ceil(env.now)] - p2.rssi[math.ceil(env.now)],2)))
    #print ("pwr: node {0.nodeid} {0.rssi:3.2f} dBm node {1.nodeid} {1.rssi:3.2f} dBm; diff {2:3.2f} dBm".format(p1, p2, round(p1.rssi[math.ceil(env.now)] - p2.rssi[math.ceil(env.now)],2)))
    if abs(p1.rssi[math.ceil(env.now)] - p2.rssi[math.ceil(env.now)]) < powerThreshold:
        print( "{:3.5f} || Collision power both node {} and node {}".format(env.now,p1.nodeid, p2.nodeid))
        # packets are too close to each other, both collide
        # return both packets as casualties
        return (p1, p2)
    elif p1.rssi[math.ceil(env.now)] - p2.rssi[math.ceil(env.now)] < powerThreshold:
        # p2 overpowered p1, return p1 as casualty
        print ("{:3.5f} || Collision pwr node {} has overpowered node {}".format(env.now,p2.nodeid, p1.nodeid))
        return (p1,)
    print ("{:3.5f} || p1 wins, p2 lost".format(env.now))
    # p2 was the weaker packet, return it as a casualty
    return (p2,)

class myNode():
    def __init__(self, nodeid, bs, avgSendTime, packetlen, total_data):
        global channel
        self.nodeid = nodeid
        self.avgSendTime = avgSendTime
        self.bs = bs
        self.dist = distance[nodeid,:]
        self.mindist = np.amin(distance[nodeid,:])
        self.mindist_pos = int(np.where(distance[nodeid,:] == np.amin(distance[nodeid,:]))[0])
        #print('node %d' %nodeid, "dist: ", self.dist[0])
        self.buffer = total_data
        self.packetlen = packetlen
        #self.ch = int(random.choice(channel)) 
        self.packet = myPacket(self.nodeid, packetlen, self.dist)
        self.sent = 0 #INITIAL SENT PACKETS
        self.totalLost = 0 #INITIAL TOTAL LOST FOR PARTICULAR NODE
        self.totalColl = 0
        self.totalRec = 0
        self.totalProc = 0
        

class myPacket():
    def __init__(self, nodeid, packetlen, dist):
        #global experiment
        global Ptx
        global Prx
        #global gamma
        #global d0
        #global var
        global Lpl
        #global freq
        #global GL
        global c
        global distance
        global channel
        global frequency

        self.nodeid = nodeid
        self.txpow = Ptx
        self.sf = 12
        self.cr = 1 ##CODING RATE
        self.bw = 125
        # for experiment 3 find the best setting
        # OBS, some hardcoded values
        #Prx = self.txpow  ## zero path loss by default    
        #Prx = self.txpow + G_device + G_sat - Lpl_node

        # transmission range, needs update XXX
        self.transRange = 150
        self.pl = packetlen
        self.symTime = (2.0**self.sf)/self.bw
        self.arriveTime = 0
        self.rssi = Prx[nodeid,:]
        self.freq = int(random.choice(frequency)) 
        
        #self.ch = int(random.choice(channel))
        # frequencies: lower bound + number of 61 Hz steps
        #self.freq = 860000000 + random.randint(0,2622950)

        # for certain experiments override these and
        # choose some random frequences
        #print ("frequency" ,self.freq, "symTime ", self.symTime)
        #print ("bw", self.bw, "sf", self.sf, "cr", self.cr, "rssi", self.rssi)
        
        self.rectime = airtime(self.sf,self.cr,self.pl,self.bw) ##RECTIME IS THE RECEPTION TIME (ie AIRTIME)
        self.proptime = distance[nodeid,:]*(1/c)
        #print ("rectime node ", self.nodeid, "  ", self.rectime)
        #print ("Airtime for node {} is {} [seconds]".format(self.nodeid,self.rectime)) #from https://www.loratools.nl/#/airtime
        # denote if packet is collided
        self.collided = 0
        self.processed = 0
        self.lost = bool


def airtime(sf,cr,pl,bw):
    H = 0        # implicit header disabled (H=0) or not (H=1)
    DE = 0       # low data rate optimization enabled (=1) or not (=0)
    Npream = 8   # number of preamble symbol (12.25  from Utz paper)

    if bw == 125 and sf in [11, 12]:
        # low data rate optimization mandated for BW125 with SF11 and SF12
        DE = 1
    if sf == 6:
        # can only have implicit header with SF6
        H = 1

    Tsym = (2.0**sf)/bw
    Tpream = (Npream + 4.25)*Tsym
    #print ("PARAMS FOR TRANSMISION: sf", sf, " cr", cr, "pl", pl, "bw", bw)
    payloadSymbNB = 8 + max(math.ceil((8.0*pl-4.0*sf+28+16-20*H)/(4.0*(sf-2*DE)))*(cr+4),0)
    Tpayload = payloadSymbNB * Tsym
    return ((Tpream + Tpayload)/1000) ##IN SECS

def sendAck (env,node):
    global sf7,sf8,sf9,sf10,sf11,sf12    
    #if node.packet.lost == 0 and node.packet.collided ==0:
    packet_len = 5 ##IS A SMALL PACKET FOR ASKING
    ack_backoff = airtime(12,1,packet_len,125)
    #yield env.timeout(ack_backoff)
    print ("{:3.5f} || Request received from node {}. Sending ACK from Sat...".format(env.now,node.nodeid))
    rssi = node.packet.rssi[math.ceil(env.now)]
    #print ("{:3.5f} || RSSI for node {} is {} dB...".format(env.now,node.nodeid,rssi))
    ###SELECTING THE SF FOR DEVICE BASED ON RSSI
    if rssi > sf7[1]:
        #print ("----Select SF7")
        node.packet.sf = 7
    elif rssi > sf8[1]:
        #print ("----Select SF8")
        node.packet.sf = 8
    elif rssi > sf9[1]:
        #print ("----Select SF9")
        node.packet.sf = 9
    elif rssi > sf10[1]:
        #print ("----Select SF10")
        node.packet.sf = 10
    elif rssi > sf11[1]:
        #print ("----Select SF11")
        node.packet.sf = 11
    else:
        #print ("----Select S12")
        node.packet.sf = 12
    return ack_backoff 

def askSF (env,node):
    packet_len = 5 ##IS A SMALL PACKET FOR ASKING
    ask_backoff = airtime(12,1,packet_len,125)
    print ("{:3.5f} || Node {} sends asking packet...".format(env.now,node.nodeid))
    #yield env.timeout(wait)
    return ask_backoff
                             

def transmit(env,node):
    #while nodes[node.nodeid].buffer > 0.0:
    global wait_min
    global wait_max
    global back_off
    global beacon_time
    #firstPacket = 0
    while node.buffer > 0.0:
        yield env.timeout(node.packet.rectime + float(node.packet.proptime[math.ceil(env.now)])) ##GIVE TIME TO RECEIVE BEACON
        node.packet.sf = 12 ##ASSIGN SF=12 TO REACH SAT FOR ASKING PACKET              
        if node in packetsAtBS:
            print ("{:3.5f} || ERROR: packet is already in...".format(env.now))
        else:
            sensibility = sensi[node.packet.sf - 7, [125,250,500].index(node.packet.bw) + 1]
            if node.packet.rssi[math.ceil(env.now)] < sensibility: #HERE WE ARE CONSIDERING RSSI AT TIME ENV.NOW
                print ("{:3.5f} || Node {}: Can not reach beacon due Lpl".format(env.now,node.nodeid))
                wait =0 ##LETS WAIT FOR NEXT BEACON
                node.packet.lost = False
                trySend = False

            else:
                wait = random.uniform(0,back_off - node.packet.rectime - float(node.packet.proptime[math.ceil(env.now)])) ##TRIGGER BACK-OFF TIME
                yield env.timeout(wait)
                print ("{:3.5f} || Node {} begins to transmit a packet".format(env.now,node.nodeid))
                trySend = True
                #global packetlen
                #node.packet.pl = packetlen
                
                if node in packetsAtBS:
                    print ("{} || ERROR: packet is already in...".format(env.now))
                else:
                    
                    sensibility = sensi[node.packet.sf - 7, [125,250,500].index(node.packet.bw) + 1] ##CHECK IF askSF WILL ARRIVE
                    
                    if node.packet.rssi[math.ceil(env.now)] < sensibility: #HERE WE ARE CONSIDERING RSSI AT TIME ENV.NOW
                        print ("{:3.5f} || Node {}: Asking packet won't arrive to Sat".format(env.now,node.nodeid))
                        node.packet.lost = True ## LOST ONLY CONSIDERING Lpl
                    else:
                        node.packet.lost = False ## LOST ONLY CONSIDERING Lpl
                        #print ("{:3.5f} || Prx for node {} is {:3.2f} dB".format(env.now, node.nodeid, node.packet.rssi[math.ceil(env.now)]))
                        #print ("Prx for node",node.nodeid, "is: ",node.packet.rssi[math.ceil(env.now)],"at time",env.now)
                        #print ("{:3.5f} || Let's try if there are collisions...".format(env.now))
                        ask_backoff = askSF(env,node)
                        
                        if (checkcollision(node.packet)==1):
                            node.packet.collided = 1
                            #node.packet.processed = 0 ##FLUSH PROCESSED
                        else:
                            node.packet.collided = 0
                            #node.packet.processed = 0 ##FLUSH PROCESSED
                            print ("{:3.5f} || ...No Collision by now for asking packet!".format(env.now))
                        
                        packetsAtBS.append(node)
                        node.packet.addTime = env.now   
                        yield env.timeout(ask_backoff)
                        
                        if (node in packetsAtBS):
                            packetsAtBS.remove(node)
                        if node.packet.collided == 0 and node.packet.processed == 1:
                            #node.packet.processed = 0 ##FLUSH PROCESSED
                            ack_backoff = sendAck(env,node) 
                            yield env.timeout(ack_backoff)
                            #print ("----NODE SF",node.packet.sf)
                            node.sent = node.sent + 1
                            node.buffer = node.buffer - node.packetlen
                            #print ("{:3.5f} || Prx for node {} is {:3.2f} dB".format(env.now, node.nodeid, node.packet.rssi[math.ceil(env.now)]))
                            print ("{:3.5f} || SF {} is asigned to node {}...".format(env.now,node.packet.sf,node.nodeid))
                            print ("{:3.5f} || Let's try if there are collisions...".format(env.now))
                            if (checkcollision(node.packet)==1):
                                node.packet.collided = 1
                            else:
                                node.packet.collided = 0
                                print ("{:3.5f} || ...No Collision by now!".format(env.now))
                            packetsAtBS.append(node)
                            node.packet.addTime = env.now
                            yield env.timeout(node.packet.rectime)
                        else:
                            node.packet.processed = 0 ##FLUSH PROCESSED
                      
                    
                        
    
        if node.packet.lost:
            global nrLost
            nrLost += 1
            node.totalLost += 1 #ONLY DUE Lpl
        if node.packet.collided == 1:
            global nrCollisions
            nrCollisions = nrCollisions +1
            node.totalColl += 1
        
        
        #if node.packet.collided == 0 and node.packet.processed == 1 and not node.packet.lost and trySend:
        if node.packet.collided == 0 and node.packet.processed == 1 and not node.packet.lost and trySend:
            global nrReceived
            nrReceived = nrReceived + 1
            node.totalRec += 1
        
        
        if node.packet.processed == 1:
            global nrProcessed
            nrProcessed = nrProcessed + 1
            node.totalProc += 1
       
        # complete packet has been received by base station
        # Let's remove from Base Station
        if (node in packetsAtBS):
            packetsAtBS.remove(node)
            # reset the packet
        node.packet.collided = 0
        node.packet.processed = 0
        node.packet.lost = False
        
        #yield env.timeout(beacon_time-wait-node.packet.rectime)
        if trySend:
            yield env.timeout(beacon_time-wait-2*node.packet.rectime)
        else:
            yield env.timeout(beacon_time-wait-node.packet.rectime)
                      
                
def beacon (env):
    global beacon_time
    i = 0
    while True:
        if i == 0:
            yield env.timeout(0)           
        else:
            yield env.timeout(beacon_time)
        i=i+1
        print ("{:3.5f} || ***A new beacon has been sended from Satellite***".format(env.now))    
    

           
env.process(beacon(env)) ##BEACON SENDER

### THIS IS GOING TO CREATE NODES AND DO TRAMSMISIONS. IS THE MAIN PROGRAM ###
for i in range(nrNodes):
    node = myNode(i,bsId, avgSendTime, packetlen, total_data)
    nodes.append(node)
    env.process(transmit(env,node))
    
env.run(until=600)

sent = sum(n.sent for n in nodes)
print ("\nRESULTS....")
print ("Number of total sent packets (sent)",sent)
print ("Number of total collided packets (nrCollisions)",nrCollisions)
print ("Number of total lost packets (due Lpl) (nrLost)",nrLost)
print ("Number of total processed packets (nrProcessed)",nrProcessed)
print ("Number of total received packets (correct demodulation on gw) (nrReceived)",nrReceived)




   

#####PLOTTING STUFF ####
if plots_bar == 1:
    #### BAR PLOTS ####
    data = [sent,nrLost,nrCollisions,nrProcessed,nrReceived]
    plt.figure(figsize=(6, 4), dpi= 80, facecolor='w', edgecolor='k')
    #fig = plt.subplot(2,2,1)
    x = ["Sent","Lost(Lpl)","Collided","Processed","Received"]
    #y_pos = np.arange(len(x))
    plt.bar(x, data, color=["darkgreen","red","red","seagreen","springgreen"], align="center",width=0.25)
    plt.title("N° nodes: %i"%nrNodes)
    plt.xlabel('Packets')
    plt.ylabel('N° of packets')
    axes = plt.gca()
    #axes.set_xlim([0,xmax])
    axes.set_ylim([0,5000])
    #plt.ylim(100)
    #yint = range(min(sf_dist), math.ceil(max(sf_dist))+1)
    #plt.yticks(yint)
    plt.grid()
    plt.show()
    plt.savefig("bar.png")

### PARAMETERS FOR PLOTTING ###
###BEST NODE
node_best = int(np.where(Prx == np.amax(Prx))[0])
###WORST NODE
node_worst =int(np.where(Prx == np.amin(Prx))[0])

#SELECT NODE TO BE PLOTED
node = 10 #SELECT NODE TO PLOT IN 3D

if plots_nodes == 1:
    
    ###2D GRAPHIC  
    plt.figure(1)
    plt.title("Distance from node %i to Leo Sat" %node)
    plt.xlabel("Time/step Leo Sat passing by")
    plt.ylabel("Distance in [km]")
    plt.grid()
    plt.plot(distance[node,:])
    plt.show()
    ###2D GRAPHIC LINK BUDGET
    plt.figure(2)
    plt.subplot(1,3,1)
    plt.title("Link Budget for Node %i"%node)
    plt.xlabel("Time/step Leo Sat passing by")
    plt.ylabel("Prx in [dB]")
    plt.grid()
    plt.plot(Prx[node,:])

    plt.subplot(1,3,2)
    plt.title("Link Budget for Best node (node %i)"%node_best)
    plt.xlabel("Time/step Leo Sat passing by")
    plt.ylabel("Prx in [dB]")
    plt.grid()
    plt.plot(Prx[node_best,:])

    plt.subplot(1,3,3)
    plt.title("Link Budget for Worst Node (node %i)"%node_worst)
    plt.xlabel("Time/step Leo Sat passing by")
    plt.ylabel("Prx in [dB]")
    plt.grid()
    plt.plot(Prx[node_worst,:])
    plt.show()

    ###3D GRAPHIC 
    plt.figure(3)
    ax = plt.axes(projection="3d")
    ax.view_init(elev=-150, azim=-45)
    ax.set_title("Node %i vs Leo Sat Passing by" %node)
    ax.set_xlabel("AXIS X")
    ax.set_ylabel("AXIS Y")
    ax.set_zlabel("AXIS Z")
    ax.grid()
    ax.plot(leo_pos[:,0],leo_pos[:,1],leo_pos[:,2],'r',label="Leo Sat passing by")
    ax.scatter(sites_pos[node,0],sites_pos[node,1],sites_pos[node,2],c="green", label="Node %i position"%node)
    #ax.scatter (leo_pos[317,0],leo_pos[317,1],leo_pos[317,2],c="red") #FOR A PARTICULAR POINT ON LEO PASS
    ax.scatter (0,0,0,c="black",marker="X", label="Origin of coordinates")
    ax.legend()
    #ax.scatter (0,0,0,c="black",marker="$O$") #IF SPECIAL MARKER USED

    ###3D GRAPHIC CONSIDERING EARTH'S SURFACE
    plt.figure(4)
    ax = plt.axes(projection="3d")
    ax.view_init(elev=27, azim=-5) ##SUCH ANOTHER VIEW
    ax.set_title("Node %i vs Leo Sat Passing by" %node)
    ax.set_xlabel("AXIS X")
    ax.set_ylabel("AXIS Y")
    ax.set_zlabel("AXIS Z")
    ax.grid()
    ax.plot(leo_pos[:,0],leo_pos[:,1],leo_pos[:,2],'r',label="Leo Sat passing by")
    ax.scatter(sites_pos[node,0],sites_pos[node,1],sites_pos[node,2],c="green", label="Node %i position"%node)
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = 6378 * np.outer(np.cos(u), np.sin(v)) ##ECUATORIAL RADIO
    y = 6378 * np.outer(np.sin(u), np.sin(v))
    z = 6357 * np.outer(np.ones(np.size(u)), np.cos(v)) #POLAR RADIO
    ax.scatter (0,0,0,c="blue",marker="X", label="Earth's surface")
    ax.plot_surface(x, y, z,color="blue", alpha=0.1 )
    ax.legend()

# =============================================================================
# ###3D GRAPHIC CONSIDERING IMAGE OF EARTH'S SURFACE
# plt.figure(5)
# bm = PIL.Image.open('bluemarble_esc3.jpg')
# bm = np.array(bm)/255
# lons = np.linspace(-180, 180, bm.shape[1]) * np.pi/180 
# lats = np.linspace(-90, 90, bm.shape[0])[::-1] * np.pi/180 
# ax = plt.axes(projection="3d")
# ax.view_init(elev=27, azim=-5) ##SUCH ANOTHER VIEW
# ax.set_title("Node %i vs Leo Sat Passing by" %node)
# ax.set_xlabel("AXIS X")
# ax.set_ylabel("AXIS Y")
# ax.set_zlabel("AXIS Z")
# x = 6378 * np.outer(np.cos(lons), np.cos(lats)).T
# y = 6378 * np.outer(np.sin(lons), np.cos(lats)).T
# z = 6357 * np.outer(np.ones(np.size(lons)), np.sin(lats)).T
# ax.plot(leo_pos[:,0],leo_pos[:,1],leo_pos[:,2],'g',label="Leo Sat passing by")
# ax.scatter(sites_pos[node,0],sites_pos[node,1],sites_pos[node,2],c="red", label="Node %i position"%node)
# ax.scatter (0,0,0,c="green",marker="X", label="Earth's surface")
# ax.plot_surface(x, y, z, rstride=4, cstride=4, facecolors = bm, alpha=0.8)
# ax.legend()
# =============================================================================



