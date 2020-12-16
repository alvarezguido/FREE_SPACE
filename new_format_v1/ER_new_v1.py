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
import os
import datetime


####WE START BY USING SF=12 ADN BW=125 AND CR=1, FOR ALL NODES AND ALL TRANSMISIONS######
####WE ALSO CONSIDER SIMPLE CHECK, WHERE TWO PACKETS COLLIDE WHEN THEY ARRIVE AT: SAME TIME, SAME FREQUENCY AND SAME SF####
#nrNodes = 10 ##NUMBER OF NODES TO BE SIMULATED (IN ORDER FROM CSV FILE)
#multi_nodes = [1400,1000,500,250,100,50,25,10,5]
RANDOM_SEED = 6
random.seed(RANDOM_SEED) #RANDOM SEED IS FOR GENERATE ALWAYS THE SAME RANDOM NUMBERS (ie SAME RESULTS OF SIMULATION)

###PLOTS ##
plots_nodes = 0 ## FOR PLOT SIMULATION
plots_bar = 0 ##FOR PLOT BARS RESULTS OF SIMULATION

full_collision = False

###GLOBAL PARAMS ####
bsId = 1 ##ID OF BASE STATION (NOT USED)
#channel = [0,1,2] ##NOT USED BY NOW
#channel = [0]
avgSendTime = 3  ## NOT USED! --> A NODE SENDS A PACKET EVERY X SECS
packetlen = 20   ##NODES SEND PACKETS OF JUST 20 Bytes
total_data = 60 ##TOTAL DATA ON BUFFER, FOR EACH NODE (IT'S THE BUFFER O DATA BEFORE START SENDING)

beacon_time = 120 ###SAT SENDS BEACON EVERY CERTAIN TIME
back_off = beacon_time * 0.95 ###BACK OFF TIME FOR SEND A PACKET
packetsAtBS = [] ##USED FOR CHEK IF THERE ARE ALREADY PACKETS ON THE SATELLITE
c = 299792.458 ###SPEED LIGHT [km/s]
Ptx = 14
G_device = 0; ##ANTENNA GAIN FOR AN END-DEVICE
G_sat = 12;   ##ANTENNA GAIN FOR SATELLITE
nodes = [] ###EACH NODE WILL BE APPENDED TO THIS VARIABLE
freq =868e6 ##USED FOR PATH LOSS CALCULATION
frequency = [868100000, 868300000, 868500000] ##FROM LORAWAN REGIONAL PARAMETERS EU863-870 / EU868
#frequency = [868100000,868100000,868100000]

maxBSReceives = 500 ##MAX NUMBER OF PACKETS THAT BS (ie SATELLITE) CAN RECEIVE AT SAME TIME

nrLost = 0 ### TOTAL OF LOST PACKETS DUE Lpl
nrCollisions = 0 ##TOTAL OF COLLIDED PACKETS
nrProcessed = 0 ##TOTAL OF PROCESSED PACKETS
nrReceived = 0 ###TOTAL OF RECEIVED PACKETS
nrNoProcessed = 0 ##TOTAL OF INTRA-PACKETS NO PROCESSED
nrIntraTot = 0
nrLostMaxRec = 0
nrCollFullPacket = 0
nrSentIntra = 0 ##TOTAL OF SENT INTRA-PACKTES
nrReceivedIntra = 0 ##TOTAL OF RECEIVED INTRA-PACKETS

multi_nodes = [5,10,25,50,100,250,500,1000,1400]
#multi_nodes = [5,5,25,50,50,50,50,50,50]





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

elev = np.degrees(np.arcsin(600/distance))

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

def simulate_scenario (nrNodes):
    env = simpy.Environment()
    


    ##FOR CHECK HEADERS
    def checkcollision2(header,replica):
        col = 0 # flag needed since there might be several collisions for packet
        processing = 0
        #print ("PROCESSINGGGG ONEEE", processing)
        #print ("MAX RECEIVE IS: ", maxBSReceives)
        #print ("PACKETs AT BS",len(packetsAtBS))
        for i in range(0,len(packetsAtBS)):
            #print ("PACKETSSS ATT BSSSS PROCESSED",packetsAtBS[i].header.processed)
            if packetsAtBS[i].header.processed == 1:
                processing = processing + 1
        #print ("PROCESSINGGGG NUMBER", processing)
        if (processing >= maxBSReceives):
            #print ("{:3.5f} || Too much packets ({}) on BS.. node {} header replica {} is lost!".format(env.now, len(packetsAtBS),header.nodeid,replica))
            header.noProcessed +=1
        else:
            header.processed = 1
      
        if packetsAtBS:
            #print ("{:3.5f} || >> FOUND header overlap... node {}, others {}".format(env.now,header.nodeid,len(packetsAtBS)))
            for other in packetsAtBS:
                if other.nodeid != header.nodeid:
                   #print ("{:3.5f} || >> node {} header replica n° {} overlapped with node {}... Let's check channel...".format(env.now,header.nodeid,replica,other.nodeid))
                   # simple collision
                   #if frequencyCollision(packet, other.packet) and sfCollision(packet, other.packet):
                   #print ("REPLICAAAA:",replica)
                   if frequencyCollision2(header, other.header,replica) and sfCollision(nodes[header.nodeid].packet,nodes[other.nodeid].packet):             
                        header.collided += 1
                        other.header.collided += 1  # other also got lost, if it wasn't lost already. OTHER IS HEADER OR INTRAPACKET, ITS THE SAME
                        col = 1
                                                        
            return col
        return 0
    
    ###FOR CHECK INTRA-PACKETS (QUITE SIMILAR TO CHECKCOLLISION2)
    def checkcollision3(intraPacket,nrIntra):
        #print ("INTRA PACKET SUB CHHHH",intraPacket.subCh)
        col = 0 # flag needed since there might be several collisions for packet
        processing = 0
        #print ("MAX RECEIVE IS: ", maxBSReceives)
        for i in range(0,len(packetsAtBS)):
            if packetsAtBS[i].intraPacket.processed == 1:
                processing = processing + 1
        if (processing >= maxBSReceives):
            #print ("{:3.5f} || Too much packets ({}) on BS.. node {} intra-packet {} is lost!".format(env.now, len(packetsAtBS),intraPacket.nodeid,nrIntra))
            intraPacket.noProcessed +=1
        else:
            intraPacket.processed = 1
    
        if packetsAtBS:
            #print ("{:3.5f} || >> FOUND intra-packet overlap... node {}, others {}".format(env.now,intraPacket.nodeid,len(packetsAtBS)))
            for other in packetsAtBS:
                if other.nodeid != intraPacket.nodeid:
                   #print ("{:3.5f} || >> node {} intra-packet n° {} overlapped with node {}... Let's check channel...".format(env.now,intraPacket.nodeid,nrIntra,other.nodeid))
                   # simple collision
                   #if frequencyCollision(packet, other.packet) and sfCollision(packet, other.packet):
                   #print ("REPLICAAAA:",replica)
                   if frequencyCollision3(intraPacket, other.intraPacket,nrIntra) and sfCollision(nodes[intraPacket.nodeid].packet,nodes[other.nodeid].packet):
                        intraPacket.collided += 1
                        #print ("OTHER CLASSSS",other)
                        other.intraPacket.collided += 1  # other also got lost, if it wasn't lost already
                        col = 1
                                                        
            return col
        return 0
    
    
    def sfCollision(p1, p2):
        if p1.sf == p2.sf:
            print ("{:3.5f} || >> COLLISION! SF coll on node {} and node {} (ie same SF)...".format(env.now,p1.nodeid, p2.nodeid))
            # p2 may have been lost too, will be marked by other checks
            return True
        print ("{:3.5f} || >> No SF Collision!".format(env.now))
        return False
    
    ###FOR HEADER COLLISON
    def frequencyCollision2(p1,p2,replica):
        if (p1.ch == p2.ch):
            #print ("{:3.5f} || >> same channel for header on node {} and node {}.. Let's check sub-channels...".format(env.now,p1.nodeid,p2.nodeid))
            #if (p1.freqHopHeader[replica] == p2.freqHopHeader[replica]):
            if (p1.subCh == p2.subCh):
                #print ("{:3.5f} || >> same sub-channel for header on node {} and node {}".format(env.now,p1.nodeid,p2.nodeid))
                print ("{:3.5f} || >> Header {} from node {} collided!!!".format(env.now,replica,p1.nodeid))
                return True
            else:
                #print ("{:3.5f} || >> No sub-channel collision".format(env.now))
                pass
        else:
            #print ("{:3.5f} || >> No header channel collision..".format(env.now))
            return False
    
    ##FOR INTRA-PACKET COLLISION (QUITE SIMILAR TO FREQUENCYCOLLISION2)
    def frequencyCollision3(p1,p2,replica):
        if (p1.ch == p2.ch):
            #print ("{:3.5f} || >> same channel for intra-packet on node {} and node {}.. Let's check sub-channels...".format(env.now,p1.nodeid,p2.nodeid))
            #if (p1.freqHopHeader[replica] == p2.freqHopHeader[replica]):
            #print ("SUBCHANELLLLL",p1.subCh)
            if (p1.subCh == p2.subCh):
                #print ("{:3.5f} || >> same sub-channel for intra-packet on node {} and node {}".format(env.now,p1.nodeid,p2.nodeid))
                print ("{:3.5f} || >> Intra-packet {} from node {} collided!!!".format(env.now,replica,p1.nodeid))
                return True
            else:
                #print ("{:3.5f} || >> No sub-channel collision".format(env.now))
                pass
        else:
            #print ("{:3.5f} || >> No intra-packet channel collision..".format(env.now))
            return False
    
    
    
    class myNode():
        def __init__(self, nodeid, bs, avgSendTime, packetlen, total_data):
            global channel
            global DR
            self.dr = random.choice(DR)
            #carriers = list(range(280))
            #random.shuffle(carriers) #TO CHOOSE THE HOPPING JUMPS
            self.nodeid = nodeid
            self.avgSendTime = avgSendTime
            self.bs = bs
            self.dist = distance[nodeid,:]
            self.elev = elev[nodeid,:]
            self.mindist = np.amin(distance[nodeid,:])
            self.mindist_pos = int(np.where(distance[nodeid,:] == np.amin(distance[nodeid,:]))[0])
            #print('node %d' %nodeid, "dist: ", self.dist[0])
            self.buffer = total_data
            self.packetlen = packetlen
            self.ch = int(random.choice(channel)) 
            self.packet = myPacket(self.nodeid, packetlen, self.dist)
            #self.freqHop = carriers[0:35]
            self.sent = 0 #INITIAL SENT PACKETS
            self.totalLost = 0 #INITIAL TOTAL LOST FOR PARTICULAR NODE
            self.totalColl = 0
            self.totalRec = 0
            self.totalProc = 0
            if self.dr == "dr8":
                carriers = list(range(280))
                random.shuffle(carriers) #TO CHOOSE THE HOPPING JUMPS
                self.freqHop = carriers[0:35]
            elif self.dr == "dr9":
                carriers = list(range(280))
                random.shuffle(carriers) #TO CHOOSE THE HOPPING JUMPS
                self.freqHop = carriers[0:35]
            elif self.dr == "dr10":
                carriers = list(range(688))
                random.shuffle(carriers) #TO CHOOSE THE HOPPING JUMPS
                self.freqHop = carriers[0:86]
            elif self.dr == "dr11":
                carriers = list(range(688))
                random.shuffle(carriers) #TO CHOOSE THE HOPPING JUMPS
                self.freqHop = carriers[0:86]
            
            self.header = myHeader(self.nodeid,self.dist,self.ch,self.freqHop, self.dr)
            self.intraPacket = myIntraPacket(self.nodeid,self.dist,self.ch,self.freqHop,self.dr)
            
            
    class myHeader ():
        def __init__(self,nodeid,dist,ch,freqHop,dr):
            global Ptx
            global Prx
            global Lpl
            global c
            global distance
            global channel
            global frequency
            self.nodeid = nodeid
            self.txpow = Ptx
            self.transRange = 150
            self.arriveTime = 0
            self.rssi = Prx[nodeid,:]
            self.rectime = 0.233
            #self.rectime = 1.5
            self.proptime = distance[nodeid,:]*(1/c)
            self.collided = 0
            self.noCollided = 0
            self.processed = 0
            self.noProcessed = 0
            self.ch = ch
            self.lost = bool
            self.subCh = 0
            self.sentIntra = 0
            if dr == "dr8":
                self.freqHopHeader = freqHop[0:3]
            elif dr == "dr9":
                self.freqHopHeader = freqHop[0:2]
            elif dr == "dr10":
                self.freqHopHeader = freqHop[0:3]
            elif dr == "dr11":
                self.freqHopHeader = freqHop[0:2]
        
    
    class myIntraPacket ():
        def __init__(self,nodeid,dist,ch,freqHop,dr):
            global Ptx
            global Prx
            global Lpl
            global c
            global distance
            global channel
            global frequency
            self.nodeid = nodeid
            self.txpow = Ptx
            self.transRange = 150
            self.arriveTime = 0
            self.rssi = Prx[nodeid,:]
            self.freqHopIntraPacket = freqHop[3:]
            self.rectime = 50e-3
            #self.rectime = 3
            self.proptime = distance[nodeid,:]*(1/c)
            self.collided = 0
            self.noCollided = 0
            self.nrColl = 0
            self.processed = 0
            self.noProcessed = 0
            self.ch = ch
            self.lost = bool
            self.subCh = 0
            self.sentIntra = 0
            if dr == "dr8":
                self.freqHopIntraPacket = freqHop[3:]
            elif dr == "dr9":
                self.freqHopIntraPacket = freqHop[2:]
            elif dr == "dr10":
                self.freqHopIntraPacket = freqHop[3:]
            elif dr == "dr11":
                self.freqHopIntraPacket = freqHop[2:]
    
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
            #SF = [7,8,9,10,11,12]
            self.nodeid = nodeid
            self.txpow = Ptx
            #self.sf = random.choice(SF)
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
    
    
    def transmit(env,node):
        #while nodes[node.nodeid].buffer > 0.0:
        global wait_min
        global wait_max
        global back_off
        global beacon_time
        global logs
        while node.buffer > 0.0:
            yield env.timeout(node.packet.rectime + float(node.packet.proptime[math.ceil(env.now)])) ##GIVE TIME TO RECEIVE BEACON
                          
            if node in packetsAtBS:
                print ("{:3.5f} || ERROR: packet is already in...".format(env.now))
            else:
                sensibility = sensi[node.packet.sf - 7, [125,250,500].index(node.packet.bw) + 1]
                if node.packet.rssi[math.ceil(env.now)] < sensibility: #HERE WE ARE CONSIDERING RSSI AT TIME ENV.NOW
                    print ("{:3.5f} || Node {}: Can not reach beacon due Lpl".format(env.now,node.nodeid))
                    wait =0 ##LETS WAIT FOR NEXT BEACON
                    node.header.lost = False
                    node.intraPacket.lost = False
                    trySend = False
                    nIntraPackets =0
    
                else:
                    wait = random.uniform(0,back_off - node.packet.rectime - float(node.packet.proptime[math.ceil(env.now)])) ##TRIGGER BACK-OFF TIME
                    yield env.timeout(wait)
                    #print ("{:3.5f} || Node {} begins to transmit a packet".format(env.now,node.nodeid))
                    trySend = True
                    node.sent = node.sent + 1
                    node.buffer = node.buffer - node.packetlen
                    if node in packetsAtBS:
                        print ("{} || ERROR: packet is already in...".format(env.now))
                    else:
                        sensibility = sensi[node.packet.sf - 7, [125,250,500].index(node.packet.bw) + 1]
                        if node.packet.rssi[math.ceil(env.now)] < sensibility: #HERE WE ARE CONSIDERING RSSI AT TIME ENV.NOW
                            print ("{:3.5f} || Node {}: The Packet will be Lost due Lpl".format(env.now,node.nodeid))
                            node.header.lost = True ## LOST ONLY CONSIDERING Lpl
                            node.intraPacket.lost = True ## LOST ONLY CONSIDERING Lpl
                            #nIntraPackets =0
                        else:
                            node.header.lost = False ## LOST ONLY CONSIDERING Lpl
                            node.intraPacket.lost = False ## LOST ONLY CONSIDERING Lpl
                            #print ("{:3.5f} || Prx for node {} is {:3.2f} dB".format(env.now, node.nodeid, node.packet.rssi[math.ceil(env.now)]))
                            #print ("Prx for node",node.nodeid, "is: ",node.packet.rssi[math.ceil(env.now)],"at time",env.now)
                           
                            for i in range(len(node.header.freqHopHeader)):
                                ###print ("{:3.5f} || Sending Header replica {} node {}...".format(env.now,i,node.nodeid))
                                ###print ("{:3.5f} || Let's try if there are collisions...".format(env.now))
                                node.header.subCh = node.header.freqHopHeader[i]
                                #print ("SUBCHANELLLL: ",node.header.subCh)
                                node.header.sentIntra +=1;
                                if (checkcollision2(node.header,i)==1):
                                    #node.packet.collided = 1
                                    print ("---{:3.5f} || Collision for Header replica {} node {} !!!".format(env.now,i,node.nodeid))
                                    #node.packet.collided = 1
                                    #node.header.collided +=1 #ALREADY COUNTED ON FUNCTION                                
                                else:
                                    ###print ("{:3.5f} || ...No Collision for Header replica {} node {}!".format(env.now,i,node.nodeid))
                                    #node.packet.collided = 0
                                    node.header.noCollided = 1 ##ALMOST ONE HEADER IS OK, THEN HEADER IS OK
                                packetsAtBS.append(node)
                                node.packet.addTime = env.now
                            
                                yield env.timeout(node.header.rectime)
                                if (node in packetsAtBS):
                                    packetsAtBS.remove(node)
                            ##CALCULATE N OF INTRAPACKETS BASED ON PACKETLEN
                            payloadTime = airtime(12,1,node.packetlen,125)
                            nIntraPackets = math.ceil(payloadTime / 50e-3)
                            #print ("NUMBER OF INTRA PACKETSSSS",nIntraPackets)
                            
                            for j in range (nIntraPackets):
                                ###print ("{:3.5f} || Sending intra-packet {} of {} for node {}...".format(env.now,j,nIntraPackets-1,node.nodeid))
                                ###print ("{:3.5f} || Let's try if there are collisions...".format(env.now))
                                node.intraPacket.subCh = node.intraPacket.freqHopIntraPacket[j]
                                node.intraPacket.sentIntra +=1
                                #print ("INTRA-PACKT SUB CHANNELLLL", node.intraPacket.subCh)
                                if (checkcollision3(node.intraPacket,j)==1):
                                    print ("---{:3.5f} || Collision for intra-packet {} for node {} !!!".format(env.now,j,node.nodeid))
                                    #node.intraPacket.collided+=1 #ALREADY COUNTED ON FUNCTION
                                else:
                                    ###print ("{:3.5f} || ...No Collision for intra-packet {} for node {}!".format(env.now,j,node.nodeid))
                                    node.intraPacket.noCollided +=1
                                    pass
                                packetsAtBS.append(node)
                                node.packet.addTime = env.now
                                yield env.timeout(node.intraPacket.rectime)
                                if (node in packetsAtBS):
                                    packetsAtBS.remove(node)
                                #print ("INTRA-PACKET NO-PROCESEDDD",node.intraPacket.noProcessed)
                            
                            
    # =============================================================================
    #                         if (checkcollision(node.packet)==1):
    #                             node.packet.collided = 1
    #                         else:
    #                             node.packet.collided = 0
    #                             print ("{:3.5f} || ...No Collision by now!".format(env.now))
    #                         packetsAtBS.append(node)
    #                         node.packet.addTime = env.now
    #                         yield env.timeout(node.packet.rectime)
    # =============================================================================
            
            #print ("Intrapacket NOOO PROCESEDD",node.intraPacket.noProcessed)
    # =============================================================================
    #         global nrLostMaxRec
    #         nrLostMaxRec = nrLostMaxRec + node.header.noProcessed + node.intraPacket.noProcessed                       
    # =============================================================================
    # =============================================================================
    #         if node.header.noProcessed == 3 or node.intraPacket.noProcessed > (1/3)*nIntraPackets:
    #             global nrLostMaxRec
    #             nrLostMaxRec = nrLostMaxRec + 1
    # =============================================================================
            if trySend == 1:
                if node.header.lost or node.intraPacket.lost:
                    logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},PL,#{},#{},#{},#{}".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.dr,nIntraPackets,node.intraPacket.noCollided,len(node.header.freqHopHeader),node.header.noCollided))
                
                else:
                    if node.dr =="dr8" or node.dr=="dr10":
                        if node.header.collided == 3:
                            logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},PCh,#{},#{},#{},#{}".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.dr,nIntraPackets,node.intraPacket.noCollided,len(node.header.freqHopHeader),node.header.noCollided))
                        elif node.intraPacket.collided > (1/3)*nIntraPackets:
                            logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},PCp,#{},#{},#{},#{}".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.dr,nIntraPackets,node.intraPacket.noCollided,len(node.header.freqHopHeader),node.header.noCollided))
                        elif node.header.noProcessed == 3:
                            logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},NP,#{},#{},#{},#{}".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.dr,nIntraPackets,node.intraPacket.noCollided,len(node.header.freqHopHeader),node.header.noCollided))
                        elif node.intraPacket.noProcessed > (1/3)*nIntraPackets:
                            logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},NP,#{},#{},#{},#{}".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.dr,nIntraPackets,node.intraPacket.noCollided,len(node.header.freqHopHeader),node.header.noCollided))
                        else:
                            logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},PE,#{},#{},#{},#{}".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.dr,nIntraPackets,node.intraPacket.noCollided,len(node.header.freqHopHeader),node.header.noCollided))
                                   
                    elif node.dr=="dr9" or node.dr=="dr11":
                        if node.header.collided == 2:
                            logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},PCh,#{},#{},#{},#{}".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.dr,nIntraPackets,node.intraPacket.noCollided,len(node.header.freqHopHeader),node.header.noCollided))
                        elif node.intraPacket.collided > (2/3)*nIntraPackets:
                            logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},PCp,#{},#{},#{},#{}".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.dr,nIntraPackets,node.intraPacket.noCollided,len(node.header.freqHopHeader),node.header.noCollided))
                        elif node.header.noProcessed == 2:
                            logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},NP,#{},#{},#{},#{}".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.dr,nIntraPackets,node.intraPacket.noCollided,len(node.header.freqHopHeader),node.header.noCollided))
                        elif node.intraPacket.noProcessed > (2/3)*nIntraPackets:
                            logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},NP,#{},#{},#{},#{}".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.dr,nIntraPackets,node.intraPacket.noCollided,len(node.header.freqHopHeader),node.header.noCollided))
                        else:
                            logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},PE,#{},#{},#{},#{}".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.dr,nIntraPackets,node.intraPacket.noCollided,len(node.header.freqHopHeader),node.header.noCollided))
            
            if trySend == 1:
                
                global nrLost 
                if node.header.lost ==1 or node.intraPacket.lost == 1:
                    nrLost +=1
                    
                global nrSentIntra
                nrSentIntra = nrSentIntra + node.header.sentIntra + node.intraPacket.sentIntra
                    ####ALL PACKET IS LOST...
               
                global nrCollisions
                nrCollisions = nrCollisions + node.header.collided + node.intraPacket.collided
                
                global nrCollFullPacket
                if node.dr =="dr8" or node.dr=="dr10":
                    if node.header.collided == 3:
                        nrCollFullPacket +=1
                    elif node.intraPacket.collided > (1/3)*nIntraPackets:
                        nrCollFullPacket +=1
                elif node.dr=="dr9" or node.dr=="dr11":
                    if node.header.collided == 2:
                        nrCollFullPacket +=1
                    elif node.intraPacket.collided > (2/3)*nIntraPackets:
                        nrCollFullPacket +=1
               
                ##RECEIVED FULL PACKETS
                global nrReceived
                if node.dr == "dr8" or node.dr=="dr10":
                    if node.header.noProcessed <3 and node.header.collided <3 and node.intraPacket.noProcessed < (1/3)*nIntraPackets and node.intraPacket.collided < (1/3)*nIntraPackets:
                        nrReceived +=1
                elif node.dr == "dr9" or node.dr=="dr11":
                    if node.header.noProcessed <2 and node.header.collided <2 and node.intraPacket.noProcessed < (2/3)*nIntraPackets and node.intraPacket.collided < (2/3)*nIntraPackets:
                        nrReceived +=1
                
                ##NO PROCESSED PACKETS (too much intra-packets on BS)
                global nrNoProcessed
                nrNoProcessed = nrNoProcessed + node.header.noProcessed + node.intraPacket.noProcessed
                ##TOTAL OF RECEIVED INTRA-PACKETS
                global nrReceivedIntra
                nrReceivedIntra = nrReceivedIntra + node.header.noCollided + node.intraPacket.noCollided
            
            ##RESET
            node.header.collided = 0
            node.header.processed = 0
            node.header.noProcessed = 0
            node.header.lost = False
            node.header.noCollided =0
            node.intraPacket.nrColl = 0
            node.intraPacket.collided = 0
            node.intraPacket.processed = 0
            node.intraPacket.noProcessed = 0
            node.intraPacket.lost = False
            node.intraPacket.noCollided = 0
            node.header.sentIntra = 0
            node.intraPacket.sentIntra = 0
            if trySend:
                #print ("BEACON TIMEEE",beacon_time)
                #print ("WAITTT",wait)
                #print ("NODE HEADER TIME",node.header.rectime)
                #print ("ONE INTRA-PACKET TIMEE",node.intraPacket.rectime)
                #yield env.timeout(beacon_time-wait)
                yield env.timeout(beacon_time-wait-2*3*node.header.rectime-2*nIntraPackets*node.intraPacket.rectime)
            else:
                yield env.timeout(beacon_time-wait-3*node.header.rectime-nIntraPackets*node.intraPacket.rectime)
                #yield env.timeout(beacon_time-wait)
    # =============================================================================
    #         if node.packet.lost: ##LOST DUE LPL
    #             global nrLost
    #             nrLost += 1
    #             node.totalLost += 1 #ONLY DUE Lpl
    # =============================================================================
            #print ("INTRA PACKETSS COLLL",node.intraPacket.nrColl)
            #print ("N INTRA PACKETSSS",(1/3)*nIntraPackets)
            
            #node.totalColl += 1
            
            ###CHECKS FULL PACKET COLLISION
            #if trySend: 
            
            
    # =============================================================================
    #         if node.header.collided ==3 or node.intraPacket.collided > (1/3)*nIntraPackets:
    #             global nrCollFullPacket
    #             nrCollFullPacket +=1
    #         elif trySend:
    #             global nrReceived
    #             global nrIntraTot
    #             nrIntraTot = nrIntraTot + 3 + 27
    #             nrReceived = nrReceived + 1
    #             node.totalRec += 1
    # =============================================================================
            #if not node.header.noProcessed == 3 or not node.intraPacket.noProcessed > (1/3)*nIntraPackets:
    # =============================================================================
    #         if node.header.noProcessed < 3 or node.intraPacket.noProcessed < (1/3)*nIntraPackets:
    #             global nrProcessed
    #             nrProcessed = nrProcessed + 1
    #             node.totalProc += 1
    # =============================================================================
            
           
           
            
    # =============================================================================
    #         if node.header.noProcessed <3 or node.header.collided !=3:
    #             nrReceived +=1
    #         elif node.intraPacket.noProcessed < (1/3)*nIntraPackets or node.intraPacket.collided > (1/3)*nIntraPackets:
    #             nrReceived +=1
    # =============================================================================
                
            
            
         
            # complete packet has been received by base station
            # Let's remove from Base Station
    # =============================================================================
    #         if (node in packetsAtBS):
    #             packetsAtBS.remove(node)
    # =============================================================================
                # reset the packet
           
            
            #yield env.timeout(beacon_time-wait-node.packet.rectime)
            
                          
                    
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
        
    
    ###PLEASE UNCOMMENT FOR REGULAR FUNCIONALITY           
    env.process(beacon(env)) ##BEACON SENDER
    
    ### THIS IS GOING TO CREATE NODES AND DO TRAMSMISIONS. IS THE MAIN PROGRAM ###
    for i in range(nrNodes):
        node = myNode(i,bsId, avgSendTime, packetlen, total_data)
        nodes.append(node)
        env.process(transmit(env,node))
        
    env.run(until=600)
    
    sent = sum(n.sent for n in nodes)
    return ([sent,nrCollFullPacket,None,None,nrReceived],logs)

#multi_nodes = [1400]
DR = ["dr8","dr9","dr10","dr11"]
# =============================================================================
# ###SCENARIO 1 CHANNEL###
# channel = [0]
# 
# nodes = [] ###EACH NODE WILL BE APPENDED TO THIS VARIABLE
# nrLost = 0 ### TOTAL OF LOST PACKETS DUE Lpl
# nrCollisions = 0 ##TOTAL OF COLLIDED PACKETS
# nrProcessed = 0 ##TOTAL OF PROCESSED PACKETS
# nrReceived = 0 ###TOTAL OF RECEIVED PACKETS
# nrNoProcessed = 0 ##TOTAL OF INTRA-PACKETS NO PROCESSED
# nrIntraTot = 0
# nrLostMaxRec = 0
# nrCollFullPacket = 0
# nrSentIntra = 0 ##TOTAL OF SENT INTRA-PACKTES
# nrReceivedIntra = 0 ##TOTAL OF RECEIVED INTRA-PACKETS
# 
# i =0
# scenario_1ch = np.zeros((len(multi_nodes),5))
# results = []
# ## WHERE:
#     ## scenario_1ch[i,j]:
#         ## i --> the node i
#         ## j --> [sent, nrCollisions, nrLost, nrProcessed, nrReceived]
# 
# for nrNodes in multi_nodes:
#     print ("\n\n***NEW SCENARIO BEGINS***\n")
#     logs = []
#     results,logs = simulate_scenario(nrNodes)
#     print ("\n**Results for scenario {}**".format(i))
#     print ("Number of total sent full-packets (sent)",results[0])
#     print ("Number of total collided full-packets (nrCollFullPacket)",results[1])
#     print ("Number of total lost full-packets (due Lpl) (None)",results[2])
#     print ("Number of total processed full-packets (None)",results[3])
#     print ("Number of total received full-packets (correct demodulation on gw) (nrReceived)",results[4])
#     scenario_1ch[i,:] = results
#     fname = "./ER_1CH/" + str("ER_"+str(nrNodes)+"_1CH_"+str(maxBSReceives)) + ".csv"
#     with open(fname,"w") as myfile:
#         myfile.write("\n".join(logs))
#     myfile.close()
#     i=i+1
#     nodes = [] ###EACH NODE WILL BE APPENDED TO THIS VARIABLE
#     nrLost = 0 ### TOTAL OF LOST PACKETS DUE Lpl
#     nrCollisions = 0 ##TOTAL OF COLLIDED PACKETS
#     nrProcessed = 0 ##TOTAL OF PROCESSED PACKETS
#     nrReceived = 0 ###TOTAL OF RECEIVED PACKETS
#     nrNoProcessed = 0 ##TOTAL OF INTRA-PACKETS NO PROCESSED
#     nrIntraTot = 0
#     nrLostMaxRec = 0
#     nrCollFullPacket = 0
#     nrSentIntra = 0 ##TOTAL OF SENT INTRA-PACKTES
#     nrReceivedIntra = 0 ##TOTAL OF RECEIVED INTRA-PACKETS
# =============================================================================

# =============================================================================
# fname = str("ER") + ".txt"
# #print (fname)
# date = "ER_1CH"+"\n"+str(datetime.datetime.now())+"\n"
# nods = str(multi_nodes[0])+","+str(multi_nodes[1])+","+str(multi_nodes[2])+","+str(multi_nodes[3])+","+str(multi_nodes[4])\
#      +","+str(multi_nodes[5])+","+str(multi_nodes[6])+","+str(multi_nodes[7])+","+str(multi_nodes[8])+"\n"
# header = "sent,nrCollFullPacket,None,None,nrReceived"
# #if os.path.isfile(fname):
# res = "\n"+ str(scenario_1ch[0,0]) + "," + str(scenario_1ch[0,1]) + "," + str(scenario_1ch[0,2]) + "," + str(scenario_1ch[0,3]) + "," + str(scenario_1ch[0,4]) \
#     + "\n"+ str(scenario_1ch[1,0]) + "," + str(scenario_1ch[1,1]) + "," + str(scenario_1ch[1,2]) + "," + str(scenario_1ch[1,3]) + "," + str(scenario_1ch[1,4]) \
#     + "\n"+ str(scenario_1ch[2,0]) + "," + str(scenario_1ch[2,1]) + "," + str(scenario_1ch[2,2]) + "," + str(scenario_1ch[2,3]) + "," + str(scenario_1ch[2,4]) \
#     + "\n"+ str(scenario_1ch[3,0]) + "," + str(scenario_1ch[3,1]) + "," + str(scenario_1ch[3,2]) + "," + str(scenario_1ch[3,3]) + "," + str(scenario_1ch[3,4]) \
#     + "\n"+ str(scenario_1ch[4,0]) + "," + str(scenario_1ch[4,1]) + "," + str(scenario_1ch[4,2]) + "," + str(scenario_1ch[4,3]) + "," + str(scenario_1ch[4,4]) \
#     + "\n"+ str(scenario_1ch[5,0]) + "," + str(scenario_1ch[5,1]) + "," + str(scenario_1ch[5,2]) + "," + str(scenario_1ch[5,3]) + "," + str(scenario_1ch[5,4]) \
#     + "\n"+ str(scenario_1ch[6,0]) + "," + str(scenario_1ch[6,1]) + "," + str(scenario_1ch[6,2]) + "," + str(scenario_1ch[6,3]) + "," + str(scenario_1ch[6,4]) \
#     + "\n"+ str(scenario_1ch[7,0]) + "," + str(scenario_1ch[7,1]) + "," + str(scenario_1ch[7,2]) + "," + str(scenario_1ch[7,3]) + "," + str(scenario_1ch[7,4]) \
#     + "\n"+ str(scenario_1ch[8,0]) + "," + str(scenario_1ch[8,1]) + "," + str(scenario_1ch[8,2]) + "," + str(scenario_1ch[8,3]) + "," + str(scenario_1ch[8,4]) 
# 
# newres ="\n\n"+ date+"nodes\n"+nods+header+res
# #print (newres)
# with open(fname, "a") as myfile:
#     myfile.write(newres)
# myfile.close()       
# =============================================================================


###SCENARIO 3 CHANNELS###
channel = [0,1,2]
nodes = [] ###EACH NODE WILL BE APPENDED TO THIS VARIABLE
nrLost = 0 ### TOTAL OF LOST PACKETS DUE Lpl
nrCollisions = 0 ##TOTAL OF COLLIDED PACKETS
nrProcessed = 0 ##TOTAL OF PROCESSED PACKETS
nrReceived = 0 ###TOTAL OF RECEIVED PACKETS
nrNoProcessed = 0 ##TOTAL OF INTRA-PACKETS NO PROCESSED
nrIntraTot = 0
nrLostMaxRec = 0
nrCollFullPacket = 0
nrSentIntra = 0 ##TOTAL OF SENT INTRA-PACKTES
nrReceivedIntra = 0 ##TOTAL OF RECEIVED INTRA-PACKETS
i =0
scenario_3ch = np.zeros((len(multi_nodes),5))
results = []

for nrNodes in multi_nodes:
    print ("\n\n***NEW SCENARIO BEGINS***\n")
    logs = []
    results,logs = simulate_scenario(nrNodes)
    print ("\n**Results for scenario {}**".format(i))
    print ("Number of total sent full-packets (sent)",results[0])
    print ("Number of total collided full-packets (nrCollFullPacket)",results[1])
    print ("Number of total lost full-packets (due Lpl) (None)",results[2])
    print ("Number of total processed full-packets (None)",results[3])
    print ("Number of total received full-packets (correct demodulation on gw) (nrReceived)",results[4])
    scenario_3ch[i,:] = results
    fname = "./ER_3CH/" + str("ER_"+str(nrNodes)+"_3CH_"+str(maxBSReceives)) + ".csv"
    with open(fname,"w") as myfile:
        myfile.write("\n".join(logs))
    myfile.close()
    i=i+1
    nodes = [] ###EACH NODE WILL BE APPENDED TO THIS VARIABLE
    nrLost = 0 ### TOTAL OF LOST PACKETS DUE Lpl
    nrCollisions = 0 ##TOTAL OF COLLIDED PACKETS
    nrProcessed = 0 ##TOTAL OF PROCESSED PACKETS
    nrReceived = 0 ###TOTAL OF RECEIVED PACKETS
    nrNoProcessed = 0 ##TOTAL OF INTRA-PACKETS NO PROCESSED
    nrIntraTot = 0
    nrLostMaxRec = 0
    nrCollFullPacket = 0
    nrSentIntra = 0 ##TOTAL OF SENT INTRA-PACKTES
    nrReceivedIntra = 0 ##TOTAL OF RECEIVED INTRA-PACKETS

# =============================================================================
# fname = str("ER") + ".txt"
# #print (fname)
# date = "ER_3CH"+"\n"+str(datetime.datetime.now())+"\n"
# nods = str(multi_nodes[0])+","+str(multi_nodes[1])+","+str(multi_nodes[2])+","+str(multi_nodes[3])+","+str(multi_nodes[4])\
#      +","+str(multi_nodes[5])+","+str(multi_nodes[6])+","+str(multi_nodes[7])+","+str(multi_nodes[8])+"\n"
# header = "sent,nrCollFullPacket,None,None,nrReceived"
# #if os.path.isfile(fname):
# res = "\n"+ str(scenario_3ch[0,0]) + "," + str(scenario_3ch[0,1]) + "," + str(scenario_3ch[0,2]) + "," + str(scenario_3ch[0,3]) + "," + str(scenario_3ch[0,4]) \
#     + "\n"+ str(scenario_3ch[1,0]) + "," + str(scenario_3ch[1,1]) + "," + str(scenario_3ch[1,2]) + "," + str(scenario_3ch[1,3]) + "," + str(scenario_3ch[1,4]) \
#     + "\n"+ str(scenario_3ch[2,0]) + "," + str(scenario_3ch[2,1]) + "," + str(scenario_3ch[2,2]) + "," + str(scenario_3ch[2,3]) + "," + str(scenario_3ch[2,4]) \
#     + "\n"+ str(scenario_3ch[3,0]) + "," + str(scenario_3ch[3,1]) + "," + str(scenario_3ch[3,2]) + "," + str(scenario_3ch[3,3]) + "," + str(scenario_3ch[3,4]) \
#     + "\n"+ str(scenario_3ch[4,0]) + "," + str(scenario_3ch[4,1]) + "," + str(scenario_3ch[4,2]) + "," + str(scenario_3ch[4,3]) + "," + str(scenario_3ch[4,4]) \
#     + "\n"+ str(scenario_3ch[5,0]) + "," + str(scenario_3ch[5,1]) + "," + str(scenario_3ch[5,2]) + "," + str(scenario_3ch[5,3]) + "," + str(scenario_3ch[5,4]) \
#     + "\n"+ str(scenario_3ch[6,0]) + "," + str(scenario_3ch[6,1]) + "," + str(scenario_3ch[6,2]) + "," + str(scenario_3ch[6,3]) + "," + str(scenario_3ch[6,4]) \
#     + "\n"+ str(scenario_3ch[7,0]) + "," + str(scenario_3ch[7,1]) + "," + str(scenario_3ch[7,2]) + "," + str(scenario_3ch[7,3]) + "," + str(scenario_3ch[7,4]) \
#     + "\n"+ str(scenario_3ch[8,0]) + "," + str(scenario_3ch[8,1]) + "," + str(scenario_3ch[8,2]) + "," + str(scenario_3ch[8,3]) + "," + str(scenario_3ch[8,4]) 
#         
# #else:
#     ##NRTRANSMISSIONS IS TOTAL OF SENT PACKAGES, OR SENT VARIABLE
#  #   res = "#randomseed, collType, nrNodes, DataSize, nrTransmissions, nrCollisions, nrlost, nrLostError, nrnoack, nracklost, CollectionTIme, DER1, DER2, OverallEnergy, nodefair1, nodefair2, sfdistribution, slotlengths, framelengths, Guards\n" + str(sys.argv[4]) + ", " + str(full_collision) + ", " + str(nrNodes) + ", " + str(datasize) + ", " + str(sent) + ", "  + str(nrCollisions) + ", "  + str(nrLost) + ", "  + str(nrLostError) + ", " +str(nrNoACK) + ", " +str(nrACKLost) + ", " + str(env.now)+ ", " + str(der1) + ", " + str(der2)  + ", " + str(energy) + ", "  + str(nodefair1) + ", "  + str(nodefair2) + ", "  + str(SFdistribution) + ", "  + str(Slotlengths) + ", "  + str(Framelengths) + ", "  + str(Guards)
# #newres=re.sub('[^#a-zA-Z0-9 \n\.]','',res)
# newres ="\n\n"+ date+"nodes\n"+nods+header+res
# #print (newres)
# with open(fname, "a") as myfile:
#     myfile.write(newres)
# myfile.close()
# =============================================================================

