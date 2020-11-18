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
nrNodes = 150 ##NUMBER OF NODES TO BE SIMULATED (IN ORDER FROM CSV FILE)
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
packetsAtBS = [] ##USED FOR CHEK IF THERE ARE ALREADY PACKETS ON THE SATELLITE
c = 299792.458 ###SPEED LIGHT [km/s]
Ptx = 14
G_device = 0; ##ANTENNA GAIN FOR AN END-DEVICE
G_sat = 12;   ##ANTENNA GAIN FOR SATELLITE
nodes = [] ###EACH NODE WILL BE APPENDED TO THIS VARIABLE
freq =868e6 ##USED FOR PATH LOSS CALCULATION
frequency = [868100000, 868300000, 868500000] ##FROM LORAWAN REGIONAL PARAMETERS EU863-870 / EU868
#frequency = [868100000,868100000,868100000]

maxBSReceives = 8 ##MAX NUMBER OF PACKETS THAT BS (ie SATELLITE) CAN RECEIVE AT SAME TIME

nrLost = 0 ### TOTAL OF LOST PACKETS DUE Lpl
nrCollisions = 0 ##TOTAL OF COLLIDED PACKETS
nrProcessed = 0 ##TOTAL OF PROCESSED PACKETS
nrReceived = 0 ###TOTAL OF RECEIVED PACKETS
nrNoProcessed = 0 ##TOTAL OF INTRA-PACKETS NO PROCESSED
nrIntraTot = 0
nrLostMaxRec = 0
nrCollFullPacket = 0


multi_nodes = [5,10,25,50,100,250,500,1000,1400]




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

Prx = Ptx + G_sat + G_device - Lpl #DISTANCE IS CONVERTED TO METERS
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
            print ("{:3.5f} || Too much packets ({}) on BS.. node {} header replica {} is lost!".format(env.now, len(packetsAtBS),header.nodeid,replica))
            header.noProcessed +=1
        else:
            header.processed = 1
      
        if packetsAtBS:
            print ("{:3.5f} || >> FOUND header overlap... node {}, others {}".format(env.now,header.nodeid,len(packetsAtBS)))
            for other in packetsAtBS:
                if other.nodeid != header.nodeid:
                   print ("{:3.5f} || >> node {} header replica n° {} overlapped with node {}... Let's check channel...".format(env.now,header.nodeid,replica,other.nodeid))
                   # simple collision
                   #if frequencyCollision(packet, other.packet) and sfCollision(packet, other.packet):
                   #print ("REPLICAAAA:",replica)
                   if frequencyCollision2(header, other.header,replica):               
                        header.collided = 1
                        other.header.collided = 1  # other also got lost, if it wasn't lost already
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
            print ("{:3.5f} || Too much packets ({}) on BS.. node {} intra-packet {} is lost!".format(env.now, len(packetsAtBS),intraPacket.nodeid,nrIntra))
            intraPacket.noProcessed +=1
        else:
            intraPacket.processed = 1
    
        if packetsAtBS:
            print ("{:3.5f} || >> FOUND intra-packet overlap... node {}, others {}".format(env.now,intraPacket.nodeid,len(packetsAtBS)))
            for other in packetsAtBS:
                if other.nodeid != intraPacket.nodeid:
                   print ("{:3.5f} || >> node {} intra-packet n° {} overlapped with node {}... Let's check channel...".format(env.now,intraPacket.nodeid,nrIntra,other.nodeid))
                   # simple collision
                   #if frequencyCollision(packet, other.packet) and sfCollision(packet, other.packet):
                   #print ("REPLICAAAA:",replica)
                   if frequencyCollision3(intraPacket, other.intraPacket,nrIntra):
                        intraPacket.collided = 1
                        #print ("OTHER CLASSSS",other)
                        other.intraPacket.collided = 1  # other also got lost, if it wasn't lost already
                        col = 1                                 
            return col
        return 0
   
        ###FOR HEADER COLLISON
    def frequencyCollision2(p1,p2,replica):
        if (p1.ch == p2.ch):
            print ("{:3.5f} || >> same channel for header on node {} and node {}.. Let's check sub-channels...".format(env.now,p1.nodeid,p2.nodeid))
            #if (p1.freqHopHeader[replica] == p2.freqHopHeader[replica]):
            if (p1.subCh == p2.subCh):
                print ("{:3.5f} || >> same sub-channel for header on node {} and node {}".format(env.now,p1.nodeid,p2.nodeid))
                print ("{:3.5f} || >> Header {} from node {} collided!!!".format(env.now,replica,p1.nodeid))
                return True
            else:
                print ("{:3.5f} || >> No sub-channel collision".format(env.now))
        else:
            print ("{:3.5f} || >> No header channel collision..".format(env.now))
            return False
    
    ##FOR INTRA-PACKET COLLISION (QUITE SIMILAR TO FREQUENCYCOLLISION2)
    def frequencyCollision3(p1,p2,replica):
        if (p1.ch == p2.ch):
            print ("{:3.5f} || >> same channel for intra-packet on node {} and node {}.. Let's check sub-channels...".format(env.now,p1.nodeid,p2.nodeid))
            #if (p1.freqHopHeader[replica] == p2.freqHopHeader[replica]):
            #print ("SUBCHANELLLLL",p1.subCh)
            if (p1.subCh == p2.subCh):
                print ("{:3.5f} || >> same sub-channel for intra-packet on node {} and node {}".format(env.now,p1.nodeid,p2.nodeid))
                print ("{:3.5f} || >> Intra-packet {} from node {} collided!!!".format(env.now,replica,p1.nodeid))
                return True
            else:
                print ("{:3.5f} || >> No sub-channel collision".format(env.now))
        else:
            print ("{:3.5f} || >> No intra-packet channel collision..".format(env.now))
            return False
    
    class myNode():
        def __init__(self, nodeid, bs, avgSendTime, packetlen, total_data):
            global channel
            carriers = list(range(280))
            random.shuffle(carriers) #TO CHOOSE THE HOPPING JUMPS
            self.nodeid = nodeid
            self.avgSendTime = avgSendTime
            self.bs = bs
            self.dist = distance[nodeid,:]
            self.mindist = np.amin(distance[nodeid,:])
            self.mindist_pos = int(np.where(distance[nodeid,:] == np.amin(distance[nodeid,:]))[0])
            #print('node %d' %nodeid, "dist: ", self.dist[0])
            self.buffer = total_data
            self.packetlen = packetlen
            self.ch = int(random.choice(channel)) 
            self.packet = myPacket(self.nodeid, packetlen, self.dist)
            self.freqHop = carriers[0:35]
            self.header = myHeader(self.nodeid,self.dist,self.ch,self.freqHop)
            self.intraPacket = myIntraPacket(self.nodeid,self.dist,self.ch,self.freqHop)
            self.sent = 0 #INITIAL SENT PACKETS
            self.totalLost = 0 #INITIAL TOTAL LOST FOR PARTICULAR NODE
            self.totalColl = 0
            self.totalRec = 0
            self.totalProc = 0
    
    class myHeader ():
        def __init__(self,nodeid,dist,ch,freqHop):
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
            self.freqHopHeader = freqHop[0:3]
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
    
    class myIntraPacket ():
        def __init__(self,nodeid,dist,ch,freqHop):
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
    
    def transmit(env,node):
        #while nodes[node.nodeid].buffer > 0.0:
        global wait_min
        global wait_max
        global back_off
        global beacon_time
        while node.buffer > 0.0:
            yield env.timeout(node.packet.rectime + float(node.packet.proptime[math.ceil(env.now)])) ##GIVE TIME TO RECEIVE BEACON
                          
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
                    node.sent = node.sent + 1
                    node.buffer = node.buffer - node.packetlen
                    if node in packetsAtBS:
                        print ("{} || ERROR: packet is already in...".format(env.now))
                    else:
                        sensibility = sensi[node.packet.sf - 7, [125,250,500].index(node.packet.bw) + 1]
                        if node.packet.rssi[math.ceil(env.now)] < sensibility: #HERE WE ARE CONSIDERING RSSI AT TIME ENV.NOW
                            print ("{:3.5f} || Node {}: The Packet will be Lost due Lpl".format(env.now,node.nodeid))
                            node.packet.lost = True ## LOST ONLY CONSIDERING Lpl
                        else:
                            node.packet.lost = False ## LOST ONLY CONSIDERING Lpl
                            #print ("{:3.5f} || Prx for node {} is {:3.2f} dB".format(env.now, node.nodeid, node.packet.rssi[math.ceil(env.now)]))
                            #print ("Prx for node",node.nodeid, "is: ",node.packet.rssi[math.ceil(env.now)],"at time",env.now)
                           
                            for i in range(3):
                                ###print ("{:3.5f} || Sending Header replica {} node {}...".format(env.now,i,node.nodeid))
                                ###print ("{:3.5f} || Let's try if there are collisions...".format(env.now))
                                node.header.subCh = node.header.freqHopHeader[i]
                                #print ("SUBCHANELLLL: ",node.header.subCh)
                                if (checkcollision2(node.header,i)==1):
                                    #node.packet.collided = 1
                                    print ("---{:3.5f} || Collision for Header replica {} node {} !!!".format(env.now,i,node.nodeid))
                                    #node.packet.collided = 1
                                    node.header.collided +=1
                                else:
                                    ###print ("{:3.5f} || ...No Collision for Header replica {} node {}!".format(env.now,i,node.nodeid))
                                    #node.packet.collided = 0
                                    node.header.noCollided = 1 ##ALMOST ONE HEADER IS OK, THEN HEADER IS OK
                                packetsAtBS.append(node)
                                node.packet.addTime = env.now
                            
                                yield env.timeout(node.header.rectime)
                                if (node in packetsAtBS):
                                    packetsAtBS.remove(node)
                            payloadTime = airtime(12,1,node.packetlen,125)
                            nIntraPackets = math.ceil(payloadTime / 50e-3)
                            #print ("NUMBER OF INTRA PACKETSSSS",nIntraPackets)
                            
                            for j in range (nIntraPackets):
                                ###print ("{:3.5f} || Sending intra-packet {} of {} for node {}...".format(env.now,j,nIntraPackets-1,node.nodeid))
                                ###print ("{:3.5f} || Let's try if there are collisions...".format(env.now))
                                node.intraPacket.subCh = node.intraPacket.freqHopIntraPacket[j]
                                #print ("INTRA-PACKT SUB CHANNELLLL", node.intraPacket.subCh)
                                if (checkcollision3(node.intraPacket,j)==1):
                                    print ("---{:3.5f} || Collision for intra-packet {} for node {} !!!".format(env.now,j,node.nodeid))
                                else:
                                    ###print ("{:3.5f} || ...No Collision for intra-packet {} for node {}!".format(env.now,j,node.nodeid))
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
            if node.header.noProcessed == 3 or node.intraPacket.noProcessed > (1/3)*nIntraPackets:
                global nrLostMaxRec
                nrLostMaxRec = nrLostMaxRec + 1
            if node.packet.lost: ##LOST DUE LPL
                global nrLost
                nrLost += 1
                node.totalLost += 1 #ONLY DUE Lpl
            #print ("INTRA PACKETSS COLLL",node.intraPacket.nrColl)
            #print ("N INTRA PACKETSSS",(1/3)*nIntraPackets)
                ####ALL PACKET IS LOST...
            global nrCollisions
            nrCollisions = nrCollisions + node.header.collided + node.intraPacket.collided
            #node.totalColl += 1
            
            if node.header.collided ==3 or node.intraPacket.collided > (1/3)*nIntraPackets:
                global nrCollFullPacket
                nrCollFullPacket +=1
            elif trySend:
                global nrReceived
                global nrIntraTot
                nrIntraTot = nrIntraTot + 3 + 27
                nrReceived = nrReceived + 1
                node.totalRec += 1
            #if not node.header.noProcessed == 3 or not node.intraPacket.noProcessed > (1/3)*nIntraPackets:
            if node.header.noProcessed < 3 or node.intraPacket.noProcessed < (1/3)*nIntraPackets:
                global nrProcessed
                nrProcessed = nrProcessed + 1
                node.totalProc += 1
            
            ##NO PROCESSED PACKETS (too much intra-packets on BS)
            global nrNoProcessed
            nrNoProcessed = nrNoProcessed + node.header.noProcessed + node.intraPacket.noProcessed
            # complete packet has been received by base station
            # Let's remove from Base Station
    # =============================================================================
    #         if (node in packetsAtBS):
    #             packetsAtBS.remove(node)
    # =============================================================================
                # reset the packet
            node.header.collided = 0
            node.header.processed = 0
            node.header.noProcessed = 0
            node.header.lost = False
            node.intraPacket.nrColl = 0
            node.intraPacket.collided = 0
            node.intraPacket.processed = 0
            node.intraPacket.noProcessed = 0
            node.intraPacket.lost = False
            
            #yield env.timeout(beacon_time-wait-node.packet.rectime)
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
    
    for i in range(nrNodes):
        node = myNode(i,bsId, avgSendTime, packetlen, total_data)
        nodes.append(node)
        env.process(transmit(env,node))
    env.run(until=600)
    sent = sum(n.sent for n in nodes)
    totalLost = nrLost+nrLostMaxRec+nrCollFullPacket
    nrProcessed = sent - totalLost
    nrIntraProcessed = nrIntraTot-nrNoProcessed
    #return ([sent,nrCollisions,totalLost,nrReceived])
    return ([sent,nrCollisions,nrLost,totalLost,nrIntraProcessed,nrProcessed])




multi_nodes = [5,10,15,20]
    
#############################################################
###SCENARIO 1 CHANNEL###
channel = [0] 

nodes = [] ###EACH NODE WILL BE APPENDED TO THIS VARIABLE
nrLost = 0 ### TOTAL OF LOST PACKETS DUE Lpl
nrCollisions = 0 ##TOTAL OF COLLIDED PACKETS
nrProcessed = 0 ##TOTAL OF PROCESSED PACKETS
nrReceived = 0 ###TOTAL OF RECEIVED PACKETS
nrNoProcessed = 0 ##TOTAL OF INTRA-PACKETS NO PROCESSED
nrIntraTot = 0
nrLostMaxRec = 0
nrCollFullPacket = 0
totalLost = 0
nrIntraProcessed = 0

i =0
scenario_1ch = np.zeros((len(multi_nodes),6))
results = []
## WHERE:
    ## scenario_1ch[i,j]:
        ## i --> the node i
        ## j --> [sent, nrCollisions, nrLost, nrProcessed, nrReceived]

for nrNodes in multi_nodes:
    print ("\n\n***NEW SCENARIO BEGINS***\n")
    results = simulate_scenario(nrNodes)
    print ("\n**Results for scenario {}**".format(i))
    print ("Number of total sent packets (sent)",results[0])
    print ("Number of total collided intra-packets (nrCollisions)",results[1])
    print ("Number of total lost packets (due Lpl,MaxRecBS,coll full-packet) (totalLost)",results[3])
    print ("Number of total processed packets (nrProcessed)",results[5])
    #print ("Number of total received packets (correct demodulation on gw) (nrReceived)",results[4])
    scenario_1ch[i,:] = results
    i=i+1

processed_1ch = scenario_1ch[:,5] / scenario_1ch[:,0]
totalLost_1ch = scenario_1ch[:,3] / scenario_1ch[:,0]
#############################################################

#############################################################
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
totalLost = 0
nrIntraProcessed = 0

i =0
scenario_3ch = np.zeros((len(multi_nodes),6))
results = []
## WHERE:
    ## scenario_3ch[i,j]:
        ## i --> the node i
        ## j --> [sent, nrCollisions, nrLost, nrProcessed, nrReceived]
for nrNodes in multi_nodes:
    print ("\n\n***NEW SCENARIO BEGINS***\n")
    results = simulate_scenario(nrNodes)
    print ("\n**Results for scenario {}**".format(i))
    print ("Number of total sent packets (sent)",results[0])
    print ("Number of total collided intra-packets (nrCollisions)",results[1])
    print ("Number of total lost packets (due Lpl,MaxRecBS,coll full-packet) (totalLost)",results[3])
    print ("Number of total processed packets (nrProcessed)",results[5])
    #print ("Number of total received packets (correct demodulation on gw) (nrReceived)",results[4])
    scenario_3ch[i,:] = results
    i=i+1


processed_3ch = scenario_1ch[:,5] / scenario_1ch[:,0]
totalLost_3ch = scenario_1ch[:,3] / scenario_1ch[:,0]
#############################################################

plot_curves =1 #FOR PLOT CURVES AMONG SEVERAL SCENARIOS
if plot_curves == 1:
    plt.figure(figsize=(6, 4), dpi= 80, facecolor='w', edgecolor='k')
    #plt.figure(figsize=(16, 10), dpi= 80, facecolor='w', edgecolor='k')
    x = multi_nodes
    plt.title("Simulation Results")
    plt.xlabel("N° of nodes")
    plt.ylabel("count")
    plt.grid()
    
    plt.plot(x,processed_1ch,'c',label="nrProcessed/sent - 1 ch",marker="^")
    plt.plot(x,totalLost_1ch,'c',label="nrtotLost/sent - 1 ch",marker="v")
    
    #plt.figure(2)
    plt.plot(x,processed_3ch,'b',label="nrProcessed/sent - 3 ch",marker="^")
    plt.plot(x,totalLost_3ch,'b',label="nrtotLost/sent - 3 ch",marker="v")
    plt.xlim(0,800)
    #plt.xlim(0,200)
    plt.show()
    plt.legend()
    plt.savefig("results_simulation.png")
    
    

# =============================================================================
# if plots_bar == 1:
#     #### BAR PLOTS ####
#     data = [sent,nrLost,nrCollisions,nrProcessed,nrReceived]
#     plt.figure(figsize=(6, 4), dpi= 80, facecolor='w', edgecolor='k')
#     #fig = plt.subplot(2,2,1)
#     x = ["Sent","Lost(Lpl)","Collided","Processed","Received"]
#     #y_pos = np.arange(len(x))
#     plt.bar(x, data, color=["darkgreen","red","red","seagreen","springgreen"], align="center",width=0.25)
#     plt.title("N° nodes: %i"%nrNodes)
#     plt.xlabel('Packets')
#     plt.ylabel('N° of packets')
#     axes = plt.gca()
#     #axes.set_xlim([0,xmax])
#     axes.set_ylim([0,5000])
#     #plt.ylim(100)
#     #yint = range(min(sf_dist), math.ceil(max(sf_dist))+1)
#     #plt.yticks(yint)
#     plt.grid()
#     plt.show()
#     plt.savefig("bar.png")
# =============================================================================





