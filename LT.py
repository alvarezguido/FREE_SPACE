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
G_sat = 5;   ##ANTENNA GAIN FOR SATELLITE
nodes = [] ###EACH NODE WILL BE APPENDED TO THIS VARIABLE
freq =868e6 ##USED FOR PATH LOSS CALCULATION
frequency = [868100000, 868300000, 868500000] ##FROM LORAWAN REGIONAL PARAMETERS EU863-870 / EU868
#frequency = [868100000,868100000,868100000]
maxBSReceives = 8 ##MAX NUMBER OF PACKETS THAT BS (ie SATELLITE) CAN RECEIVE AT SAME TIME
nrLost = 0 ### TOTAL OF LOST PACKETS DUE Lpl
nrCollisions = 0 ##TOTAL OF COLLIDED PACKETS
nrProcessed = 0 ##TOTAL OF PROCESSED PACKETS
nrReceived = 0 ###TOTAL OF RECEIVED PACKETS

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
Prx = Ptx + G_sat + G_device -20*np.log10(distance*1000) - 20*np.log10(freq) + 147.55 #DISTANCE IS CONVERTED TO METERS
## WHERE:
    ## Prx[i,j]:
        ## i --> the node i
        ## j --> the step time in sat pass 
c7=0
c8=0
c9=0
c10=0
c11=0
c12=0
def simulate_scenario (nrNodes):
    env = simpy.Environment()
    
    def checkcollision(packet):
        col = 0 # flag needed since there might be several collisions for packet
        processing = 0
        #print ("MAX RECEIVE IS: ", maxBSReceives)
        for i in range(0,len(packetsAtBS)):
            if packetsAtBS[i].packet.processed == 1:
                processing = processing + 1
        if (processing > maxBSReceives):
            print ("{:3.5f} || Too much packets on Base Sattion.. Packet will be lost!", len(packetsAtBS))
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
    
    def selectSF (env, node):
        global sf7,sf8,sf9,sf10,sf11,sf12 
        global c7,c8,c9,c10,c11,c12
        rssi = node.packet.rssi[math.ceil(env.now)]
        #print ("{:3.5f} || RSSI for node {} is {} dB...".format(env.now,node.nodeid,rssi))
        if rssi > sf7[1]:
            #print ("----Select SF7")
            node.packet.sf = 7
            c7+=1
        elif rssi > sf8[1]:
            #print ("----Select SF8")
            node.packet.sf = 8
            c8+=1
        elif rssi > sf9[1]:
            #print ("----Select SF9")
            node.packet.sf = 9
            c9+=1
        elif rssi > sf10[1]:
            #print ("----Select SF10")
            node.packet.sf = 10
            c10+=1
        elif rssi > sf11[1]:
            #print ("----Select SF11")
            node.packet.sf = 11
            c11+=1
        else:
            #print ("----Select S12")
            node.packet.sf = 12
            c12+=1
        return c7,c8,c9,c10,c11,c12

    
    def transmit(env,node):
        #while nodes[node.nodeid].buffer > 0.0:
        global wait_min
        global wait_max
        global back_off
        global beacon_time
        while node.buffer > 0.0:
            node.packet.sf = 12
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
                    c7,c8,c9,c10,c11,c12 = selectSF(env,node)
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
                            print ("{:3.5f} || Prx for node {} is {:3.2f} dB".format(env.now, node.nodeid, node.packet.rssi[math.ceil(env.now)]))
                            #print ("Prx for node",node.nodeid, "is: ",node.packet.rssi[math.ceil(env.now)],"at time",env.now)
                            print ("{:3.5f} || Let's try if there are collisions...".format(env.now))
                            if (checkcollision(node.packet)==1):
                                node.packet.collided = 1
                            else:
                                node.packet.collided = 0
                                print ("{:3.5f} || ...No Collision by now!".format(env.now))
                            packetsAtBS.append(node)
                            node.packet.addTime = env.now
                            yield env.timeout(node.packet.rectime)
            
            if node.packet.lost:
                global nrLost
                nrLost += 1
                node.totalLost += 1 #ONLY DUE Lpl
            if node.packet.collided == 1:
                global nrCollisions
                nrCollisions = nrCollisions +1
                node.totalColl += 1
            
            
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
            node.packet.sf = 12
            
            #yield env.timeout(beacon_time-wait-node.packet.rectime)
            if trySend:
                yield env.timeout(beacon_time-wait-2*node.packet.rectime)
            else:
                yield env.timeout(beacon_time-wait-node.packet.rectime)
        #print(c7,c8,c9,c10,c11,c12)
                          
                    
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
    
    return ([sent,nrCollisions,nrLost,nrProcessed,nrReceived])


#multi_nodes = [5,10,15,20]
multi_nodes = [100]
#############################################################
###SCENARIO 1 CHANNEL###
frequency = [868100000] #1 CH

nodes = [] ###EACH NODE WILL BE APPENDED TO THIS VARIABLE
nrLost = 0 ### TOTAL OF LOST PACKETS DUE Lpl
nrCollisions = 0 ##TOTAL OF COLLIDED PACKETS
nrProcessed = 0 ##TOTAL OF PROCESSED PACKETS
nrReceived = 0 ###TOTAL OF RECEIVED PACKETS

i =0
scenario_1ch = np.zeros((len(multi_nodes),5))
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
    print ("Number of total collided packets (nrCollisions)",results[1])
    print ("Number of total lost packets (due Lpl) (nrLost)",results[2])
    print ("Number of total processed packets (nrProcessed)",results[3])
    print ("Number of total received packets (correct demodulation on gw) (nrReceived)",results[4])
    scenario_1ch[i,:] = results
    i=i+1
   #nodes = [] ###EACH NODE WILL BE APPENDED TO THIS VARIABLE
    nrLost = 0 ### TOTAL OF LOST PACKETS DUE Lpl
    nrCollisions = 0 ##TOTAL OF COLLIDED PACKETS
    nrProcessed = 0 ##TOTAL OF PROCESSED PACKETS
    nrReceived = 0 ###TOTAL OF RECEIVED PACKETS

###EXPORT scenario_1ch

# save experiment data into a dat file that can be read by e.g. gnuplot
# name of file would be:  exp0.dat for experiment 0
# =============================================================================
# fname = str("LT") + ".csv"
# #print (fname)
# date = "LT_1CH"+"\n"+str(datetime.datetime.now())+"\n"
# nods = str(multi_nodes[0])+","+str(multi_nodes[1])+","+str(multi_nodes[2])+","+str(multi_nodes[3])+","+str(multi_nodes[4])\
#      +","+str(multi_nodes[5])+","+str(multi_nodes[6])+","+str(multi_nodes[7])+","+str(multi_nodes[8])+"\n"
# header = "sent,nrCollisions,nrLost,nrProcessed,nrReceived"
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

#received_1ch = scenario_1ch[:,4] / scenario_1ch[:,0]
#collided_1ch = scenario_1ch[:,1] / scenario_1ch[:,0]
#############################################################

#############################################################

# =============================================================================
# ###SCENARIO 3 CHANNELS###
# frequency = [868100000, 868300000, 868500000] ##FROM LORAWAN REGIONAL PARAMETERS EU863-870 / EU868
# 
# nodes = [] ###EACH NODE WILL BE APPENDED TO THIS VARIABLE
# nrLost = 0 ### TOTAL OF LOST PACKETS DUE Lpl
# nrCollisions = 0 ##TOTAL OF COLLIDED PACKETS
# nrProcessed = 0 ##TOTAL OF PROCESSED PACKETS
# nrReceived = 0 ###TOTAL OF RECEIVED PACKETS
# 
# i =0
# scenario_3ch = np.zeros((len(multi_nodes),5))
# results = []
# ## WHERE:
#     ## scenario_3ch[i,j]:
#         ## i --> the node i
#         ## j --> [sent, nrCollisions, nrLost, nrProcessed, nrReceived]
# for nrNodes in multi_nodes:
#     print ("\n\n***NEW SCENARIO BEGINS***\n")
#     results = simulate_scenario(nrNodes)
#     print ("\n**Results for scenario {}**".format(i))
#     print ("Number of total sent packets (sent)",results[0])
#     print ("Number of total collided packets (nrCollisions)",results[1])
#     print ("Number of total lost packets (due Lpl) (nrLost)",results[2])
#     print ("Number of total processed packets (nrProcessed)",results[3])
#     print ("Number of total received packets (correct demodulation on gw) (nrReceived)",results[4])
#     scenario_3ch[i,:] = results
#     i=i+1
#     nodes = [] ###EACH NODE WILL BE APPENDED TO THIS VARIABLE
#     nrLost = 0 ### TOTAL OF LOST PACKETS DUE Lpl
#     nrCollisions = 0 ##TOTAL OF COLLIDED PACKETS
#     nrProcessed = 0 ##TOTAL OF PROCESSED PACKETS
#     nrReceived = 0 ###TOTAL OF RECEIVED PACKETS
# 
# fname = str("LT") + ".csv"
# #print (fname)
# date = "LT_3CH"+"\n"+str(datetime.datetime.now())+"\n"
# nods = str(multi_nodes[0])+","+str(multi_nodes[1])+","+str(multi_nodes[2])+","+str(multi_nodes[3])+","+str(multi_nodes[4])\
#      +","+str(multi_nodes[5])+","+str(multi_nodes[6])+","+str(multi_nodes[7])+","+str(multi_nodes[8])+"\n"
# header = "sent,nrCollisions,nrLost,nrProcessed,nrReceived"
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
#received_3ch = scenario_3ch[:,4] / scenario_3ch[:,0]
#collided_3ch = scenario_3ch[:,1] / scenario_3ch[:,0]
