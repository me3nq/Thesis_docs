#!/usr/bin/env python
# coding: utf-8

# In[1]:


'''Run this file in the directory with the .TXT files to be decoded. This will output
a CSV file for each channel with data for Analysis mode, and a single CSV file
for Verify mode, where the data for each event is stored in new row. Note, there
must be at least 2 events for the Analysis mode portion not to throw an error'''

import numpy as np
import struct
import os, glob
import scipy.io
import bitarray as b
import matplotlib.pyplot as plt
import csv
import re



#for the Analysis mode, there are 32*(4*8+4)=1152 bits per event
def Header_data(chunk):
    serial_num=int(chunk[5:14].to01(),2)  #with each event get the serial num, should be 000
    trig_num=int(chunk[14:].to01(),2)
        
    return trig_num,serial_num

def Trig_Time(chunk):
    word1=chunk[0:32]
    word2=chunk[32:]
    t=word2[8:]
    t.extend(word1[8:])
    time=int(t.to01(),2)
    return time

def single_chan_params(chunk):
    word1=chunk[0:32] #first 32 bits of each piece
    word2=chunk[32:64] #second ....
    word3=chunk[64:96] #third...
    word4=chunk[96:128] #...

    chan_num=int(word1[4:8].to01(),2)
    PulseInPed=word1[17]
    MaxPedDetect=word1[19]
    Pedestal=int(word1[20:].to01(),2)

    ovrflow=word2[5] #one or more samples in integral is overflow
    underflow=word2[6] #one or more samples in integral is underflow
    Area=int(word2[7:].to01(),2)

    pulse_too_long=word3[5] #True if the pulse extended beyond next trigger
    Tf=int(word3[6:19].to01(),2) #time relative to trigger when last sample <= TET
    Tr=int(word3[19:].to01(),2) #time relative to trigger when first sample >= TET

    Tp=int(word4[6:19].to01(),2) # time relative to trigger of Vp
    Vp=int(word4[20:].to01(),2) # Peak amplitude of pulse within integration window

    if ovrflow == True or underflow==True or pulse_too_long==True:
        flag=True #flag is true if there is a problem
    else:
        flag=False
    #chan_info=[chan_num,Area,Vp,Tp,Tf,Tr,flag]
    #chan_info=[chan_num,Area,Vp,Tp,Tf,Tr,Pedestal,MaxPedDetect,ovrflow,underflow,pulse_too_long]
    chan_info=[chan_num,Area,Vp,Tp,Tf,Tr,Pedestal,ovrflow,underflow] 
    return chan_info

def Pulse_Params(file):
    skipsize=1152 #size for analysis mode
    events=round(len(file)/skipsize)
    params=[]
    for i in range(events):
        params_chunk=file[i*skipsize+32*3:i*skipsize+32*3+32*32] #8 channels, each with four 32-bit words
        chan_info=[]
        for k in range(8):
            chan_info.append(single_chan_params(params_chunk[4*k*32:4*k*32+128]))
        params.append(chan_info)
    return params

def Trailer(file):
    skipsize=1152 #size for analysis mode
    events=round(len(file)/skipsize)
    end=[] #just gives the 32 bits at the end of each event
    for i in range(events):
        end.append(file[i*1152+1120:1152*(i+1)])
    return end

def datasaver(file_list,event_data,outvals,file_count):         
    if len(event_data)>0:
        with open(file_list[z][:-4]+'verify_'+str(file_count)+'.csv','w') as csvfile:
            datawriter=csv.writer(csvfile,lineterminator = '\n')                
            datawriter.writerow(['Channel Number','Raw Samples'])
            for i in range(len(event_data)):
                verifyoutput=event_data[i]
                CHAN=event_data[i][0]
                datawriter.writerow(verifyoutput)

    for k in range(8):
        CHAN=k
        field_headings=['Trigger time','Channel number', 'Area', 'Peak Height', 'Peak Time','End Time','Start Time', 
                        'Pedestal', 'overflow', 'underflow']
        if len(outvals)>0:
            outputdata=np.asarray(outvals)[:,CHAN,:]
            if len(np.nonzero(outputdata)[0])>0:
                with open(file_list[z][:-4]+'_'+str(file_count)+'Ch'+str(CHAN)+'.csv','w') as csvfile:
                    datawriter=csv.writer(csvfile,lineterminator = '\n')                
                    datawriter.writerow(field_headings)
                    datawriter.writerows(outputdata[outputdata[:,0]!=0])
                    

file_list=[]
for filename in glob.iglob(os.getcwd()+'\\**', recursive=True):
    if os.path.isfile(filename): # filter dirs
        if filename[-4:]=='.bin':
            file_list.append(filename)
            
for z in range(len(file_list)):
    check_start=[] #the first 4 bits should be 1001
    serial_num=[] #with each event get the serial num, should be 000
    trig_num=[]
    timestamp=[]
    outvals=[]
    event_data=[]
    data=b.bitarray()
    t=0
    chunksize=262144
    trigcount=0
    file_count=0
    
    fsize=os.path.getsize(file_list[z])
    file= open(file_list[z], mode='rb')    
    if fsize>chunksize:
        data.fromfile(file,chunksize)
    else:
        data.fromfile(file)

    event=False
    while len(data)>=32:
        if len(data)/8<round(chunksize/2) and file.tell() != fsize:
            print(file.tell())
            if fsize-file.tell()>chunksize:
                data.fromfile(file,chunksize)
            else:
                data.fromfile(file)
        indicator=data[0:4].to01()
        if indicator=='1001': #header word always there
            head=Header_data(data[0:32])
            t=Trig_Time(data[32:96])
            trigcount+=1
            #serial_num.append(head[1])
            trig_num.append(head[0])
            timestamp.append(t)
            data=data[96:]
        elif indicator=='1011': #window raw data - happens in verify mode 1 & 2? 
            chan_num=int(data[4:8].to01(),2)
            num_samples=int(data[19:32].to01(),2)
            data=data[32:]
            sampledata=[chan_num]
            for i in range(round(num_samples/2)):
                if len(data)>=32:
                    s1=int(data[7:19].to01(),2)
                    s2=int(data[20:32].to01(),2)
                    sampledata.append(s1)
                    sampledata.append(s2)
                    data=data[32:]
            event_data.append(sampledata)
        elif indicator=='1100': #pulse parameters happens in all modes
            event=True
            chandata=np.zeros((8,10))
            while data[0:4].to01()=='1100':
                chan_params=single_chan_params(data[0:128])
                chan_params.insert(0,t)
                chandata[chan_params[1]]=chan_params
                buffer=data[96:128]
                data=data[128:]
            outvals.append(chandata)
        elif indicator=='1000': #trailer always there
            check=int(data[1:32].to01())
            if check != 0:
                print('Problem - trailer not all zeros')
                break
            if event==True:
                event=False
            else:
                outvals.append(np.zeros((8,10)))
            if trigcount>=2.5e5 or len(data)<=32: #Output to csv every 250,000 triggers
                datasaver(file_list,event_data,outvals,file_count)
                trigcount=0
                trig_num=[]
                timestamp=[]
                outvals=[]
                event_data=[]
                file_count+=1
                
            data=data[32:]
        elif indicator=='0100':
            print('Error: Data saving problems, repeate entries. Skipping identical entries....')
            skipcount=0
            if data[0:32]==buffer:
                while data[0:32]==buffer:
                    skipcount=skipcount+4
                    data=data[32:]
                print('Skipped ',skipcount, ' bytes of data. Error hopefully mitigated.')
            else:
                data=data[32:]        
        else:
            print('Unknown Error, next bit indicator is: ', indicator)
            next_trailer=re.search('{0:32b}'.format(2**31),data[:10000].to01())
            if not next_trailer:
                data=data[10000:]
                print('At least 10,000 bits before next trailer...')
            else:
                data=data[next_trailer.end():]
                print('Trailer found after ', next_trailer.end(), ' bits, continuing file...')
    file.close()
    
    datasaver(file_list,event_data,outvals,file_count)
    

