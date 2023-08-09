#!/usr/bin/python
# -- coding: utf-8 --


import sys

pos = open(sys.argv[1],"r")

#pos=open("GSM2585795-No_MMS.sacCer3.pos.info",'r')
poslist = []
for line1 in pos.readlines():
    line1=line1.replace(':',',').replace("\n","").replace('-',',').replace('(',',').replace(')','').replace('>','')
    # print(line)
    poslinelist=line1.split(',')
    new_poslinelist = [int(n) if n.isdigit() else n for n in poslinelist]
    if len(new_poslinelist)==5 :
        new_poslinelist.pop(-1)
        new_poslinelist.pop(-1)
        new_poslinelist.append("-")
    poslist.append(new_poslinelist)

#print(poslist[:30])

NT = open(sys.argv[2],"r")
#NT=open("GSM2585795-No_MMS.sacCer3.NT.info",'r')
NTlist = []
for line2 in NT.readlines():
    line2=line2.replace("\n","")
    NTlist.append(line2)
    
#print(NTlist[:30])

#a = 0
filterlist = []

for i in range(len(poslist)):
    if NTlist[i-1]=="A" or NTlist[i-1]=="G":
        poslist[i-1].append(NTlist[i-1])
        filterlist.append(poslist[i-1])

for i in filterlist:
    print ("\t".join(str(j) for j in i)+"\t"+i[3])

