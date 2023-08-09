#!/usr/bin/python
# -- coding: utf-8 --
import sys

pos = open(sys.argv[1],"r")
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
# print(poslist)

NT = open(sys.argv[2],"r")
NTlist = []
for line2 in NT.readlines():
    line2=line2.replace("\n","")
    NTlist.append(line2)
# print(NTlist)

a = 0
filterlist = []
while a < len(poslist):
    if poslist[a][3] == "+" and (NTlist[a] == "CC" or NTlist[a] == "TT" or NTlist[a] == "CT" or NTlist[a] == "TC"):
        poslist[a].append(NTlist[a])
        filterlist.append(poslist[a])
        a = a + 1
    elif poslist[a][3] == "-" and (NTlist[a] == "CC" or NTlist[a] == "TT" or NTlist[a] == "CT" or NTlist[a] == "TC"):
        poslist[a].append(NTlist[a])
        filterlist.append(poslist[a])
        a = a + 1
    else:
        a = a + 1
# print(filterlist)
# print(len(filterlist))
for i in filterlist:
    print "\t".join(str(j) for j in i)+"\t"+i[3]
