#!/usr/bin/python 
import re
import sys
import copy

if len(sys.argv) > 1:
    f=open(sys.argv[1])
else:
    sys.exit()

space=re.compile(r'(^ *$)')
comm1=re.compile(r'(^\*--)')
comm2=re.compile(r'(^\*==)')
comm3=re.compile(r'(^\!--)')
cont=re.compile(r'(^     [&,$,x])')
statement=re.compile(r'(^      [A-Z][A-Z]*)')
line=f.readline()
i=0
a=[]
b=[]
while line:
    a.append(line)
    i = i + 1
    line=f.readline()

sys.stdout=open(sys.argv[1],"w")

for i in range(len(a)):
    if cont.search(a[i]):
        sa=a[i-1][:-1]
        while len(sa) < 80:
            sa=sa+' '

        sa=sa+'&'
        b.pop(i-1)
        b.insert(i-1,sa)
        b.append(a[i])
    elif comm1.search(a[i]):
        ga=comm1.sub('!!$--',a[i])
        b.append(ga)        
    elif comm2.search(a[i]):
        ga=comm2.sub('!!$--',a[i])
        b.append(ga)        
    elif comm3.search(a[i]):
        ga=comm3.sub('!!$--',a[i])
        b.append(ga)        
    else:
        b.append(a[i][:-1])

for i in range(len(b)):
    print "%-s\n" % (b[i]),
