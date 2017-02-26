from pylab import *

hdis=[0.168009,0.119536,0.0736471,0.0532812,0.0341399,0.0235517] #maxlength
edis=[0.0232143,0.0184129,0.00918331,0.00766159,0.00531869,0.00279231] #relmaxerr
eocdis=[]
for i in arange(len(hdis)-1):
    eocdis.append(log(edis[i+1]/edis[i])/log(hdis[i+1]/hdis[i]))
print eocdis
print log(edis[-1]/edis[0])/log(hdis[-1]/hdis[0])
print sum(eocdis)/len(eocdis)

h=[0.131932,0.0943725,0.057498,0.041338,0.0264607,0.018636]
e=[0.00897851,0.00458192,0.00170132,0.000879471,0.000361231,0.000178942]
eoc=[]
for i in arange(len(h)-1):
    eoc.append(log(e[i+1]/e[i])/log(h[i+1]/h[i]))
print eoc
print log(e[-1]/e[0])/log(h[-1]/h[0])
print sum(eoc)/len(eoc)
