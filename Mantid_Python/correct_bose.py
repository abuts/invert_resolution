import numpy

file5K = 'MER31588_Ei150.00meV_One2One.nxspe'
file250K = 'MER31588_Ei150.00meV_One2One.nxspe'

if 'ws5K' in mtd:
    ws5K = mtd['ws5K']
else:
    ws5K = Load(Filename=file5K)
if 'ws250K' in mtd:
    ws250K = mtd['ws250K']
else:
    ws250K = Load(Filename=file250K)
#----------------------------------------------------    
en = ws5K.readX(0)
en_scale        = 0.5*abs(en[:-1]+en[1:])
Inv_bose5K     = 1-numpy.exp(-11.6044/5*en_scale)
Inv_bose250K  = 1-numpy.exp(-11.6044/250*en_scale)
#*Bose(5K)/Bose(250K)
bose_corr = Inv_bose250K/Inv_bose5K
corWS = ExtractSpectra(ws5K,StartWorkspaceIndex=0,EndWorkspaceIndex=0)
for i in xrange(0,bose_corr.size):
        corWS.dataY(0)[i]=bose_corr[i]
        
#tws =   ws250K*corWS
diff_ws = ws5K - ws250K*corWS


slw = plotSlice(diff_ws)
#slw.setXYDim(1,0)
slw.setColorScale(0,10000,0)
SaveNXSPE(diff_ws,'MER_Ei150_5vs250K_diffr.nxspe')