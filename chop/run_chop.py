#import PyChop2
import PyChop as pch

if __name__ == '__main__':
    ch = pch.ISISFermi.ISISFermi('MERLIN','G',500)
    tau = ch.getChopWidth(150)
    flux  = ch.getFlux(150)
    print 'tau = ',tau,' flux = ',flux
    