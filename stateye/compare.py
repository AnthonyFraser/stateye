############################################################
# main stateye analysis script
########################################################################################################################

import sys
import pdb
import getopt
sys.path += ['../v5_2']
import warnings
warnings.filterwarnings('ignore')

print "Tested with v5_5"

from pylab import *
#from time import *

from extractsignal 	import extractsignal
from touchstone 	import touchstone
from stateye 		import stateye

############################################################
# create objects required
############################################################

eye1 		= stateye()
eye2 		= stateye()
extract 	= extractsignal()
channel1	= touchstone()
channel2	= touchstone()

############################################################
# default parameter than can be overwritten by command line 
# be aware than each object definition contains a further number of default parameters not visable at
# this hierarchy
############################################################

analysisMode		= 'channel compliance'	# 'channel compliance','tx compliance', 'rx compliance'

# stateye engine analysis parameters 
baudrate 		= 6.0e9
scrambled		= 1			# is analysis based on scrambled or 8b10b data
deemphasis		= 0.0			# which de-emphasis, in dB
dfetaps			= 3			# how many DFE taps
plotenable		= 1
BER			= -12
channel1.filename	= 'halfmeter_A5A6B5B6.s4p'
channel2.filename	= 'halfmeter_A5A6A2A3NE.s4p'

############################################################
# main analysis routine
############################################################

channel1.samplePerBit	= eye1.nomUI
channel1.resolution	= 1.0/baudrate/eye1.nomUI
channel1.loadFile()
channel1.map()
channel1.cascadeSimple()
channel1.extractTransfer()
channel1.addDC()
channel1.calcStep([], [])
eye1.loadStep(channel1.step)
	
channel2.samplePerBit	= eye2.nomUI
channel2.resolution	= 1.0/baudrate/eye2.nomUI
channel2.loadFile()
channel2.map()
channel2.cascadeSimple()
channel2.extractTransfer()
channel2.addDC()
channel2.calcStep([], [])
eye2.loadStep(channel2.step)

############################################################
# figure outputs
############################################################

figure()
subplot(2,1,1)
hold(0)
x = []
y = []
y2 = []
for cursor in arange(eye1.startCursor,eye1.lastCursor,0.01) :
	x += [cursor]
	y += [eye1.pulse[eye1.cursor2index(0,cursor)]]
	y2 += [eye1.inputStep[eye1.cursor2index(0,cursor+2)]]
		
plot(x,y)
hold(1)
plot(x,y2)
grid(1)
xlabel('Time [UI]')
ylabel('Amplitude [V]')
	
subplot(2,1,2)
x = []
y = []
y2 = []
for cursor in arange(eye2.startCursor,eye2.lastCursor,0.01) :
	x += [cursor]
	y += [eye2.pulse[eye2.cursor2index(0,cursor)]]
		
plot(x,y)
hold(1)
grid(1)
xlabel('Time [UI]')
ylabel('Amplitude [V]')
	
