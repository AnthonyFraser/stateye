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

eye 		= stateye()
extract 	= extractsignal()
channel		= touchstone()

############################################################
# default parameter than can be overwritten by command line 
# be aware than each object definition contains a further number of default parameters not visable at
# this hierarchy
############################################################

analysisMode		= 'channel compliance'	# 'channel compliance','tx compliance', 'rx compliance'

# stateye engine analysis parameters 
baudrate 		= 10.0e9
scrambled		= 1			# is analysis based on scrambled or 8b10b data
deemphasis		= 10.0			# which de-emphasis, in dB
dfetaps			= 0			# how many DFE taps
plotenable		= 1
BER			= -12
channel.filename	= 'halfmeter_A5A6B5B6.s4p'
extract.filename 	= '0m-prbs10-12_5ps-ag.csv'
############################################################
# check command line and overwrite parameters accordingly
############################################################

def usage():
    print "Arguments: -s <s4p filename> [-R N][-a N][-d N][-t N][-b N][-o <prefix>][-8]"
    print "-s <s4p filename>   input .s4p file (mandatory argument)"
    print "-e <filename>       measurement input file"
    print "-x N		       timestep of input file, if time column not included"	
    print "-a N                transmitter amplitude (in Vpp); default is 1 (ampN added to output filenames)"
    print "-d N                transmitter deemphasis (in dB); default is 0 (deempN added to output filenames)"
    print "-t N                receiver dfe taps; default is 3 (.dfeN added to output filename)"
    print "-c 8|R              coding - 8b10b (.8b10b added to output filenames) or random (.rand added to output filename)"
    print "-b N                baud rate (in bps); default is 6.0e9 (.baudN added to output filename)"
    print "-m N                analysis mode; 'rx compliance', 'tx compliance', 'channel compliance'"
    print "-R N		       random jitter"
    print "-z N                zero emphasis for extraction (80.0 for 10m channel, 1.0 for transmitter measurement)"
    print "-o <prefix>         prefix for result files"
    print "                    default is random.  8b10b simulations take longer."
    print "An optional argument only contribues to the output filename if used - default settings do not."

try :
	opts, args = getopt.getopt(sys.argv[1:], "a:c:d:e:s:t:m:b:q:i:R:x:z:o:", ["help"])
except getopt.GetoptError:
	usage()
    	sys.exit(2)

ion()
pauseMode = 1

outfilemodifier = ''
for opt, arg in opts:
    if opt in ("-h", "--help"):
        usage()
        sys.exit()                  
    elif opt == "-a":  # amplitude
        channel.txamp = float(arg)
        outfilemodifier += '_amp' + arg
	print "set transmit amplitude to %f"%channel.txamp
    elif opt == "-d":  # deemp
        deemphasis = float(arg)
        outfilemodifier += '_deemp' + arg
	print "set transmit de-emphasis to %f"%deemphasis
    elif opt == "-t":  # dfetaps
        dfetaps = int(arg)
        outfilemodifier += '_dfe' + arg
	print "set rx dfe taps to %f"%dfetaps
    elif opt == "-o":
	print "setting prefix to %s"%arg
	resultPrefix = arg    
    elif opt == "-c":  # coding scrambled/8b10b
        if arg == "8":
            scrambled = 0
            outfilemodifier += '_8b10b'
	    print "set scrambled"
        elif arg == "r":
            scrambled = 1
            outfilemodifier += '_rand'
	    print "set 8b10b"
    elif opt == "-b":  # baudrate
        baudrate = float(arg)
        outfilemodifier += '_baud' + arg
	print "set baudrate to %e"%baudrate
    elif opt == "-s":  # --s4p filename
        channel.filename = arg
	print "set channel filename to %s"%channel.filename
    elif opt == "-e":  
        extract.filename = arg
	print "set extract filename to %s"%extract.filename
    elif opt == "-q" : 	# quiet
        ioff()
	pauseMode = 0
	print "set graphics off"
    elif opt == "-m" :
	analysisMode = arg
	print "set analysis mode to %s"%analysisMode
    elif opt == "-x" :
	extract.timestep = float(arg)
	extract.timecol = []
	extract.sigcol  = 0
	print "set timestep to %e"%extract.timestep
    elif opt == "-R" :
 	eye.rj = float(arg)	   
    elif opt == "-z" :
	extract.alpha = float(arg)


############################################################
# main analysis routine
############################################################

if analysisMode == 'rx compliance' or analysisMode == 'tx compliance' :
	extract.samplePerBit	= eye.nomUI
	extract.resolution	= 1.0/baudrate/eye.nomUI	# this is in ps
	extract.findtaps()
	if analysisMode == 'rx compliance' :
		# defaults for jitter are used
		eye.loadStep(extract.step)
		eye.dj 			= extract.DJ
	if analysisMode == 'tx compliance' :
		eye.dj 			= extract.DJ

if analysisMode == 'channel compliance' or analysisMode == 'tx compliance' :
	# parameters for loading touchstone file
	channel.samplePerBit	= eye.nomUI
	channel.resolution	= 1.0/baudrate/eye.nomUI

	channel.loadFile()
	channel.map()
	channel.cascadeSimple()
	channel.extractTransfer()
	channel.addDC()
	if analysisMode == 'tx compliance' :
		channel.calcStep(extract.step, [])
	else :
		channel.calcStep([], [])
	eye.loadStep(channel.step)
	
fir 		= -(1.0 - 10**(-deemphasis / 20)) / 2.0

if scrambled :
        eye.create2TapFIR( [1.0+fir, fir], dfetaps)
else :
        eye.create8b10b_2TapFIR( [1.0+fir, fir], dfetaps)

eye.calcpdf()

############################################################
# figure outputs
############################################################

if plotenable :
	# pulse and step response
	figure()
	hold(0)
	x = []
	y = []
	y2 = []
	for cursor in arange(eye.startCursor,eye.lastCursor,0.01) :
		x += [cursor]
		y += [eye.pulse[eye.cursor2index(0,cursor)]]
		y2 += [eye.inputStep[eye.cursor2index(0,cursor+2)]]
		
	plot(x,y)
	hold(1)
	plot(x,y2)
	grid(1)
	xlabel('Time [UI]')
	ylabel('Amplitude [V]')
	title('Post channel Pulse and Step Response')
	#savefig(resultPrefix + outfilemodifier + '_pulse_step.png')
	
	# jitter statistical eye
	figure()
	pdf_pj_log = ( log10(transpose(eye.pdf_pj)+1.0e-17) )
	contour(eye.sweep, eye.binaxis , pdf_pj_log,arange(-15,0,0.5))
	#contourf(eye.sweep, eye.binaxis , pdf_pj_log,arange(-15,0,0.5))
	axis([-0.6,0.6,-1,1])
	title('Measurement function currently disabled')
#	maxamp = 0.0
#	minamp = eye.midBin
#		
#	for _a in transpose(pdf_pj_log)[2:-2] :
#		_maxamp = min(find(_a[eye.midBin:] > BER ))	
#		_minamp = max(find(_a[:eye.midBin] > BER ))	
#		if _maxamp > maxamp :
#			maxamp = _maxamp
#		if _minamp < minamp :
#			minamp = _minamp
#	amplitude = eye.binaxis[eye.midBin + maxamp]-eye.binaxis[minamp] 
#	print 'Amplitude is %0.3f'%(amplitude)
#	Jmin = max(find( pdf_pj_log[eye.midBin][:len(eye.sweep)/2] > BER ))
#	Jmax = min(find( pdf_pj_log[eye.midBin][len(eye.sweep)/2:] > BER )) + len(eye.sweep)/2
#
#	jitter = 1.0 - (eye.sweep[Jmax] - eye.sweep[Jmin])
#	print 'Jitter is %0.3f'%(jitter)
#
#	amplitude = eye.binaxis[eye.midBin + maxamp]-eye.binaxis[minamp] 
#
#	print 'Amplitude is %0.3f'%(amplitude)
#	print 'Jitter is %0.3f'%(jitter)
#
#	grid(1)
#	eyeLeft = find( max( pdf_pj_log[eye.midBin][:len(eye.sweep)/2] ) == pdf_pj_log[eye.midBin][:len(eye.sweep)/2] )[0]
#	eyeRght = find( max( pdf_pj_log[eye.midBin][len(eye.sweep)/2:] ) == pdf_pj_log[eye.midBin][len(eye.sweep)/2:] )[0] + len(eye.sweep)/2
#	axisMax = eye.binaxis[ max(find(transpose(pdf_pj_log)[0]>-17)) ]
#	axisMin = eye.binaxis[ min(find(transpose(pdf_pj_log)[0]>-17)) ]
#	axis([eye.sweep[eyeLeft],eye.sweep[eyeRght],axisMin*1.2,axisMax*1.2])
#	xlabel('Time [UI]')
#	ylabel('Amplitude [V]')
#	#title('Eye Opening %0.3fV, Jitter %0.3fUIpp\nTx=%0.3fmV, BER=10%d'%(amplitude, jitter,TxNorm,BER))
#	savefig(resultPrefix + outfilemodifier + '_stateye.png')
	
	print "prefix %s"%resultPrefix
	for f in range(8) :
		figure(f)
		print "saving %s.%d.png"%(resultPrefix,f)
		savefig('%s.%d.png'%(resultPrefix,f))

if pauseMode  :
	raw_input('Press anykey to continue')



