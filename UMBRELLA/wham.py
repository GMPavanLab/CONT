import math
import sys

# arguments read from command line
# name of input file
FILENAME_ = sys.argv[1]
# number of BIAS
NBIAS_ = int(sys.argv[2])
# temperature
KBT_ = float(sys.argv[3])
nwham = int(sys.argv[4])
deltaout = int(sys.argv[5])
# optional arguments
if len(sys.argv) > 6:
   nwham_i=int(sys.argv[6])
   filez=sys.argv[7]
else :
   nwham_i=0
# default parameters for WHAM

# convergence thresold
THRES_ = 1.0e-10


def get_wham_weights(nbias, nframes, bias_ts, deltaout, nwham_i, nwham, filez=None, thres=THRES_ ):
    # find minimum bias
    fconv = open('convergence.dat', 'w+')
    min_bias = min(bias_ts)

    # initialize weights
    w = []
    for i in range(0, nframes): w.append(1.0)
    # offset and exponential of the bias
    expv = []
    for i in range(0, len(bias_ts)): expv.append(math.exp((-bias_ts[i]+min_bias)/KBT_))

    Z = []
    if nwham_i==0:
        # initialize Z
        for j in range(0, nbias): Z.append(1.0)
    else:
        for lines in open(filez, "r").readlines():
            Z.append(float(lines))
             
    # WHAM iterations
    for iii in range(nwham_i, nwham):
        # store Z
        Z_old = Z[:]
        # recompute weights
        norm = 0.0
        for i in range(0, len(w)):
            ew = 0.0
            for j in range(0, len(Z)): ew += expv[i*len(Z)+j] / Z[j]
            w[i] = 1.0 / (ew+0.000000000001)
            norm += w[i]

        # normalize weights
        for i in range(0, len(w)):
            w[i] /= norm
        #print intermediate output 
        if iii%deltaout == 0:
            outfile="weights_"+str(iii)+".dat"
            fout = open(outfile, 'w+')
            for i in range(0, len(w)): fout.write("%32.30lf\n" % w[i])
            fout.close()
            #print intermediate z
            outz="z_"+str(iii)+".dat"
            fz = open(outz, 'w+')
            for i in range(0, len(Z)): fz.write("%32.30lf\n" % Z[i])
            fz.close()

        # recompute Z
        for j in range(0, len(Z)): Z[j] = 0.0
        for i in range(0, len(w)):
            for j in range(0, len(Z)): Z[j] += w[i]*expv[i*len(Z)+j]
        # normalize Z
        norm = sum(Z)
        for j in range(0, len(Z)): Z[j] /= norm
        # compute change in Z
        eps = 0.0
        for j in range(0, len(Z)):
            d = math.log(Z[j]/Z_old[j])
            eps += d*d
        # check convergence
        fconv.write("%d\t%.4e\n" %  (iii,eps))
        fconv.flush()
        
        if(eps<thres): break
    # return weights
    fconv.close()
    #print final z
    outz="z_"+str(iii)+".dat"
    fz = open(outz, 'w+')
    for i in range(0, len(Z)): fz.write("%32.30lf\n" % Z[i])
    fz.close()
    
    # compute final weights
    norm = 0.0
    for i in range(0, len(w)):
        ew = 0.0
        for j in range(0, len(Z)): ew += expv[i*len(Z)+j] / Z[j]
        w[i] = 1.0 / (ew+0.000000000001)
        norm += w[i]

    # normalize weights
    for i in range(0, len(w)):
        w[i] /= norm

    return w

# read FILENAME_

bias_ts=[]
for lines in open(FILENAME_, "r").readlines():
    riga=lines.strip().split()
    # skip comment lines
    if(riga[0]=="#!"): continue
    # read bias values
    # umbrella-sampling typical format
    if(len(riga) == NBIAS_+1):
       i0 = 1
       i1 = NBIAS_+1
    # bias exchange typical format
    elif(len(riga) == 2*NBIAS_+1):
       i0 = NBIAS_+1
       i1 = 2*NBIAS_+1
    # unknown format
    else:
       print(FILENAME_,"format is unknown!")
       exit()
    for i in range(i0, i1):
        bias_ts.append(float(riga[i]))

# number of frames
nframes = len(bias_ts) / NBIAS_

# printout
print("Number of frames::", nframes)
print("Number of bias entries::", len(bias_ts))

# get wham weights
if(nwham_i>=nwham):
    print("Asked for 0 wham iterations?")
    print("Relaunch changing number of cycles from starting one") 
else:
    if(nwham_i==0):
        ws = get_wham_weights(NBIAS_, int(nframes), bias_ts, deltaout, nwham_i, nwham)
    else:
        print("Restarting wham from iteration::",nwham_i)
        ws = get_wham_weights(NBIAS_, int(nframes), bias_ts, deltaout, nwham_i, nwham, filez)
     
    # printout WHAM weights to file
    print("Weights have been written to weights.dat")
    log = open("weights.dat", "w")
    for i in range(0, len(ws)): log.write("%32.30lf\n" % ws[i])
    log.close()

