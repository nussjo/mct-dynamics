#!/usr/bin/python
import sys

# check command line
if (len(sys.argv)!=2):
	print "Usage: "+sys.argv[0]+" <filename>"
	sys.exit(0)

inputfile = open(sys.argv[1], "r")

# columns for entries; start with 0
msdentry = 1

# read file, save all Ns and xs into arrays
t=[]
msd=[]
lines = inputfile.readlines()
t.append(0.0)
msd.append(0.0)
for line in lines:
	# ignore comment lines:
	if (line[0]=='#'):
		continue
	# ignore empty lines
	if (len(line)==0):
		continue
	columns = line.split()
	# ignore undetected columns
	if (len(columns)==0):
		continue
	# extract values
	tcur = float(columns[0])
	msdcur = float(columns[msdentry])
        # append to array
	t.append(tcur)
	msd.append(msdcur)
inputfile.close()

# then central diff for the remaining values
print "#1:t 2:diff_msd " 
for i in range(1, len(t)-1):
    dt1 = t[i]-t[i-1];
    dt2 = t[i+1]-t[i];
    if dt1 <= 0: continue
    if dt2 <= 0: continue
    diff_msd = ((msd[i+1] - msd[i])/dt2 + (msd[i]- msd[i-1])/dt1)/12.0 # 2*6.0
    print (t[i+1]+t[i-1])/2.0, diff_msd
