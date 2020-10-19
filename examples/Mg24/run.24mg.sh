################################################################################
# Example calculation of the triaxial (Q,gamma) plane for 24Mg using the 
# USDB interaction.
#
# We have organized the calculation according to deformation: for every 
# deformation we start from beta=0.4 and descend to beta=10.0. Hence, a file is
# produced for every deformation, and the user will need to do some 
# postprocessing to obtain a final figure.
#
# This is a somewhat difficult calculation because a) the system is so light and
# b) pairing correlations are almost always very weak. Obtaining the minimum is 
# straightforward with standard parameters of the code, but obtaining the 
# full surface is not.
#
# At the end of execution, the code will have produced a set of data files 
# 
#   tables/Mg.$q1.$q2.out
# 
# with all tabulated output for that particular combination of q1,q2, for 
# four different values of the inverse temperature beta.
################################################################################

inter='usdb.24'

if [ ! -d "work" ]; then
  mkdir work
fi

if [ ! -d "tables" ]; then
  mkdir tables
fi

cp ../../hf_shell.exe work/

cd work

Z=4
N=4


for q1 in  `seq 0 2 19`
do

maxq2=12

if(( q1 > 5 )); then
maxq2=11
fi

if(( q1 > 8 )); then
maxq2=7
fi

if(( q1 > 12 )); then
maxq2=7
fi

if(( q1 > 16 )); then
maxq2=6
fi

if(( q1 > 17 )); then
maxq2=4
fi

if(( q1 > 19 )); then
maxq2=1
fi


################################################################################
# First all calculations with q2 = 0
q2=0
cat << EOF > in.data
&modelspace
spsfile   ="../$inter/pn.sps", 
interfile ="../$inter/$inter.int", 
qfile     =""
maxiter   =4000, printiter=100
ptype     ='BCS'
outfile="Mg24.$q1.$q2.out"
e_prec=1e-12
q_prec=1e-6
/
EOF

for beta in 0.4 0.6 0.8 10.0 
do

# note that for this example, we manually set stepsize and momentum parameters
# as well as a slowdown for the mixing of the density. 
# since there is no pairing, the surface is discontinous at low temperatures.
# this means the lagrange parameters of the quadrupole constraints can grow very
# large during the iterative process, and convergence is very difficult.
#
# this is not necessary for high temperatures, but we have left this in to 
# make the scripting easier.
cat << EOF >> in.data
&config
neutrons   = $N,
protons    = $Z,
inversetemp= $beta
constrainttype = 2
q1target =$q1
q2target =$q2
moreconfigs=.true.
stepsize=0.0075
momentum=0.0
denmix=0.1
lambda20=0
lambda22=0
/
EOF

done
echo "Running q1=$q1, q2=$q2"
./hf_shell.exe < in.data > HF.out
mv Mg24.$q1.$q2.out ../tables

################################################################################
# Then all calculations with q2 != 0
for q2 in `seq 1 2 $maxq2`
do
cat << EOF > in.data
&modelspace
spsfile   ="../$inter/pn.sps", 
interfile ="../$inter/$inter.int", 
qfile     =""
maxiter   =2000, printiter=100
ptype     ='BCS'
outfile="Mg24.$q1.$q2.out"
e_prec=1e-9
q_prec=1e-3
/
EOF

for beta in 0.4 0.6 0.8 10.0 
do
cat << EOF >> in.data
&config
neutrons   = $N,
protons    = $Z,
inversetemp= $beta
constrainttype = 2
q1target =$q1
q2target =$q2
moreconfigs=.true.
stepsize=0.0075
momentum=0.0
denmix=0.05
lambda20=0
lambda22=0
/
EOF

done
echo "Running q1=$q1, q2=$q2"
./hf_shell.exe < in.data > HF.out
mv Mg24.$q1.$q2.out ../tables/
done
done

cat << EOF > in.data
&modelspace
spsfile   ="../$inter/pn.sps", 
interfile ="../$inter/$inter.int", 
qfile     =""
maxiter   =2000, printiter=100
ptype     ='HFB'
outfile="Mg24.uncon.out"
e_prec=1e-9
q_prec=1e-3
/
EOF

cat << EOF >> in.data
&config
neutrons   = $N,
protons    = $Z,
inversetemp= 0.4
moreconfigs=.true.
/
&config
neutrons   = $N,
protons    = $Z,
inversetemp= 0.6
constraintiter=10
constrainttype=2
q1target =4.0
q2target =0.0
moreconfigs=.true.
/
&config
neutrons   = $N,
protons    = $Z,
inversetemp= 0.8
constraintiter=10
q1target =10.0
q2target =4.0
moreconfigs=.true.
/
&config
neutrons   = $N,
protons    = $Z,
inversetemp= 10.0
constrainttype = 2
constraintiter=10
q1target =10.0
q2target =4.0
moreconfigs=.true.
/
EOF

./hf_shell.exe < in.data > HF.out
mv Mg24.uncon.out ../tables/
rm hf_shell.exe
