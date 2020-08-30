#-------------------------------------------------------------------------------
# Example script for HF_SHELL
# This script demonstrates the calculation of the pairing transition for 144Nd.
#
# At the end of execution, the code will have produced 
# a)  std.out  : all printing of the code to STDOUT.
# b)  Nd144.out: tabulated output of all calculations, suitable for 
#                plotting afterwards. 
#-------------------------------------------------------------------------------
if [ ! -d "work" ]; then
  mkdir work
fi

inter='Nd144'

cp ../../hf_shell.exe work/

cd work

#-------------------------------------------------------------------------------
# Creating the input data.
echo "Creating input data"
cat << EOF > in.data

&modelspace
spsfile   ="../inter/pn.sps", 
interfile ="../inter/${inter}.int", 
qfile     ="../inter/r2.red"
maxiter   =1000, printiter=100
outfile   ="Nd144.out"
ptype='HFB'
e_prec=1e-11
q_prec=1e-5
/
EOF

cat << EOF >> in.data
&config
neutrons   = 14,
protons    = 10,
inversetemp= -1
moreconfigs=.true.
/
EOF

for beta in  6.1 6.0 5.5 5.0 4.75 4.5 4.25 4.0 3.9 3.8 3.7 3.6 3.5 3.4 3.3 3.25 \
             3.2 3.15 3.1 3.0 2.9 2.8 2.7 2.6 2.5 2.4 2.3 2.2 2.1 2.05 2.0 1.975\
             1.95 1.9 1.75 1.50 1.25 
do

cat << EOF >> in.data
&config
neutrons   = 14,
protons    = 10,
inversetemp= $beta
moreconfigs=.true.
/
EOF
done
cat << EOF >> in.data
&config
neutrons   = 14,
protons    = 10,
inversetemp= 1.0
moreconfigs=.false.
/
EOF

#-------------------------------------------------------------------------------
# Run the code.
echo "Running the code."
time ./hf_shell.exe < in.data > std.out
#-------------------------------------------------------------------------------
echo "Storing output"
mv Nd144.out .. 
mv std.out   ..
#rm hf_shell.exe in.data
echo "Done"
#-------------------------------------------------------------------------------
