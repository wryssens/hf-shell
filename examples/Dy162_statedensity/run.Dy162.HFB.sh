#-------------------------------------------------------------------------------
# Example script for HF_SHELL
# This script demonstrates the calculation of the state density for 162Dy.
# We perform calculations for 493 values of the inverse temperature, which are
# contained in the file betas.dat.
#
# Note that this is a somewhat longer calculation than the other examples.
# On the modest laptop of one of the authors, this script executes in roughly 
# ten minutes.
#
# At the end of execution, the code will have produced 
# a)  std.out  : all printing of the code to STDOUT.
# b)  Dy162.out: tabulated output of all calculations, suitable for 
#                plotting afterwards. 
#-------------------------------------------------------------------------------
if [ ! -d "work" ]; then
  mkdir work
fi

inter='Dy162'

cp ../../hf_shell.exe work/
cp betas.dat          work/

cd work

#-------------------------------------------------------------------------------
# Creating the input data.
echo "Creating input data"
cat << EOF > in.data

&modelspace
spsfile   ="../inter/Dy162.sps", 
interfile ="../inter/${inter}.int", 
qfile     ="../inter/r2.red"
maxiter   =1000, printiter=100
outfile   ="Dy162.out"
ptype='HFB'
e_prec=1e-9
q_prec=1e-5
/
EOF
# Note how we set constraintiter to 10, to use a temporary constraint such that
# the nucleus is able to break rotational symmetry at the start of the iterative 
# process and we can find the deformed configuration.
cat << EOF >> in.data
&config
neutrons   = 26,
protons    = 16,
inversetemp= 30
moreconfigs=.true.
Q20target=100
constraintiter=10
constrainttype=2
moreconfigs=.true.
/
EOF
# We add a &config namelist for every inverse temperature considered.
while read p; do
cat << EOF >> in.data
&config
inversetemp= $p
moreconfigs=.true.
/
EOF
done < betas.dat
# And a final one
cat << EOF >> in.data
&config
inversetemp= 0.003906
moreconfigs=.false.
/
EOF

#-------------------------------------------------------------------------------
# Run the code.
echo "Running the code."
time ./hf_shell.exe < in.data |tee std.out
#-------------------------------------------------------------------------------
echo "Storing output"
mv Dy162.out .. 
mv std.out   ..
rm hf_shell.exe in.data
echo "Done"
#-------------------------------------------------------------------------------
