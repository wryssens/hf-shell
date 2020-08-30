#-------------------------------------------------------------------------------
# Example script for HF_SHELL
# This script demonstrates the calculation of the shape transition for 162Dy.
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
inversetemp= -1
moreconfigs=.true.
Q20target=100
constraintiter=10
constrainttype=2
moreconfigs=.true.
/
EOF
# We add a &config namelist for every inverse temperature considered.
for p in 0.90 0.89 0.87 0.86 0.85 0.84 0.839 0.838 0.837 0.836 0.8365 0.8360 0.8355 \
         0.835 0.8349 0.8348 0.8347 0.8346 0.8345 0.834 0.833 0.832 0.831 0.83 \
         0.82 0.81
do 

cat << EOF >> in.data
&config
neutrons   = 26,
protons    = 16,
inversetemp= $p
moreconfigs=.true.
/
EOF
done

cat << EOF >> in.data
&config
neutrons   = 26,
protons    = 16,
inversetemp= 0.80
moreconfigs=.false.
/
EOF
#-------------------------------------------------------------------------------
# Run the code.
echo "Running the code."
./hf_shell.exe < in.data > std.out
#-------------------------------------------------------------------------------
echo "Storing output"
mv Dy162.out .. 
mv std.out   ..
rm hf_shell.exe in.data
echo "Done"
#-------------------------------------------------------------------------------
