################################################################################
# Example calculation to find the triaxial mean-field minimum of 24Mg using the
# USD interaction.  At the end of execution, this script diplays the total
# energy and quadrupole deformation calculated by HF-SHELL, to be compared to 
# the values published in Gao et al. PRC 92, 064310 (2015). 
#
# The script produces as well an output file ""minimum.out".
#
################################################################################
inter='usdb'

if [ ! -d "work" ]; then
  mkdir work
fi

if [ ! -d "tables" ]; then
  mkdir tables
fi

cp ../../hf_shell.exe work/
#cp $HOME/Documents/Codes/hf_shell/hf_shell.exe work/

cd work

# Note that we scale the TBMEs with the prescribed USDB scaling.
cat << EOF > in.data
&modelspace
spsfile   ="../$inter/pn.sps", 
interfile ="../$inter/$inter.int", 
qfile     =""
maxiter   =2000, printiter=100
ptype     ='HFB'
outfile="minimum.out"
e_prec=1e-12
q_prec=1e-6
TBME_A    = 24
TBME_Aref = 18
TBME_x = 0.3
/
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# To break rotational and axial symmetry at the start of the iterations, we 
# impose a small constraint on both Q20 and Q22. Note that we use the 
# constraintiter keyword to use this constraint only during the first two
# iterations, after that the calculation proceeds unconstrained.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
&config
neutrons   = 4,
protons    = 4,
inversetemp= -1.0
constrainttype= 2
q1target   =1
q2target   =2
constraintiter=2
moreconfigs=.false.
/
EOF

./hf_shell.exe < in.data > std.out

echo "------------------------------------------------------"
echo "HF-SHELL result"
grep "Total energy" std.out| tail -1
echo "                      proton      neutron      total"
grep "Q     " std.out| tail -1
grep "Gamma " std.out| tail -1
echo "------------------------------------------------------"
echo "Result of Gao et al. PRC 92, 064310 (2015)"
echo " Total Energy     (MeV)  =  -80.965"
echo " (Q, gamma)              =  (18.7, 12.0)"
echo "------------------------------------------------------"
#rm std.out     
mv minimum.out ../
