# This script assume that outputs for different simulations are stored in 
# different folder called "rep_N". Where N is the number ID of each simulation.

### 1. Extract population sizes ###

for i in rep_*
do
grep 'Pop_size at' $i/F_*.txt | cut -d ':' -f2 > $i/pop_size.txt
done

paste rep_*/pop_size.txt > Pop_sizes.txt

###
## Get genetic loads

# parse files:
for i in rep*
do
cd $i
for j in gen_*.txt
do
sed -n '/Mutations:/,/Individuals:/p' $j | grep -v 'Mutations:' | grep -v 'Individuals' > ${j%.txt}_parsed.txt
done
rm gen_????.txt
cd ../
done

# Genetic load, realized genetic load, masked genetic load, absolute numbers

for i in rep*
do
cd $i
ls -1 *_parsed.txt > file_list.txt
cd ../
done

for i in rep_*
do
cd $i
Rscript /PATH/TO/Parse_simulations_output.R pop_size.txt
Rscript /PATH/TO/Drift_Load.R pop_size.txt
cd ../
done

####################################################

# Merge outputs from different simulations
paste rep_*/GL.txt > GL_all.txt
paste rep_*/RGL.txt > RGL_all.txt
paste rep_*/MGL.txt > MGL_all.txt
paste rep_*/DL.txt > DL_all.txt
