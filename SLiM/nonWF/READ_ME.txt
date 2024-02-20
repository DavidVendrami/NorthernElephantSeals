SLiM input files are in the 'SLiM_inputs' folder. 
There are two 'Master' files (one for the complete model and one for the neutral model) that can be used to generate the various
replicates.

To replicate the SLiM workflow one need to create 100 'rep_' folders (rep_1, rep_2, ..., rep_100) and
create 100 copies of the relevant input file (Make sure to set the bottleneck carrying capacity to the right value!
make also sure to change the output location in each file! file 1 will
output in rep_1, file 2 will output in rep_2 and so on. This can be quickly done in bash as follows:
for i in {1..100}
do
sed 's/BOTT_SIZE/DESIRED_SIZE/g' Wout_complete_Master.txt > Wout_complete_MasterP2.txt # DESIRED_SIZE should eb 50, 100, 250, 500 or 1000
sed "s/PATH_TO_repN_directory/YOUR_PATH_TO_\/rep_${i}/g" Wout_complete_MasterP2.txt > Wout_complete_rep$i.txt
done
)

Then SLiM can be run as for the WF models, but supplying enough ram (in my case 24G).

SLiM outputs then need to be processed following the 'Parse_nonWF_andMakeSpace.bash' script. 

After that stuff can be plotted using the various plotting r scripts.

