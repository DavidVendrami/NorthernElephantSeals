SLiM input files are in the 'SLiM_files' folder. 
Files:
- WF_delet_along_whole_genome_RAD_rep1: based on RAD Ne estimates with deleterious mutations arising anywhere
- WF_delet_exons_RAD_rep1: based on RAD Ne estimates with deleterious mutations arising only in exons
- WF_WGS_rep1: based on WGS Ne estimates with deleterious mutations arising anywhere

To replicate the SLiM workflow one need to create 100 'rep_' folders (rep_1, rep_2, ..., rep_100) and
create 100 copies of the relevant input file (make sure to change the output location in each file! file 1 will
output in rep_1, file 2 will output in rep_2 and so on. This can be quickly done in bash as follows:
for i in {2..100}
do
sed "s/\/rep_1\//\/rep_${i}\//g" INPUT_rep1.txt > INPUT_rep$i
done
)

Then SLiM can be run using the 'run_SLiM.sh' file.

SLiM outputs then need to be processed following the 'Process_output_v2.sh' script. Make sure to use the right Ne 
table (see details within the script).

After that stuff can be plotted using the various plotting r scripts.

The 'Purge.r' take the outputs from 'Process_output_v2.sh', generates number of lost deleterious mutations and 
then plot them.

--------------------------------------------------------------------------------------------------

To generate equivalent plots using the Ne estimates obtained from the 100 bootstrapped SFS, use the material in 
the 'Boot_estimates' folder.

First, generate SLiM input files and run them straight away with the 'GenerateInputs_and_RunSlim.sh' file. This will
take the 'WF_delet_along_whole_genome_F_repN.txt' file and use it as a template to create the actual files with the various
bootstrap estimates.

Then, first generate the Ne Table using the 'Generate_Ne_Table.r' script and then process the slim output with 
the 'Process_output_Boot_NeTable' script. 

Plots with plot scripts.