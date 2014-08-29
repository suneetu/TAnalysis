# Script to merge HFTs
inputDIR="/gdata/atlas/suneetu/Documents/LFV_Higgs2014/output/R1_August_16/"

numero=(
"110805" 
"110806" 
"110807" 
"110808" 
"110817" 
"110818" 
"110819" 
"110820" 
"110809" 
"110810" 
"110811" 
"110812" 
"110821" 
"110822" 
"110823" 
"110824" 
"110813" 
"110814" 
"110815" 
"110816" 
"110825" 
"110826" 
"110827" 
"110828" 
)


for j in "${!numero[@]}"
do
  mv -v ${inputDIR}'CENTRAL_'${numero[$j]}.root ${inputDIR}'NOM_'${numero[$j]}.root
done
