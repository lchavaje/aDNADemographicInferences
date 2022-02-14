######################################
# Simulation configuration variables #
######################################

# CHANGE TO WHEREVER YOU'RE RUNNING THIS
g_dir="/mnt/Timina/dortega/lchavaje/AncientDNAProject/AncientDNAProject/Scripts/JazepsPipeline/Paleogenomic-Datasim-master/data_simulation"
size_freq_path="/cm/shared/apps/gargammel/6oct2017/src/sizefreq.size.gz"

# Parameter file
parameter_file=$1

#Number of ancient and modern samples
ancients=5
moderns=30
# Reads will be simulated for all combinations of coverages and contaminations,
# coverages=( "31" "30" "5" "1" "0.1" "0.01" ) indicate 30X, 5X, 1X, 0.1 and 0.01X coverages, 31 indicates 30x without damage 
# contaminations=( "0" ) indicate 0% modern contamination
coverages=( "31" "30" "5" "1" "0.1" "0.01" )
contaminations=( "0" )

# Set to false if you want to skip the SHAPEIT phasing step
#should_phase=true

###############################
# End configuration variables #
###############################


# Transform coverage and contamination into comma separated strings,
# to simplify passing them as arguments into python scripts
printf -v cov_string "%s," "${coverages[@]}"
cov_string=${cov_string%?}
printf -v con_string "%s," "${contaminations[@]}"
con_string=${con_string%?}

# Fresh directory structure
rm -rf "$g_dir"/logs
mkdir "$g_dir"/logs
mkdir "$g_dir"/logs/err/
mkdir "$g_dir"/logs/std/
rm -rf "$g_dir"/present/
rm -rf "$g_dir"/reference/
rm -rf "$g_dir"/cases/

echo "OK before pre_gargammel"

# pre_gargammel.sh: simulated sequences and directory structure required for gargammel, reference panel creation
qsub -o "$g_dir"/logs/std/pre_gargammel -e "$g_dir"/logs/err/pre_gargammel \
	"$g_dir"/scripts/pre_gargammel_NewPipeline.sh $g_dir/$parameter_file $g_dir \
	$cov_string $con_string

echo "OK after pre_gargammel"

echo "OK before gargammel_anc"

# gargammel_anc.sh: ancient samples simulated reads for all combinations of coverage and contamination
qsub -t 1-$ancients -o "$g_dir"/logs/std/gargammel_anc -e "$g_dir"/logs/err/gargammel_anc \
	"$g_dir"/scripts/gargammel_anc.sh $g_dir $size_freq_path "${#coverages[@]}" "${coverages[@]}" \
	"${#contaminations[@]}" "${contaminations[@]}"

echo "OK after gargammel_anc"

echo "OK before gargammel_pre"

# gargammel_pre.sh: present samples simulated reads for all combinations of coverage and contamination
qsub -t 1-$moderns -o "$g_dir"/logs/std/gargammel_pre -e "$g_dir"/logs/err/gargammel_pre \
        "$g_dir"/scripts/gargammel_pre.sh $g_dir 

echo "OK after gargammel_pre"

echo "OK before pre_map"

# pre_map.sh: indexing of reference sequence
qsub -o "$g_dir"/logs/std/pre_map -e "$g_dir"/logs/err/pre_map \
	"$g_dir"/scripts/pre_map.sh $g_dir

echo "Ok after pre_map"

echo "Ok before map_anc"

# map_anc.sh: mapping of all ancient simulated samples to reference sequence
qsub -t 1-$ancients -o "$g_dir"/logs/std/map_anc -e "$g_dir"/logs/err/map_anc \
	"$g_dir"/scripts/map_anc.sh $g_dir "${#coverages[@]}" "${coverages[@]}" \
	"${#contaminations[@]}" "${contaminations[@]}"

echo "Ok after map_anc"

echo "Ok before map_pre"

# map_pre.sh: mapping of all present simulated samples to reference sequence
qsub -t 1-$moderns -o "$g_dir"/logs/std/map_pre -e "$g_dir"/logs/err/map_pre \
       "$g_dir"/scripts/map_pre.sh $g_dir 

echo "Ok after map_pre"

# call_variants.sh: creation and filtering of VCF files for all simulated samples
# qsub -t 1-$ancients -o "$g_dir"/logs/std/call_variants -e "$g_dir"/logs/err/call_variants \
#	"$g_dir"/scripts/call_variants.sh $g_dir "${#coverages[@]}" "${coverages[@]}" \
#	"${#contaminations[@]}" "${contaminations[@]}"

# if [ "$should_phase" = true ]
# then

	# phase.sh: phasing of all samples using SHAPEIT, accuracy and switch error reports
#	qsub -t 1-$ancients -o "$g_dir"/logs/std/phase -e "$g_dir"/logs/err/phase \
#		"$g_dir"/scripts/phase.sh $g_dir "${#coverages[@]}" "${coverages[@]}" \
#		"${#contaminations[@]}" "${contaminations[@]}"

	# results.sh: compilation of phasing results for all samples
#	qsub -o "$g_dir"/logs/std/results -e "$g_dir"/logs/err/results \
#		"$g_dir"/scripts/results.sh $g_dir $anc_age $ancients $references $bases
# fi
