######################################
# Simulation configuration variables #
######################################

# CHANGE TO WHEREVER YOU'RE RUNNING THIS
g_dir="/mnt/Timina/dortega/lchavaje/AncientDNAProject/AncientDNAProject/Scripts/JazepsPipeline/Paleogenomic-Datasim-master/data_simulation"
size_freq_path="../../../../Results/Paleogenomic-Datasim/size_freq"

# Age in generations of samples with simulated reads
anc_age=$1
# Amount of simulated ancient samples
ancients=$2
# Amount of present day sequences used to create reference panel
references=$3
# Length of sequences to simulate
bases=$4
# Branch scaling factor when passing from msprime to seq-gen
branch_scale=$5
# OPTIONAL Demographic event to simulate, currently supported: [split, bottleneck]
event=$6
# OPTIONAL Time in generations from present when event happened
event_time=$7

# Reads will be simulated for all combinations of coverages and contaminations,
# coverages=( "10" "5" "1" ) indicate 10X, 5X, and 1X coverage,
# contaminations=( "0" "2" "5" ) indicate 0%, 2%, and 5% modern contamination
coverages=( "10" "5" "1" )
contaminations=( "0" "2" "5" "10" )

# Set to false if you want to skip the SHAPEIT phasing step
should_phase=true

###############################
# End configuration variables #
###############################

# event and event_time should be mutually inclusive
if [ -z "$6" ] || [ -z "$7" ]
then
	event="null"
	event_time="null"
fi

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

# pre_gargammel.sh: simulated sequences and directory structure required for gargammel, reference panel creation
qsub -o "$g_dir"/logs/std/pre_gargammel -e "$g_dir"/logs/err/pre_gargammel \
	"$g_dir"/scripts/pre_gargammel.sh $g_dir $anc_age $ancients $references $bases $branch_scale \
	$event $event_time $cov_string $con_string

# gargammel.sh: simulated reads for all combinations of coverage and contamination
qsub -t 1-$ancients -o "$g_dir"/logs/std/gargammel -e "$g_dir"/logs/err/gargammel \
	"$g_dir"/scripts/gargammel.sh $g_dir $size_freq_path "${#coverages[@]}" "${coverages[@]}" \
	"${#contaminations[@]}" "${contaminations[@]}"

# pre_map.sh: indexing of reference sequence
qsub -o "$g_dir"/logs/std/pre_map -e "$g_dir"/logs/err/pre_map \
	"$g_dir"/scripts/pre_map.sh $g_dir

# map.sh: mapping of all simulated samples to reference sequence
qsub -t 1-$ancients -o "$g_dir"/logs/std/map -e "$g_dir"/logs/err/map \
	"$g_dir"/scripts/map.sh $g_dir "${#coverages[@]}" "${coverages[@]}" \
	"${#contaminations[@]}" "${contaminations[@]}"

# call_variants.sh: creation and filtering of VCF files for all simulated samples
qsub -t 1-$ancients -o "$g_dir"/logs/std/call_variants -e "$g_dir"/logs/err/call_variants \
	"$g_dir"/scripts/call_variants.sh $g_dir "${#coverages[@]}" "${coverages[@]}" \
	"${#contaminations[@]}" "${contaminations[@]}"

if [ "$should_phase" = true ]
then

	# phase.sh: phasing of all samples using SHAPEIT, accuracy and switch error reports
	qsub -t 1-$ancients -o "$g_dir"/logs/std/phase -e "$g_dir"/logs/err/phase \
		"$g_dir"/scripts/phase.sh $g_dir "${#coverages[@]}" "${coverages[@]}" \
		"${#contaminations[@]}" "${contaminations[@]}"

	# results.sh: compilation of phasing results for all samples
	qsub -o "$g_dir"/logs/std/results -e "$g_dir"/logs/err/results \
		"$g_dir"/scripts/results.sh $g_dir $anc_age $ancients $references $bases
fi
