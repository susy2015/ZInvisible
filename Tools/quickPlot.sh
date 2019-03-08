# quickPlots.sh
# Caleb J. Smith
# October 11, 2018

# Intended for quick use and multiple plots (not on condor)

# Overview
# - Compile MakePlots
# - Remove existing plots
# - Make plots for different samples (saved to root files)
# - hadd the results
# - Make plots of results

combineResults=true
year="2017"

# Compile MakePlots

echo " - Compile MakePlots"
make -j8

retVal=$?
if [ $retVal -ne 0 ]; then
    echo "ERROR: There was an error when compiling MakePlots"
    exit $retVal;
fi

# Remove existing plots

echo " - Remove existing plots"
rm plots/*

# Run on all samples
# declare an array of samples
# Z to LL and photon samples for Data and MC
# don't run a sample twice 
# the samples TTZ, Diboson, Rare are used for both lepton and photon
# declare -a samples=(
#                     "Data_SingleMuon_2016"
#                     "DYJetsToLL_HT_400to600"
#                     "TTbarNoHad"
#                     "SingleTopZinv"
#                     "Rare"
#                     "TTZ"
#                     "Diboson"
#                     "Data_SinglePhoton_2016"
#                     "GJets_HT-400To600"
#                     "QCD_HT500to700"
#                     "WJetsToLNu_HT_400to600"
#                     "TTbarAll"
#                     "tW"
#                     "ZJetsToNuNu_HT_400to600"
#                    )
#declare -a samples=(
#                    "Data_SingleMuon_2016"
#                    "DYJetsToLL_HT_400to600"
#                   )
#declare -a samples=(
#                    "Data_SingleMuon_"$year""
#                    "DYJetsToLL_HT_400to600_"$year""
#                    "Data_SinglePhoton_"$year""
#                    "GJets_HT-400To600_"$year""
#                    "ZJetsToNuNu_HT_400to600_"$year""
#                   )
declare -a samples=(
                    "DYJetsToLL_HT_400to600_"$year""
                    "GJets_HT-400To600_"$year""
                    "ZJetsToNuNu_HT_400to600_"$year""
                   )

outputFiles=
n_events=2000

# loop through samples array
for sample in "${samples[@]}"
do
    output=""$sample".root"
    outputFiles="$outputFiles $output"
    echo " - Running makePlots to create $output"
    echo "./makePlots -D $sample -E $n_events -I $output -Y $year | grep -v LHAPDF"
    ./makePlots -D $sample -E $n_events -I $output -Y $year | grep -v LHAPDF
done

# hadd the results
if [ "$combineResults" = true ]; then

    # IMPORTANT: You need quotes around each variable (due to underscores).
    #            Otherwise bash interprets the underscores as part of the variable names.
    #result=""$sample1"_and_"$sample2".root"
    result="quickResult.root"
    echo " - hadd the results to create $result"
    # Use -f to overwrite target file if it already exists
    #hadd -f $result $output1 $output2
    hadd -f $result $outputFiles
    
    # Make plots of results (containing both MC)
    
    echo " - Make plots of results (containing both MC)"
    echo "./makePlots -f -I $result | grep -v LHAPDF"
    ./makePlots -f -I $result -Y $year | grep -v LHAPDF
fi
