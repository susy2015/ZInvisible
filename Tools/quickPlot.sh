# quickPlots.sh
# Caleb J. Smith
# October 11, 2018

# Intended for quick use and multiple plots (not on condor)

# Outline

# Compile MakePlots
# Remove existing plots
# Save root files for different MC
# hadd the results
# Make plots of results (containing both MC)


combineResults=true

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
declare -a samples=(
                    "Data_SingleMuon_2016"
                    "DYJetsToLL_HT_400to600"
                    "IncDY"
                    "TTbarNoHad"
                    "SingleTopZinv"
                    "TTZ"
                    "Diboson"
                    "Rare"
                    "Data_SinglePhoton_2016"
                    "GJets_HT-400To600"
                    "QCD_HT500to700"
                    "TTGJets"
                    "WJetsToLNu_HT_400to600"
                    "TTbarAll"
                    "tW"
                    "ZJetsToNuNu_HT_400to600"
                   )

outputFiles=
n_events=1000

# loop through samples array
for sample in "${samples[@]}"
do
    output=""$sample".root"
    outputFiles="$outputFiles $output"
    echo " - Running makePlots to create $output"
    echo "./makePlots -D $sample -E $n_events -I $output | grep -v LHAPDF"
    ./makePlots -D $sample -E $n_events -I $output | grep -v LHAPDF
done

# --- for only two samples --- #

#sample1=GJets_HT-400To600
#sample1=DYJetsToLL_HT_400to600
#sample2=ZJetsToNuNu_HT_400to600
#output1=""$sample1".root"
#output2=""$sample2".root"

#echo " - Running makePlots to create $output1"
#echo "./makePlots -D $sample1 -E $n_events -I $output1 | grep -v LHAPDF"
#./makePlots -D $sample1 -E $n_events -I $output1 | grep -v LHAPDF
#echo " - Running makePlots to create $output2"
#echo "./makePlots -D $sample2 -E $n_events -I $output2 | grep -v LHAPDF"
#./makePlots -D $sample2 -E $n_events -I $output2 | grep -v LHAPDF

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
    ./makePlots -f -I $result | grep -v LHAPDF
fi
