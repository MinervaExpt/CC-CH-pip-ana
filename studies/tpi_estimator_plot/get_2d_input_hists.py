from ROOT import TFile, gDirectory, PlotUtils

# Open the input ROOT file
input_file = TFile.Open("scaled_MC.root", "READ")

# Create the output ROOT file
output_file = TFile("out.root", "RECREATE")

# Loop over all keys in the input file
for key in input_file.GetListOfKeys():
    # Get the object (assumed to be a histogram)
    obj = key.ReadObj()

    # Check if the object is a TH2D and its name contains "Pion_Range_vs_Pion"
    if "Pion_Range_vs_Pion" in obj.GetName():
        # Write the histogram to the output file
        obj.Write()

    if "POTUsed" in obj.GetName():
        obj.Write()

# Close the output file
output_file.Close()

# Close the input file
input_file.Close()
