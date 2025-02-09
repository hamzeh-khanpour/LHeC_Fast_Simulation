import os

# ✅ List of LHE files to merge
lhe_files = [
    "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_tautau_SM_NP_2_SMEFTsim_top_alphaScheme_UFO/Events/run_06/aa_tautau_SM_NP_2_SMEFTsim_top_alphaScheme_UFO_1.lhe",
    "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_tautau_SM_NP_2_SMEFTsim_top_alphaScheme_UFO/Events/run_07/aa_tautau_SM_NP_2_SMEFTsim_top_alphaScheme_UFO_2.lhe",
    "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_tautau_SM_NP_2_SMEFTsim_top_alphaScheme_UFO/Events/run_08/aa_tautau_SM_NP_2_SMEFTsim_top_alphaScheme_UFO_3.lhe",
    "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_tautau_SM_NP_2_SMEFTsim_top_alphaScheme_UFO/Events/run_09/aa_tautau_SM_NP_2_SMEFTsim_top_alphaScheme_UFO_4.lhe",
    "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_tautau_SM_NP_2_SMEFTsim_top_alphaScheme_UFO/Events/run_10/aa_tautau_SM_NP_2_SMEFTsim_top_alphaScheme_UFO_5.lhe",
]

# ✅ Output file name
output_file = "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_tautau_SM_NP_2_SMEFTsim_top_alphaScheme_UFO/merged_aa_tautau_SM_NP_2_SMEFTsim_top_alphaScheme_UFO.lhe"

# ✅ Function to merge LHE files
def merge_lhe_files(lhe_files, output_file):
    with open(output_file, "w") as outfile:
        first_file = True
        for lhe_file in lhe_files:
            with open(lhe_file, "r") as infile:
                for line in infile:
                    # ✅ Write header only once (from the first file)
                    if "<LesHouchesEvents" in line and not first_file:
                        continue  # Skip duplicate headers
                    if "</LesHouchesEvents>" in line:
                        continue  # Skip duplicate footers
                    outfile.write(line)
            first_file = False
        
        # ✅ Add closing tag at the end
        outfile.write("</LesHouchesEvents>\n")

# ✅ Run the merging function
merge_lhe_files(lhe_files, output_file)

print(f"✅ Merged LHE file saved as: {output_file}")

