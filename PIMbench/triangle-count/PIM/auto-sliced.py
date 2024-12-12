import subprocess
import re
import csv

# Define the base command and parameters
base_command = "./tc.out -c ../../../configs/iiswc/PIMeval_BitSerial_Rank4.cfg -v t -i \"../Dataset/v18772_symetric\" -a SLICED -s"
parameters = [4, 8, 16, 32, 64, 128, 256]

# Output CSV file
csv_file = "pim_results_sliced.csv"

# Regex patterns to extract data
data_copy_pattern = r"Data Copy Stats:.*?TOTAL --------- :\s+[\d,]+ bytes\s+([\d\.]+) ms Estimated Runtime\s+([\d\.]+) mj Estimated Energy"
pim_command_pattern = r"PIM Command Stats:.*?TOTAL --------- :\s+\d+\s+([\d\.]+)\s+([\d\.]+)"


# Initialize the CSV file
with open(csv_file, mode="w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["Num_of_Slices", "Data_Copy_Runtime_ms", "Data_Copy_Energy_mj", "PIM_Command_Runtime_ms", "PIM_Command_Energy_mj"])

# Execute commands and extract data
for param in parameters:
    command = f"{base_command} {param} > result_sliced.txt"
    print(f"Executing: {command}")

    # Run the command
    subprocess.run(command, shell=True)

    # Read the result_sliced file
    output_file = "result_sliced.txt"
    with open(output_file, "r") as file:
        content = file.read()

    # Extract data using regex
    data_copy_match = re.search(data_copy_pattern, content, re.DOTALL)
    pim_command_match = re.search(pim_command_pattern, content, re.DOTALL)

    if data_copy_match and pim_command_match:
        data_copy_runtime = float(data_copy_match.group(1))
        data_copy_energy = float(data_copy_match.group(2))
        pim_command_runtime = float(pim_command_match.group(1))
        pim_command_energy = float(pim_command_match.group(2))

        # Write data to CSV
        with open(csv_file, mode="a", newline="") as file:
            writer = csv.writer(file)
            writer.writerow([param, data_copy_runtime, data_copy_energy, pim_command_runtime, pim_command_energy])

        print(f"Data extracted for -s {param}")
    else:
        print(f"Error extracting data for -s {param}")

print(f"Results saved to {csv_file}")
