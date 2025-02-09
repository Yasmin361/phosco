import json
import os


def parse_fasta_blocks(file_path):
    """ Reads a FASTA-like file containing multiple three-line blocks and returns a list of components. """
    with open(file_path, "r") as file:
        lines = file.readlines()

    if len(lines) < 3:
        print(f"Error: File {file_path} does not contain enough lines.")
        return []

    parsed_entries = []
    for i in range(0, len(lines) - 2, 3):  # Process every three-line block
        fasta_header = lines[i].strip()  # Line 1: FASTA header
        sequence = lines[i + 1].strip()  # Line 2: Amino acid sequence
        topology = lines[i + 2].strip()  # Line 3: Topology prediction
        parsed_entries.append((fasta_header, sequence, topology))

    return parsed_entries


def extract_inner_membrane_sequences(sequence, topology):
    """ Extracts continuous amino acid stretches corresponding to 'I' regions in the topology. """
    if not sequence or not topology:
        return ["", "", "", ""], [0, 0, 0, 0]

    inner_membrane_regions = []
    lengths = []
    current_region = ""

    for i in range(len(topology)):
        if topology[i] == "I":
            current_region += sequence[i]  # Add the corresponding amino acid
        else:
            if current_region:  # If we have built an 'I' stretch, save it
                inner_membrane_regions.append(current_region)
                lengths.append(len(current_region))
                current_region = ""  # Reset for next stretch

    if current_region:  # Catch any trailing 'I' regions
        inner_membrane_regions.append(current_region)
        lengths.append(len(current_region))

    while len(inner_membrane_regions) < 4:
        inner_membrane_regions.append("")  # Fill missing regions with empty strings
        lengths.append(0)

    return inner_membrane_regions[:4], lengths[:4]


def append_to_json(output_json, new_entry):
    """ Appends a new entry to an existing JSON file or creates a new one if it doesn't exist. """
    if os.path.exists(output_json):
        with open(output_json, "r") as json_file:
            try:
                data = json.load(json_file)
            except json.JSONDecodeError:
                print(f"Warning: {output_json} was not properly formatted. Resetting file.")
                data = []
    else:
        data = []

    data.append(new_entry)

    with open(output_json, "w") as json_file:
        json.dump(data, json_file, indent=4)


def process_fasta_files(input_folder, output_json):
    """ Process all FASTA-like files in the input folder and save results to a JSON file. """
    for input_file in os.listdir(input_folder):
        file_path = os.path.join(input_folder, input_file)

        if os.path.isfile(file_path):
            print(f"Processing: {input_file}")
            parsed_entries = parse_fasta_blocks(file_path)

            for fasta_header, sequence, topology in parsed_entries:
                inner_membrane_regions, lengths = extract_inner_membrane_sequences(sequence, topology)
                total_sequence_length = len(sequence)
                total_inner_membrane_length = sum(lengths)
                ratio_inner_to_total = total_inner_membrane_length / total_sequence_length if total_sequence_length > 0 else 0

                new_entry = {
                    "file": input_file,
                    "header": fasta_header,
                    "ICL1": inner_membrane_regions[0],
                    "ICL2": inner_membrane_regions[1],
                    "ICL3": inner_membrane_regions[2],
                    "CTER": inner_membrane_regions[3],
                    "len_ICL1": lengths[0],
                    "len_ICL2": lengths[1],
                    "len_ICL3": lengths[2],
                    "len_CTER": lengths[3],
                    "total_sequence_length": total_sequence_length,
                    "ratio_inner_to_total": ratio_inner_to_total
                }

                append_to_json(output_json, new_entry)

    print(f"Processing complete. Data saved in '{output_json}'.")


input_folder = r"C:\Users\Yasmin\PycharmProjects\phosco\preprocessed_full\predictions"
output_json = "loop_lengths_full.json"
process_fasta_files(input_folder, output_json)
