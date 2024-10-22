def create_subset_files(input_file, lines_125M, lines_1M, output_file_125M, output_file_1M):
    """
    Creates two subset files from the input file: one with 125 million lines and another with 1 million lines.
    
    :param input_file: The path to the large input file
    :param lines_125M: Number of lines for the 125M file (125,000,000)
    :param lines_1M: Number of lines for the 1M file (1,000,000)
    :param output_file_125M: The path for the output file containing 125M lines
    :param output_file_1M: The path for the output file containing 1M lines
    """
    with open(input_file, 'r') as infile:
        # Open the two output files
        with open(output_file_125M, 'w') as outfile_125M, open(output_file_1M, 'w') as outfile_1M:
            for line_number, line in enumerate(infile):
                # Write to the 125M file until we hit the desired number of lines
                if line_number < lines_125M:
                    outfile_125M.write(line)
                
                # Write to the 1M file until we hit the desired number of lines
                if line_number < lines_1M:
                    outfile_1M.write(line)
                
                # If we have written enough lines to both files, break out of the loop
                if line_number >= max(lines_125M, lines_1M):
                    break

if __name__ == "__main__":
    # File paths
    input_file = "/home/afloresep/work/chelombus/data/650M.cxsmiles"  # Input file
    output_file_125M = "/home/afloresep/work/chelombus/data/125M.cxsmiles"  # Output file for 125M lines
    output_file_1M = "/home/afloresep/work/chelombus/data/1M.cxsmiles"  # Output file for 1M lines
    
    # Number of lines to write to each file
    lines_125M = 125_000_000
    lines_1M = 1_000_000
    
    # Create the subset files
    create_subset_files(input_file, lines_125M, lines_1M, output_file_125M, output_file_1M)

