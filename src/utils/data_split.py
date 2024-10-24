def create_subset_files(input_file, lines_1M, outfile):
    """
    Creates two subset files from the input file: one with 125 million lines and another with 1 million lines.
    
    :param input_file: The path to the large input file
    :param lines_125M: Number of lines for the 125M file (125,000,000)
    :param lines_1M: Number of lines for the 1M file (1,000,000)
    :param output_file_125M: The path for the output file containing 125M lines
    :param output_file_1M: The path for the output file containing 1M lines
    """
    with open (input_file, 'r') as infile:
       with open(output_file, 'w') as outfile:
        for line_number, line in enumerate(infile):
                if line_number < lines_1M:
                    outfile.write(line)
               

if __name__ == "__main__":
    # File paths
    input_file = "/home/afloresep/work/chelombus/data/230M.cxsmiles"  # Input file
    output_file= "/home/afloresep/work/chelombus/data/30M.cxsmiles"  # Output file for 125M lines
    
    # Number of lines to write to each file
    lines_30M = 30_000_000
    create_subset_files(input_file, lines_30M, output_file)

