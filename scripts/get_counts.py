import csv

# Input and output file paths
input_file = "viralcount/telescope-telescope_report.tsv"  # Replace with actual path
output_file = "viralcount/counts.tsv"  # Replace with actual path

# Open input and output files
with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
    # Initialize CSV reader and writer
    reader = csv.reader(infile, delimiter='\t')
    writer = csv.writer(outfile, delimiter='\t')

    # Skip the first two lines (headers and comments)
    next(reader)  # Skip "RunInfo" line
    next(reader)  # Skip the "transcript" header line

    # Process each row
    for row in reader:
        # Remove quotes from the first column and extract the first and third columns
        transcript = row[0].replace('"', '')  # Removing double quotes from the transcript name
        final_count = row[2]  # Extract the third column (final_count)

        # Write the cleaned-up data to the output file
        writer.writerow([transcript, final_count])

print(f"Processing complete. Results saved to {output_file}")