import re

import csv


def process_line():
    # Replace quote marks within the line with a space

    f_out = open("annot_modified.txt", "w")
    # Use the csv module to split the line by tabs
    with open("annot.txt", "r") as file:
        for line in file:
            line = line.replace('"', ' ')
            data = line.split("\t")
            gene = data[0]
            annot = data[1]
            annot = annot.rstrip()
            if annot == ("\t"):
                annot = 'NA'
            if annot == "":
                annot = 'NA'
            if len(annot.split()) < 2:
                annot = 'NA'
            #print(annot)
            annot = annot.replace('"', "")
            annot = annot.replace('*', "")
            annot = annot.replace("'", "")
            out_data = "%s\t%s\n" % (gene, annot)
            f_out.write(out_data)




process_line()




def check_quoted_strings(file_path):
    with open(file_path, 'r') as file:
        for line_number, line in enumerate(file, start=1):
            try:
                # Attempt to parse the line
                elements = line.strip().split('\t')
            except ValueError as e:
                # Print lines with improperly formatted quoted strings
                print(f"Line {line_number}: {line.strip()}")

# Check for improperly formatted quoted strings
check_quoted_strings("annot_modified.txt")



def print_problematic_lines(file_path, num_lines=10):
    with open(file_path, 'r') as file:
        for line_number, line in enumerate(file, start=1):
            elements = line.strip().split('\t')
            if len(elements) != 2:
                print(f"Line {line_number}: {line.strip()}")
                num_lines -= 1
                if num_lines == 0:
                    break

# Print the first 10 problematic lines
print_problematic_lines("annot_modified.txt", num_lines=10)






