import os
import pandas as pd
from wordcloud import WordCloud
import matplotlib.pyplot as plt

def generate_word_cloud(files_path, output_path, ignore_list):
    # Create an empty string to store all GO terms
    all_go_terms = ""

    # Loop through each file in the specified folder
    for file_name in os.listdir(files_path):
        if file_name.endswith("enriched_P0.05.txt"):
            file_path = os.path.join(files_path, file_name)

            # Read the file using pandas
            df = pd.read_csv(file_path, sep='\t')

            # Check if the DataFrame has any rows
            if not df.empty and 'go_term' in df:
                # Extract and clean GO terms
                go_terms = df['go_term'].astype(str).str.replace('|', '').str.replace(',', '').str.replace('\s+', ' ')

                # Remove substrings specified in the ignore_list
                for ignore in ignore_list:
                    go_terms = go_terms.str.replace(ignore, '')

                # Remove extra spaces and concatenate to the string
                filtered_terms = [term.strip() for term in go_terms if term.strip()]
                all_go_terms += " ".join(filtered_terms) + " "

    # Check if there are any GO terms
    if not all_go_terms:
        print("No GO terms found. Exiting.")
        return

    # Generate word cloud
    wordcloud = WordCloud(width=800, height=400, background_color='white').generate(all_go_terms)

    # Display the generated word cloud using matplotlib
    plt.figure(figsize=(10, 5))
    plt.imshow(wordcloud, interpolation='bilinear')
    plt.axis('off')

    # Save the plot as a PDF file
    plt.savefig(output_path, format='pdf', bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    # Get the current working directory
    current_directory = os.getcwd()

    # Set the output path for the PDF file
    output_file_path = os.path.join(current_directory, 'GO_wordcloud.pdf')

    # Specify the list of substrings to ignore
    ignore_list = ['BP', 'Process', 'metabolic', 'CC', 'none', 'MF', 'activity']

    # Use the current working directory as the folder path
    generate_word_cloud(current_directory, output_file_path, ignore_list)
