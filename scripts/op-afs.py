#!/usr/bin/env python3
import pandas as pd
import sys
import matplotlib.pyplot as plt


def read_file_to_matrix(file_path):
    try:
        df = pd.read_csv(file_path, sep='\t') #, index_col=0)
        return df
    except Exception as e:
        print(f"Failed to read file {file_path}: {e}")
        sys.exit(1)

def are_elements_equal(vector):
    """
    Check if all elements in the vector are equal.
    Args: - vector: A list or vector of elements.
    Returns:- True if all elements are equal, False otherwise.
    """
    if not vector:  # Check if the vector is empty
        return True  # Can decide based on context; here, an empty vector is considered to have "equal" elements
    # Check if all elements in the vector are equal to the first element
    return all(element == vector[0] for element in vector)

def allele_freq(values, label):
    """
    Process a list of values to count occurrences of each unique value.does not consider lists in which all elements are equal 
    Args: - values: A list of values; - label: A label associated with the values.
    Returns: - None, but prints counts of each unique value prefixed by the label.
    """
    if not are_elements_equal(values): 
        total_values = len(values)
        value_counts = {}
        for value in values:
            if value not in value_counts:
                value_counts[value] = 1
            else:
                value_counts[value] += 1
        # Format the output
        for value, count in value_counts.items():
            #print(f"{label} {value} {count}")
            freq=count/total_values
            return label, value, count, freq 
    else: pass 

def plot_histogram(vector, bins='auto', title='Histogram', xlabel='Values', ylabel='Frequency',  save_path='histogram.png'):
    """
    Plots a histogram from a vector of values.

    Args:
    - vector: A list or array of values to plot.
    - bins: The number of bins or the binning strategy (default 'auto').
    - title: Title of the histogram (default 'Histogram').
    - xlabel: Label for the X-axis (default 'Values').
    - ylabel: Label for the Y-axis (default 'Frequency').
    """
    plt.figure(figsize=(8, 6))  # Optional: Adjust figure size
    plt.hist(vector, bins=bins, color='skyblue', edgecolor='black', alpha=0.7)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(axis='y', alpha=0.75)  # Optional: Add grid lines for better readability
  # Save the figure to a PNG file
    plt.savefig(save_path)
    plt.close()  # Close the figure to free memory

def plot_histograms_from_dict(data_dict, bins='auto', title='Histograms from Vectors', xlabel='Values', ylabel='Frequency', save_path='histograms.png'):
    """
    Plots multiple histograms from vectors stored in a dictionary. Each vector is plotted in a separate panel.
    
    Args:
    - data_dict: A dictionary where each key-value pair corresponds to a label and its associated vector of values.
    - bins: The number of bins or the binning strategy for the histograms (default 'auto').
    - title: Overall title for the plot (default 'Histograms from Vectors').
    - xlabel: Label for the X-axis (default 'Values').
    - ylabel: Label for the Y-axis (default 'Frequency').
    - save_path: Path where the plot will be saved as a PNG file (default 'histograms.png').
    """
    num_vectors = len(data_dict)
    fig, axs = plt.subplots(num_vectors, 1, figsize=(8, 6 * num_vectors))  # Adjust the figure size as needed
    
    if num_vectors == 1:
        axs = [axs]  # Make sure axs is iterable for a single subplot
    
    for ax, (label, vector) in zip(axs, data_dict.items()):
        ax.hist(vector, bins=bins, color='skyblue', edgecolor='black', alpha=0.7)
        ax.set_title(f'{label}')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(axis='y', alpha=0.75)
    
    plt.tight_layout(pad=3.0)
    plt.suptitle(title, va='top')
    plt.subplots_adjust(top=0.95)  # Adjust top spacing to accommodate the main title
    plt.savefig(save_path)


def main():
    # Check for command-line argument for file path
    if len(sys.argv) < 2:
        print("Usage: script.py <file_path>")
        sys.exit(1)
        pass

    file_path = sys.argv[1]
    matrix_df = read_file_to_matrix(file_path)

   # print (matrix_df.info())
    # Process each vector (column) in the matrix starting from the 4th column
    counts_d={}; counts_f={}
    for column in matrix_df.columns[3:]:
        label=column
        values = matrix_df[column].iloc[0:].tolist()  # Rest of the column values
        label, value, count, freq= allele_freq(values, label)
        counts_d.setdefault(value, []).append(count)
        counts_f.setdefault(value, []).append(freq)
        #print (label, value, count, freq)
        #print (counts_d)

    plot_histograms_from_dict(counts_d, bins='auto', title='Histograms from Vectors', xlabel='Values', ylabel='Frequency', save_path='counts.png')
    plot_histograms_from_dict(counts_f, bins='auto', title='Histograms from Vectors', xlabel='Values', ylabel='Frequency', save_path='freqs.png')


if __name__ == "__main__":
    main()
