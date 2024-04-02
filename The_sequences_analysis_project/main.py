import requests
import random
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from typing import Dict


def parsing(url: str):
    """Obtaining nucleotide sequences from a database epd.expasy.org"""
    st_accept = 'text/html'
    st_useragent = ('Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/118.0.0.0 '
                    'YaBrowser/23.11.0.0 Safari/537.36')
    headers = {
        'Accept': st_accept,
        'User-Agent': st_useragent
    }

    return requests.get(url, headers)


def clean_sequences(file_name: str):
    """Creating clean sequences without unnecessary information"""
    with open(file_name, 'r') as f:
        lines = f.readlines()

    cleaned_lines = []
    sequence = ''

    for line in lines:
        if line.startswith('>'):
            if sequence and ('N' not in sequence):
                cleaned_lines.append(sequence)
            sequence = ''
            # sequence += line.strip()
        else:
            sequence += line.strip()

    if sequence:
        cleaned_lines.append(sequence)

    with open('data/clean_sequences.txt', 'w') as f:
        for line in cleaned_lines:
            f.write(line + '\n')


def create_reverse_complementary_sequences(clean_sequensec_file: str):
    """Creation of back-complementary nucleotide sequences"""
    with open(clean_sequensec_file, 'r') as f:
        lines = f.readlines()
        reversed_lines = []
        complementary_line = ''
        for line in lines:
            for i in line:
                if i == 'G':
                    complementary_line += 'C'
                elif i == 'A':
                    complementary_line += 'T'
                elif i == 'C':
                    complementary_line += 'G'
                elif i == 'T':
                    complementary_line += 'A'
            reversed_lines.append(complementary_line[::-1])
            complementary_line = ''

    with open('data/reverse_complementary_sequences.txt', 'w+') as f:
        for sequence in reversed_lines:
            f.write(sequence + '\n')


def create_random_sequences(str_quantity: int):
    """Creating random sequences of nucleotides"""

    with open('data/random_sequences.txt', 'w+') as f:
        nucleotides = ['A', 'T', 'G', 'C']
        sequence = ''
        for string in range(str_quantity):
            for nucl in range(81):
                sequence += random.choice(nucleotides)

            f.write(sequence + '\n')
            sequence = ''


def calculate_dinucleotide_properties(sequences_file: str) -> dict:
    """
    Calculates dinucleotide properties based on a given sequences file.

    Args:
        sequences_file (str): Path to the input sequences file.

    Returns:
        dict: A dictionary containing dinucleotide properties.
    """

    df = pd.read_csv('data/Dinucleotide Property Database.txt', sep='\t')
    properties = ['stacking_energy', 'mobility', 'roll', 'slide', 'slide_stiffness', 'roll_stiffness']
    dictionaries = {prop: {i: 0 for i in range(-50, 31) if i != 0} for prop in properties}

    with open(sequences_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            n = -50
            for nuleotide in range(len(line) - 1):
                dinucleotide = line[nuleotide] + line[nuleotide + 1]
                for i, prop in enumerate(properties):
                    try:
                        dictionaries[prop][n] += df.at[i, dinucleotide]
                    except KeyError:
                        n += 1 if n != -1 else 2
                n += 1 if n != -1 else 2

    total_sequences = len(lines)
    for dictionary in dictionaries.values():
        for key in dictionary:
            dictionary[key] /= total_sequences

    return dictionaries


def calculate_tetranucleotide_properties(sequences_file: str) -> dict:
    """
    Calculates tetranucleotide properties based on the given sequence file.

    Args:
        sequences_file (str): Path to the file containing sequences.

    Returns:
        dict: A dictionary containing tetranucleotide properties.
    """

    df = pd.read_csv('data/Table of tetranucleotides.csv', sep=';')

    dictionary_of_param = {i: 0 for i in range(-50, 29) if i != 0}

    with open(sequences_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            n = -50
            for tetranucleotide in range(len(line) - 4):
                substring = line[tetranucleotide:tetranucleotide + 4]
                try:
                    dictionary_of_param[n] += df.at[1, substring]
                    n += 1 if n != -1 else 2
                except KeyError:
                    n += 1 if n != -1 else 2

    total_sequences = len(lines)
    for key in dictionary_of_param:
        dictionary_of_param[key] /= total_sequences

    return dictionary_of_param


def calculate_ultrasonic_dinucleotide_properties(sequences_file: str) -> dict:
    """
    Calculates ultrasonic dinucleotide properties based on the given sequence file.

    Args:
        sequences_file (str): Path to the file containing sequences.

    Returns:
        dict: A dictionary containing ultrasonic dinucleotide properties.
    """

    df = pd.read_csv('data/Ultrasonic_dinucleotides.csv', sep=';')
    print(df.head(5))

    dictionary_of_param = {i: 0 for i in range(-50, 31) if i != 0}

    with open(sequences_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            n = -50
            for tetranucleotide in range(len(line) - 2):
                substring = line[tetranucleotide:tetranucleotide + 2]
                try:
                    dictionary_of_param[n] += df.at[1, substring]
                    n += 1 if n != -1 else 2
                except KeyError:
                    n += 1 if n != -1 else 2

    total_sequences = len(lines)
    for key in dictionary_of_param:
        dictionary_of_param[key] /= total_sequences

    return dictionary_of_param


def calculate_hexanucleotide_properties(sequences_file: str) -> Dict[int, float]:
    """
    Calculates hexanucleotide properties based on a given sequences file.

    Args:
        sequences_file (str): Path to the sequences file.

    Returns:
        Dict[int, float]: A dictionary containing hexanucleotide properties.
    """

    # Read the hexanucleotide properties table
    df = pd.read_csv('data/Table of hexanucleotides.csv', sep=';')

    # Initialize a dictionary to store the properties
    dictionary_of_param = {i: 0 for i in range(-50, 27) if i != 0}

    # Read sequences from the file
    with open(sequences_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            n = -50
            for tetranucleotide in range(len(line) - 6):
                substring = line[tetranucleotide:tetranucleotide + 6]
                try:
                    dictionary_of_param[n] += df.at[0, substring]
                    n += 1 if n != -1 else 2
                except KeyError:
                    n += 1 if n != -1 else 2

    # Normalize the properties
    total_sequences = len(lines)
    for key in dictionary_of_param:
        dictionary_of_param[key] /= total_sequences

    return dictionary_of_param


def create_graph(original_data, random_data, title, save_path, reverse_data=None):
    """
    Creates a graph based on data from two or three dictionaries.

    Args:
        original_data (dict): A dictionary with data for the original sequence.
        random_data (dict): A dictionary with data for a random sequence.
        title (str): The title of the graph.
        save_path (str): The path to save the graph.
        reverse_data (dict): A dictionary with data for a back-complementary sequence.
    """
    plt.figure(figsize=(40, 7))
    plt.title(title, fontweight='bold', fontsize=25, pad=10)
    plt.xticks(np.arange(-80, 30, 1))
    plt.xlabel('Position relative to TSS', fontsize=16, labelpad=10, color='#B22222')
    plt.ylabel('Average parameter value', fontsize=16, labelpad=10, color='#B22222')
    plt.grid()
    plt.plot(original_data.keys(), original_data.values(), label='Original sequences')
    if reverse_data:
        plt.plot(reverse_data.keys(), reverse_data.values(), label='Reverse complementary sequences')
    plt.plot(random_data.keys(), random_data.values(), label='Random sequences')

    plt.legend()
    plt.savefig(save_path)


# req = parsing('https://epd.expasy.org/epd/wwwtmp/hg38_C2y3h.fa')
#
# with open('data/dirty_sequences.txt', 'w+') as f:
#     f.write(req.text)

clean_sequences('data/dirty_sequences.txt')
create_reverse_complementary_sequences('data/clean_sequences.txt')
create_random_sequences(29597)

original_hex = calculate_hexanucleotide_properties('data/clean_sequences.txt')
reverse_hex = calculate_hexanucleotide_properties('data/reverse_complementary_sequences.txt')
random_hex = calculate_hexanucleotide_properties('data/random_sequences.txt')
create_graph(original_hex, random_hex,
             'Comparison of average values for high levels of hexanucleotide methylation',
             'graphics/graphic_of_hexanucleotides.png', reverse_hex)

original_ultrasonic = calculate_ultrasonic_dinucleotide_properties('data/clean_sequences.txt')
reverse_ultrasonic = calculate_ultrasonic_dinucleotide_properties('data/reverse_complementary_sequences.txt')
random_ultrasonic = calculate_ultrasonic_dinucleotide_properties('data/random_sequences.txt')
create_graph(original_ultrasonic, random_ultrasonic,
             'Comparison of Average Values for Ultrasonic Cleavage of Dinucleotides',
             'graphics/graphic_of_ultrasonic_dinucleotides.png', reverse_ultrasonic)

original_tetra = calculate_tetranucleotide_properties('data/clean_sequences.txt')
reverse_tetra = calculate_tetranucleotide_properties('data/reverse_complementary_sequences.txt')
random_tetra = calculate_tetranucleotide_properties('data/random_sequences.txt')
create_graph(original_tetra, random_tetra,
             'Comparison of Average Values for Ultrasonic Cleavage of Tetranucleotides',
             'graphics/graphic_of_tetranucleotides.png', reverse_tetra)

dict_of_dinucl_original = calculate_dinucleotide_properties('data/clean_sequences.txt')
dict_of_dinucl_random = calculate_dinucleotide_properties('data/random_sequences.txt')

dict_list = [dict_of_dinucl_original, dict_of_dinucl_random]

for prop, param in dict_of_dinucl_original.items():
    create_graph(dict_list[0][prop], dict_list[1][prop],
                 f'Comparison of Average "{prop}" property', f'graphics/graphic_of_average_{prop}_property.png')
