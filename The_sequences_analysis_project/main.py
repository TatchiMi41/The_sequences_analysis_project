import requests
import random
import pandas as pd
import matplotlib.pyplot as plt


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


def clean_sequences(file_name):
    """Creating clean sequences without unnecessary information"""
    with open(file_name, 'r') as f:
        lines = f.readlines()

    cleaned_lines = []
    sequence = ''

    for line in lines:
        if line.startswith('>'):
            if sequence:
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


def create_reverse_complementary_sequences(clean_sequensec_file):
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


def create_random_sequences(str_quantity):
    """Creating random sequences of nucleotides"""
    with open('data/random_sequences.txt', 'w+') as f:
        nucleotides = ['A', 'T', 'G', 'C']
        sequence = ''
        for string in range(str_quantity):
            for nucl in range(81):
                sequence += random.choice(nucleotides)

            f.write(sequence + '\n')
            sequence = ''


def calculate_data_dinucleotide_property_database(sequences_file):
    df = pd.read_csv('data/Dinucleotide Property Database.txt', sep='\t')
    print(df.head(7))

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

    for dictionary in dictionaries.values():
        for key in dictionary:
            dictionary[key] /= 29598

    for prop, dictionary in dictionaries.items():
        plt.figure(figsize=(10, 5))
        plt.plot(dictionary.keys(), dictionary.values())
        plt.savefig(f'graphics/{prop}.png')


def calculate_data_tetranucleotide_property_database(sequences_file):
    df = pd.read_csv('data/Table of tetranucleotides.csv', sep=';')
    print(df.head(5))

    dictionary_of_param = {i:0 for i in range(-50, 29) if i != 0}

    with open(sequences_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            n = -50
            for tetranucleotide in range(len(line) - 4):
                substring = line[tetranucleotide:tetranucleotide+4]
                try:
                    dictionary_of_param[n] += df.at[1, substring]
                    n += 1 if n != -1 else 2
                except KeyError:
                    n += 1 if n != -1 else 2

    for key in dictionary_of_param:
        dictionary_of_param[key] /= 29598

    plt.figure(figsize=(10, 5))
    plt.plot(dictionary_of_param.keys(), dictionary_of_param.values())
    plt.savefig(f'graphics/graphic_of_tetranucleotides.png')


def calculate_data_ultrasonic_dinucleotides_property_database(sequences_file):
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

    for key in dictionary_of_param:
        dictionary_of_param[key] /= 29598

    plt.figure(figsize=(10, 5))
    plt.plot(dictionary_of_param.keys(), dictionary_of_param.values())
    plt.savefig(f'graphics/graphic_of_ultrasonic_dinucleotides.png')
# req = parsing('https://epd.expasy.org/epd/wwwtmp/hg38_C2y3h.fa')
#
# with open('dirty_sequences.txt', 'w+') as f:
#     f.write(req.text)

# clean_sequences('dirty_sequences.txt')
# create_reverse_complementary_sequences('clean_sequences.txt')
# create_random_sequences(29598)
# calculate_data_dinucleotide_property_database('clean_sequences.txt')

# calculate_data_tetranucleotide_property_database('clean_sequences.txt')
calculate_data_ultrasonic_dinucleotides_property_database('clean_sequences.txt')