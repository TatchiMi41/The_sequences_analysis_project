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

    with open('clean_sequences.txt', 'w') as f:
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

    with open('reverse_complementary_sequences.txt', 'w+') as f:
        for sequence in reversed_lines:
            f.write(sequence + '\n')


def create_random_sequences(str_quantity):
    """Creating random sequences of nucleotides"""
    with open('random_sequences.txt', 'w+') as f:
        nucleotides = ['A', 'T', 'G', 'C']
        sequence = ''
        for string in range(str_quantity):
            for nucl in range(81):
                sequence += random.choice(nucleotides)

            f.write(sequence + '\n')
            sequence = ''


def calculate_data_dinucleotide_property_database(sequences_file):
    df = pd.read_csv('Dinucleotide Property Database.txt', sep='\t')
    print(df.head(7))
    with open(sequences_file, 'r') as f:
        lines = f.readlines()
        dictionary_stacking_energy = {i: 0 for i in range(-50, 31)}
        dictionary_mobility = {i: 0 for i in range(-50, 31)}
        dictionary_roll = {i: 0 for i in range(-50, 31)}
        dictionary_slide = {i: 0 for i in range(-50, 31)}
        dictionary_slide_stiffness = {i: 0 for i in range(-50, 31)}
        dictionary_roll_stiffness = {i: 0 for i in range(-50, 31)}
        for line in lines:
            n = -50
            for nuleotide in range(len(line) - 1):
                dinucleotide = ''
                dinucleotide += line[nuleotide] + line[nuleotide + 1]
                try:
                    dictionary_stacking_energy[n] += df.at[0, dinucleotide]
                    dictionary_mobility[n] += df.at[1, dinucleotide]
                    dictionary_roll[n] += df.at[2, dinucleotide]
                    dictionary_slide[n] += df.at[3, dinucleotide]
                    dictionary_slide_stiffness[n] += df.at[4, dinucleotide]
                    dictionary_roll_stiffness[n] += df.at[5, dinucleotide]
                except KeyError:
                    if n != -1:
                        n += 1
                    else:
                        n += 2
                if n != -1:
                    n += 1
                else:
                    n += 2

        for key in dictionary_stacking_energy:
            dictionary_stacking_energy[key] /= 29598

        for key in dictionary_mobility:
            dictionary_mobility[key] /= 29598

        for key in dictionary_slide:
            dictionary_slide[key] /= 29598

        for key in dictionary_roll:
            dictionary_roll[key] /= 29598

        for key in dictionary_roll_stiffness:
            dictionary_roll_stiffness[key] /= 29598

        for key in dictionary_slide_stiffness:
            dictionary_slide_stiffness[key] /= 29598

        plt.figure(figsize=(10, 5))
        plt.plot(dictionary_stacking_energy.keys(), dictionary_stacking_energy.values())
        plt.savefig('graphics/stacking_energy.png')

        plt.figure(figsize=(10, 5))
        plt.plot(dictionary_mobility.keys(), dictionary_mobility.values())
        plt.savefig('graphics/mobility.png')

        plt.figure(figsize=(10, 5))
        plt.plot(dictionary_slide.keys(), dictionary_slide.values())
        plt.savefig('graphics/slide.png')

        plt.figure(figsize=(10, 5))
        plt.plot(dictionary_roll.keys(), dictionary_roll.values())
        plt.savefig('graphics/roll.png')

        plt.figure(figsize=(10, 5))
        plt.plot(dictionary_roll_stiffness.keys(), dictionary_roll_stiffness.values())
        plt.savefig('graphics/roll_stiffness.png')

        plt.figure(figsize=(10, 5))
        plt.plot(dictionary_slide_stiffness.keys(), dictionary_slide_stiffness.values())
        plt.savefig('graphics/slide_stiffness.png')





# req = parsing('https://epd.expasy.org/epd/wwwtmp/hg38_C2y3h.fa')
#
# with open('dirty_sequences.txt', 'w+') as f:
#     f.write(req.text)

clean_sequences('dirty_sequences.txt')
# create_reverse_complementary_sequences('clean_sequences.txt')
# create_random_sequences(29598)
calculate_data_dinucleotide_property_database('clean_sequences.txt')
