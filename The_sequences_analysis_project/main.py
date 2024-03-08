import requests
import random


def parsing(url: str):
    '''Obtaining nucleotide sequences from a database epd.expasy.org'''
    st_accept = 'text/html'
    st_useragent = ('Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/118.0.0.0 '
                    'YaBrowser/23.11.0.0 Safari/537.36')
    headers = {
        'Accept': st_accept,
        'User-Agent': st_useragent
    }

    return requests.get(url, headers)


def clean_sequences(file_name):
    '''Creating clean sequences without unnecessary information'''
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
    '''Creation of back-complementary nucleotide sequences'''
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
    '''Creating random sequences of nucleotides'''
    with open('random_sequences.txt', 'w+') as f:
        nucleotides = ['A', 'T', 'G', 'C']
        sequence = ''
        for str in range(str_quantity):
            for nucl in range(81):
                sequence += random.choice(nucleotides)

            f.write(sequence + '\n')
            sequence = ''

req = parsing('https://epd.expasy.org/epd/wwwtmp/hg38_C2y3h.fa')

with open('dirty_sequences.txt', 'w+') as f:
    f.write(req.text)

clean_sequences('dirty_sequences.txt')
create_reverse_complementary_sequences('clean_sequences.txt')
create_random_sequences(29598)
