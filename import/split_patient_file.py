import argparse
import math
import pandas as pd
import os

def wc(F):
    with open(F, 'rt') as inf:
        for i, l in enumerate(inf):
            pass
    return i + 1

def split_df(df, number_of_chunks):
    chunks = []
    chunk_size = round(len(df) / number_of_chunks)
    for i in range(number_of_chunks):
        if i == number_of_chunks - 1:
            yield df[i*chunk_size:]
        else:
            yield df[i*chunk_size:(i+1)*chunk_size]

def get_digits_for_output(N):
    if N == 1:
        return 1
    if N < 1:
        raise ValueError(f"N = {N}: Number of chunks has to be greater than 0")
    return int(math.log(N - 1, 10) + 1)

def main(args):
    df = pd.read_csv(args.input, sep='\t')

    ind = 0
    number_of_digits = get_digits_for_output(args.number_of_chunks)
    print(number_of_digits)
    for i in split_df(df, args.number_of_chunks):
        outf_name = f"{os.path.splitext(args.input)[0]}_{str(ind).zfill(number_of_digits)}.tsv"
        pd
        i.to_csv(outf_name, sep="\t", index=None)
        ind += 1

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', dest='input', help='the patient data file')
    parser.add_argument('--number_of_chunks', dest='number_of_chunks', help='number of chunks', type=int)
    args = parser.parse_args()
    main(args)
