#!/usr/bin/env python3.6

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", '--fasta_file', type=str, help="input fasta file")
    parser.add_argument("-o", '--output_file', type=str, help="output file")
    parser.add_argument("-d", '--DeepT3_directory', type=str, help="DeepT3 directory")

    args = parser.parse_args()

    from keras.models import model_from_json
    from fasta_reader import readFile
    from helpers import *
    import numpy as np
    from Bio import SeqIO
    import os

    maxlen = 100
    seq_rows, seq_cols = 20, maxlen


    print('Loading model...')

    with open(os.path.join(args.DeepT3_directory, 'models_weights/DeepT3_1.json'), 'r') as f:
        json_string = f.read()

    DeepT3_1 = model_from_json(json_string)
    DeepT3_1.load_weights(os.path.join(args.DeepT3_directory,'models_weights/DeepT3_1.h5'))

    with open(os.path.join(args.DeepT3_directory,'models_weights/DeepT3_2.json'), 'r') as f:
        json_string = f.read()

    DeepT3_2 = model_from_json(json_string)
    DeepT3_2.load_weights(os.path.join(args.DeepT3_directory,'models_weights/DeepT3_2.h5'))

    print('Loading data...')
    testData = readFile(args.fasta_file, maxlen)

    print('Generating features...')
    x_test = createData(testData,"Onehot")
    x_test = x_test.reshape(x_test.shape[0],seq_rows, seq_cols,1)


    print('Predicting...')
    predicted_Probability_1 = DeepT3_1.predict(x_test)
    predicted_Probability_2 = DeepT3_2.predict(x_test)
    prediction_1 = DeepT3_1.predict_classes(x_test)
    prediction_2 = DeepT3_2.predict_classes(x_test)
    prediction = prediction_1 + prediction_2 

    faa_accessions = [i.name for i in SeqIO.parse(args.fasta_file, "fasta")]

    print('Saving the result..')
    f = open (args.output_file, "w")
    for accession, i in zip(faa_accessions, prediction):
        if i == 2:
            f.write("%s\tT3SE\n" % accession)
        else:
            f.write("%s\tnon-T3SE\n" % accession)

    f.close()










