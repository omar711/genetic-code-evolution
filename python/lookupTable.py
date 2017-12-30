
# Get the genetic code numbers
def getGeneticCodeIDs():
    return [ 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23 ]

# Get the alphabet - array of characters.
def getNucleotideAlphabet():
    return ['T', 'C', 'A', 'G']

# Get the amino acid alphabet
def getAminoAlphabet():
    return ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

# Make a list of all possible codons
def getCodons():
    alphabet = getNucleotideAlphabet()
    length = len(alphabet)
    codons = []
    for i in range(length):
        for j in range(length):
            for k in range(length):
                codons.append(str(alphabet[i]) + str(alphabet[j] + str(alphabet[k])))

    return codons

if __name__ == '__main__':
    print(getGeneticCodeIDs())
    print(getNucleotideAlphabet())
    print(getAminoAlphabet())
    print(getCodons())
