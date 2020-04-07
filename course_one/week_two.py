# An algorithm to return the reverse complment of a string. Note, the string is returned 5' to 3' 
def ReturnCompliment(Pattern):
    rc = ""
    for i in range(len(Pattern)):
        if Pattern[i] == "A":
            rc += "T"
        elif Pattern[i] == "T":
            rc += "A"
        elif Pattern[i] == "C":
            rc += "G"
        else:
            rc += "C"
    rc = rc[::-1]
    return rc

# A function to return the locations in the genome where the skew is a minimum
# The original design had a running time ~ O(3n). This new implementation is closer to ~ O(n)
def minimumSkew(genome):
    # value and location of current global minimum. first index is the value, second is it's position in genome
    minimum_vals_locs = [(0, 0)]
    curr_val = 0
    for i in range(len(genome)):
        if genome[i] == "C":
            curr_val -= 1
            # if we find a new minimum we replace ervything in the list with this new location
            if curr_val < minimum_vals_locs[-1][0]:
                minimum_vals_locs = [(curr_val, i+1)]
        elif genome[i] == "G":
            curr_val += 1
        # here we keeping adding indices as long as  they are equal to the minimum.
        elif curr_val == minimum_vals_locs[-1][0]:
                minimum_vals_locs.append((curr_val, i+1))
    return minimum_vals_locs

# returns the number of mismatches in two DNA strings. Assumes inputs are of same length
def hammingDistance(string1, string2):
    mismatch = 0
    for i in range(len(string1)):
        if string1[i] != string2[i]:
            mismatch += 1
    return mismatch

# given some DNA sequence, genome, and a string, returns all locations of string occurences when hamming distance is at most d
def patternMatches(genome, string, d):
    locations = []
    for i in range((len(genome) - len(string)) + 1):
        text = genome[i:i+len(string)]
        if hammingDistance(text, string) <= d:
            locations.append(i)
    return locations

# returns a list of k-mers with at most d mismatches
def neighbourPatterns(pattern, d):
    neighbours = []
    if d == 0:
        return [pattern]
    if len(pattern) == 1:
        return ["A", "C", "T", "G"]
    suffix_pattern = pattern[1:]
    suffix_neighbours = neighbourPatterns(suffix_pattern, d)
    for string in suffix_neighbours:
        if hammingDistance(pattern[1:], string) < d:
            for nuc in ["A", "C", "T", "G"]:
                neighbours.append(nuc+string)
        else:
            neighbours.append(pattern[0]+string)
    return neighbours

# returns frequent words with mismatches
def FrequentWordsWithMismatches(Text, k, d):
    FrequentWords = {}
    for i in range((len(Text)-k) + 1):
        word = Text[i:i+k]
        neighbours = neighbourPatterns(word, d)
        for kmer in neighbours:
            try:
                FrequentWords[kmer] += 1
            except KeyError:
                FrequentWords[kmer] = 1
    return FrequentWords

# returns frequent words with mismatches
def FrequentWordsWithMismatchesAndReverse(Text, k, d):
    FrequentWords = {}
    for i in range((len(Text)-k) + 1):
        word = Text[i:i+k]
        word_rev = ReturnCompliment(word)
        neighbours = neighbourPatterns(word, d)
        neighbours_rev = neighbourPatterns(word_rev, d)
        neighbours = neighbours + neighbours_rev
        #neighbours = set(neighbours)
        for kmer in neighbours:
            try:
                FrequentWords[kmer] += 1
            except KeyError:
                FrequentWords[kmer] = 1
    return FrequentWords

if __name__ == "__main__":
    # Main task
    out = []
    with open("data/dataset_w2_main.txt") as data:
        data = data.readlines()
        data[0] = data[0].strip("\n")
        data[1] = data[1].strip("\n")
        out = FrequentWordsWithMismatchesAndReverse(data[0], 7, 3)
    m = 0
    k = []
    for key in out:
        if out[key] >= m:
            m = out[key]
    for key in out:
        if out[key] == m:
            k. append(key)
    print(k)

    # Task one
    with open("data/dataset_w2_one.txt") as data:
        data = data.readlines()
        #print(minimumSkew(data[0]))

    # Task two
    with open("data/dataset_w2_two.txt") as data:
        data = data.readlines()
        data[0] = data[0].strip("\n")
        data[1] = data[1].strip("\n")
        #print(hammingDistance(data[0], data[1]))
    
    # Task three
    out = []
    with open("data/dataset_w2_three_ii.txt") as data:
        data = data.readlines()
        data[0] = data[0].strip("\n")
        data[1] = data[1].strip("\n")
        data[2] = data[2].strip("\n")
        out = patternMatches(data[1], data[0], int(data[2]))
        print(len(out))
    
    # Task four
    out = []
    out = neighbourPatterns("TGGAGTTAAC", 2)
    with open("data/output.txt", "w") as output:
        for strings in out:
            output.write(strings)
            output.write(" ")
