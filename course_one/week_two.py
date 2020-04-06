from course_one.week_one import *
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

# A more efficient version of the FrequentWords algorithm in O(n)
def FrequentWords_faster(Text, k):
    FrequentWords = {}
    for i in range((len(Text)-k) + 1):
        word = Text[i:i+k]
        try:
            FrequentWords[word] += 1
        except KeyError:
            FrequentWords[word] = 1
    return FrequentWords

if __name__ == "__main__":
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
    