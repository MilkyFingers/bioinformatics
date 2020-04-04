import time

# Takes a sequence of DNA and returns the number of occurences of a specified k-mer (pattern)
def PatternCount(Text, Pattern):
	count = 0
	for i in range((len(Text) - len(Pattern)) + 1):
		if Text[i:i+len(Pattern)] == Pattern:
			count += 1
	return count

# A O(n^2) algorithm for finding if there exists any k-mers of a specified length within a DNA sequence
def FrequentWords(Text, k):
    FrequentWords = []
    count = []
    for i in range((len(Text) - k) + 1):
        currentPat = Text[i:i+k]
        count.append(PatternCount(Text, currentPat))
    max_count = max(count)
    for i in range(len(count)):
        if count[i] == max_count:
            FrequentWords.append(Text[i:i+k])
    FrequentWords = set(FrequentWords)
    return FrequentWords

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

# Returns all the indices where a substring accurs in a given string
def PatternMatching(Input, Pattern):
    locations = []
    add = 0
    try:
        while True:
            x = Input.index(Pattern)
            locations.append(x + add)
            new_start_loc = x + 1
            add += new_start_loc
            Input = Input[new_start_loc::]
    except:
        return locations

# Turns some symbol into a number
def SymbolToNumber(Symbol):
    if Symbol == "A":
        return 0
    elif Symbol == "C":
        return 1
    elif Symbol == "G":
        return 2
    else:
        return 3

# Turns a number into a given symbol
def NumberToSymbol(Number):
    if Number == 0:
        return "A"
    elif Number == 1:
        return "C"
    elif Number == 2:
        return "G"
    else:
        return "T"

# A recursive algorithm to turn some pattern into a number
def PatternToNumber(Pattern):
    if Pattern == "":
        return 0
    else:
        symbol = Pattern[-1]
        prefix = Pattern[0:-1]
        return 4 * PatternToNumber(prefix) + SymbolToNumber(symbol)

# A recursive algorithm to find to pattern associated with a given number
def NumberToPattern(number, k):
    if k == 1:
        return NumberToSymbol(number)
    else:
        prefixNum = divmod(number,4)[0]
        r = divmod(number,4)[1]
        symbol = NumberToSymbol(r)
        return NumberToPattern(prefixNum, k-1) + symbol

# A function that returns a frequency array instead of dictionary
def ComputingFrequencies(Text, k):
    frequency_array = []
    for i in range(4 ** k):
        frequency_array.append(0)
    for i in range((len(Text) - k) + 1):
        pat = Text[i:i+k]
        j = PatternToNumber(pat)
        frequency_array[j] += 1
    return frequency_array

# Finding all (L, t) clumps in a genome of some k-mer
# Genone is the entire dataset, L is the window size (typically 500), t is the threshold for the number of times we need something to appear and k is the k-mer length
def ClumpFinding(Genome, L, t, k):
    frequent_patterns = [] # we will add the common k-mers to this
    clumps = [] # this will correspond to the k-mer patters. 1 for yes, they surpass t and 0 for no
    for i in range(4 ** k):
        clumps.append(0)
    Text = Genome[0:L]
    frequency_array = ComputingFrequencies(Text, k)
    for i in range(4 ** k):
        if frequency_array[i] >= t:
            clumps[i] = 1
    for i in range(1, (len(Genome) - L) + 1):
        first_pat = Genome[i - 1:(i - 1) + k]
        number = PatternToNumber(first_pat)
        frequency_array[number] -= 1
        last_pat = Genome[i + (L - k):(i + (L - k)) + k]
        number = PatternToNumber(last_pat)
        frequency_array[number] += 1
        if frequency_array[i] >= t:
            clumps[i] = 1
    for i in range(4 ** k):
        if clumps[i] == 1:
            pat = NumberToPattern(i, k)
            frequent_patterns.append(pat)
    return frequent_patterns

# A clump finding algorithm implemented using dictionaries. This removes the need to store k-mers never seen and is more simple. The course directors need a hand...
def ClumpFindingBetter(Genome, L, t, k):
    frequent_words = [] # the final output will be stored here
    Text = Genome[0:L]
    frequency_dict = FrequentWords_faster(Text, k) # a dictionary of all frequencies of k-mers in the first window
    for kmer in frequency_dict:
        if frequency_dict[kmer] >= t:
            frequent_words.append(kmer)
    # now we loop over the remaining windows, removing the last word, adding the next, or creating it's key if it is new
    for i in range(1, (len(Genome) - L) + 1):
        removed_pat = Genome[i - 1:(i - 1) + k] # the last k-mer, now out of the window
        new_pat = Genome[i + (L - k):(i + (L - k)) + k] # the new k-mer just seen
        frequency_dict[removed_pat] -= 1 # removes one
        try:
            frequency_dict[new_pat] += 1
        except KeyError:
            frequency_dict[new_pat] = 1 # accounting for the k-mer if it hasnt been seen yet
        if frequency_dict[removed_pat] >= t:
            frequent_words.append(removed_pat)
        if frequency_dict[new_pat] >= t:
            frequent_words.append(new_pat)
    return frequent_words


if __name__ == "__main__":

    # Task one
    with open("data/dataset_cc_one.txt") as data:
        data = data.readlines()
        data[0] = data[0].strip("\n")
        data[1] = data[1].strip("\n")
        print(PatternCount(data[0], data[1]))

    # Task two
    with open("data/dataset_cc_two.txt") as data:
        data = data.readlines()
        data[0] = data[0].strip("\n")
        data[1] = data[1].strip("\n")
        start = time.time()
        print(FrequentWords(data[0],int(data[1])))
        end = time.time()
        print(end-start)

        start = time.time()
        print(FrequentWords_faster(data[0],int(data[1])))
        end = time.time()
        print(end-start)

        print("---------")

     # Task three
    with open("data/dataset_cc_three.txt") as data:
        data = data.readlines()
        data[0] = data[0].strip("\n")
        print(ReturnCompliment(data[0]))

    # Task four
    x = []
    with open("data/vibrio_cholerae.txt") as data:
        data = data.readlines()
        data[0] = data[0].strip("\n")
        x = PatternMatching(data[0], "ATGATCAAG")
    with open("data/vibrio_locations.txt", "w") as output:
        for i in range(len(x)):
            output.write(str(x[i]))
            output.write(" ")



