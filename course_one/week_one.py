import time

# Takes a sequence of DNA and returns the number of occurences of a specified k-mer (pattern)
def PatternCount(Text, Pattern):
	count = 0
	for i in range(len(Text) - len(Pattern)):
		if Text[i:i+len(Pattern)] == Pattern:
			count += 1
	return count

# A O(n^2) algorithm for finding if there exists any k-mers of a specified length within a DNA sequence
def FrequentWords(Text, k):
    FrequentWords = []
    count = []
    for i in range(len(Text) - k):
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
    for i in range(len(Text)-k):
        word = Text[i:i+k]
        try:
            FrequentWords[word] += 1
        except KeyError:
            FrequentWords[word] = 0
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
def Find_substring_locations(Input, Pattern):
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
        x = Find_substring_locations(data[0], "ATGATCAAG")
    with open("data/vibrio_locations.txt", "w") as output:
        for i in range(len(x)):
            output.write(str(x[i]))
            output.write(" ")


