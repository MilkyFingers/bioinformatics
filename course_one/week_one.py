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
        print(FrequentWords(data[0],int(data[1])))

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



