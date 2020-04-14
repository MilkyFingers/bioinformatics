#from course_one.week_two import *
import itertools

def hammingDistance(string1, string2):
    mismatch = 0
    for i in range(len(string1)):
        if string1[i] != string2[i]:
            mismatch += 1
    return mismatch

"""
A poor brute-force solution to finding motifs in a set of Dna strings
"""
def motifEnumeration(dna, k, d):
    motifs = []
    for i in range((len(dna[0]) - k) + 1):
        kmer = dna[0][i:i+k]
        words = neighbourPatterns(kmer, d)
        for w in words:
            wordsPrime = neighbourPatterns(w, d)
            common = 0
            for string in dna:
                for wPrime in wordsPrime:
                    if wPrime in string:
                        common += 1
                        break
            if common == len(dna):
                motifs.append(w)
    return set(motifs)

# subroutine for medianString
def distanceP(pattern, dna):
    totalDistance = 0
    for string in dna:
        localMinDistance = float("inf")
        for i in range((len(string) - len(pattern)) + 1):
            compare = string[i:i+len(pattern)]
            if hammingDistance(compare, pattern) < localMinDistance:
                localMinDistance = hammingDistance(compare, pattern)
        totalDistance += localMinDistance
    return totalDistance

# Yet another brute-force...this course is slow
def medianString(dna, k):
    distanceString = float("inf")
    median = ""
    for combo in itertools.product("ACGT", repeat = k):
        kmer = ""
        for nuc in combo:
            kmer += nuc
        if distanceP(kmer, dna) < distanceString:
            distanceString = distanceP(kmer, dna)
            median = kmer
    return median

# profile-most probable kmer
def profileMostProbable(dna, k, profile):
    mostProbable = ""
    probability = float("-inf")
    for i in range((len(dna) - k) + 1):
        kmer = dna[i:i+k]
        profileVal = 1
        for i in range(k):
            if kmer[i] == "A":
                profileVal *= profile[0][i]
            elif kmer[i] == "C":
                profileVal *= profile[1][i]
            elif kmer[i] == "G":
                profileVal *= profile[2][i]
            else:
                profileVal *= profile[3][i]
        if profileVal > probability:
            probability = profileVal
            mostProbable = kmer
    return mostProbable

#subroutine to genrate a profile from list of kmers
def generateProfile(kmers):
    profile = []
    for i in range(4):
        app = []
        for i in range(len(kmers[0])):
            app.append(0.0)
        profile.append(app)
    for j in range(len(kmers[0])):
        column = ""
        for kmer in kmers:
            column += kmer[j]
        a = len([i for i in range(len(column)) if column.startswith("A", i)])/len(column)
        c = len([i for i in range(len(column)) if column.startswith("C", i)])/len(column)
        g = len([i for i in range(len(column)) if column.startswith("G", i)])/len(column)
        t = len([i for i in range(len(column)) if column.startswith("T", i)])/len(column)
        profile[0][j] = a
        profile[1][j] = c
        profile[2][j] = g
        profile[3][j] = t
    return profile

def getConsensus(profile):
    consensus = ""
    nuc = ["A", "C", "G", "T"]
    for i in range(len(profile[0])):
        positions = []
        positions.append(profile[0][i])
        positions.append(profile[1][i])
        positions.append(profile[2][i])
        positions.append(profile[3][i])
        consensus += nuc[positions.index(max(positions))]
    return consensus

def scoreMotifsOnConsensus(consensus, motifs):
    score = 0
    for string in motifs:
        score += hammingDistance(consensus, string)
    return score

def greedyMotifSearch(dna, t, k):
    bestMotifs = []
    for string in dna:
        bestMotifs.append(string[0:k])
    for i in range((len(dna[0]) - k) + 1):
        currentKmer = dna[0][i:i+k]
        motifs = [currentKmer]
        profile = generateProfile(motifs)
        for j in range(1, t):
            motifs.append(profileMostProbable(dna[j], k, profile))
            profile = generateProfile(motifs)
        consensus = getConsensus(profile)
        if scoreMotifsOnConsensus(consensus, motifs) < scoreMotifsOnConsensus(getConsensus(generateProfile(bestMotifs)), bestMotifs):
            bestMotifs = motifs
    return bestMotifs

dna = []
out = []
with open("data/week_three_gms.txt") as data:
    data = data.readlines()
    for i in range(len(data)):
        data[i] = data[i].strip("\n")
    for lines in data[1:]:
        dna.append(lines)
    out = greedyMotifSearch(dna, 25, 12)

with open("data/gms_output.txt", "w") as output:
    for lines in out:
        output.write(lines)
        output.write("\n")
