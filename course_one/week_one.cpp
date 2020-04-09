#include <iostream>
#include <map>
#include <string>

/*
All of the python algorithms rewritten for c++
Forgive the shoddy project structure...I'm learning
*/

// A fucntion that takes some part of a genome and a k-mer pattern and returns the no. of occurences
int kmer_frequency(std::string genome, std::string kmer);

// A fucntion to determine the frequency of some kmer lengths. how mnay 5-mers in this part of the genome...?
std::map<std::string, int> frequent_kmer(std::string genome, int kmer);

// Returns the reverse compliment of a kmer
std::string return_compliment(std::string kmer);

// Returns an array 

int kmer_frequency(std::string genome, std::string kmer) {
    int frequency = 0;
    int iter_len = ((genome.length() - kmer.length()) + 1);
    for (int i = 0; i < iter_len; i++) {
        // compares kmer to the current part of genome we are examining
        if (kmer.compare(genome.substr(i, kmer.length())) == 0) {
            frequency++;
        }
    }
    return frequency; 
}

std::map<std::string, int> frequent_kmer(std::string genome, int kmer) {
    std::map<std::string, int> frequency_map;
    int iter_len = ((genome.length() - kmer) + 1);
    std::string current_kmer;
    for (int i = 0; i < iter_len; i++) {
        current_kmer = genome.substr(i, kmer);
        if (frequency_map.count(current_kmer) == 1) {
            frequency_map[current_kmer] += 1;
        }
        else {
            frequency_map[current_kmer] = 1;
        }    
    }
    return frequency_map;
}

/*
    This is a poor solution for long genome's. It can take up to linear time to insert a char
    into a string. A better solution would be to store the data as an array, index it backwards,
    adding each char to a new array fromm start to end. This array could then be turned into
    a string at another time. This would be much closer to linear time than the ~ O(n^2) 
    approach below
*/
std::string return_compliment(std::string kmer) {
    std::string reverse_comp = "";
    char curr_nuc;
    for (int i = 0; i < kmer.length(); i++) {
        curr_nuc = kmer[i];
        switch (curr_nuc) {
            case 'A':
                reverse_comp.insert(0, "T");
                break;
            case 'T':
                reverse_comp.insert(0, "A");
                break;
            case 'G':
                reverse_comp.insert(0, "C");
                break;
            case 'C':
                reverse_comp.insert(0, "G");
                break;
        }
    }
    return reverse_comp;
}

int main(void) {

    /*
    std::string pat = "CAAGGAAGTTGCCGTTGCCTGCCCTATCAATTACGCGACGTAGGTCTCACCAACGGTTATGTTAGGCTA";
    int kmer = 1;
    std::map<std::string, int> freq = frequent_kmer(pat, kmer);
    std::map<std::string, int>::iterator it = freq.begin();
    for (std::pair<std::string, int> element : freq) {
		// Accessing KEY from element
		std::string word = element.first;
		// Accessing VALUE from element.
		int count = element.second;
		std::cout << word << " :: " << count << std::endl;
	}
    */
   std::string reved = return_compliment("ATCG");
   std::cout << reved;
   return 0;
}
