// CONTACT: yunheo1@illinois.edu

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <unordered_map>
#include <sstream>
#include <iterator>

int main (int argc, char** argv) {
   // check the number of arguments
   if (argc != 4) {
      std::cout << std::endl << "USAGE: " << argv[0] << " <original fastq file> <fastq file to be mapped> <output map file>" << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // open the input fastq file
   std::ifstream f_in_mapped;
   f_in_mapped.open(argv[2]);

   if (f_in_mapped.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << argv[2] << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   //--------------------------------------------------
   // construct a hash table
   //--------------------------------------------------
   std::unordered_map<std::string, short int> map_read_occurrence;

   // iterate reads
   std::string line_header;

   bool already_warned(false);

   getline(f_in_mapped, line_header);

   while (!f_in_mapped.eof()) {
      // count the number of words in the header
      if (already_warned == false) {
         std::size_t num_words(std::distance(std::istream_iterator<std::string>(std::istringstream(line_header) >> std::ws), std::istream_iterator<std::string>()));
         // multiple words in the header
         if (num_words > 1) {
            std::cout << std::endl << "WARNING: Multiple words in the header line_header. Only the 1st word will be used.\n\n";

            // do not warn it any more
            already_warned = true;
         }
      }

      // get the 1st word from the header
      std::istringstream iss_header(line_header);
      std::string word_1st;
      iss_header >> word_1st;

      // remove "@"
      word_1st.erase(0, 1);

      // word_1st already exists in the hash table
      if (map_read_occurrence.find(word_1st) != map_read_occurrence.end()) {
         map_read_occurrence[word_1st]++;
      }
      // word_1st does not exist in the hash table
      // add it to the table
      else {
         map_read_occurrence[word_1st] = 1;
      }

      // read remaining lines of the read
      getline(f_in_mapped, line_header);
      getline(f_in_mapped, line_header);
      getline(f_in_mapped, line_header);

      // new header
      getline(f_in_mapped, line_header);
   }

   f_in_mapped.close();

   //--------------------------------------------------
   // write a map file
   //--------------------------------------------------
   // open the original read file
   std::ifstream f_in_original;
   f_in_original.open(argv[1]);

   if (f_in_original.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << argv[1] << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // open the output file
   std::ofstream f_out;
   f_out.open(argv[3]);

   if (f_out.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << argv[3] << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   already_warned = false;

   getline(f_in_original, line_header);

   while (!f_in_original.eof()) {
      // count the number of words in the header
      if (already_warned == false) {
         std::size_t num_words(std::distance(std::istream_iterator<std::string>(std::istringstream(line_header) >> std::ws), std::istream_iterator<std::string>()));
         // multiple words in the header
         if (num_words > 1) {
            std::cout << std::endl << "WARNING: Multiple words in the header line_header. Only the 1st word will be used.\n\n";

            // do not warn it any more
            already_warned = true;
         }
      }

      // get the 1st word from the header
      std::istringstream iss_header(line_header);
      std::string word_1st;
      iss_header >> word_1st;

      // remove "@"
      word_1st.erase(0, 1);

      // word_1st does not exist in the hash table
      if (map_read_occurrence.find(word_1st) == map_read_occurrence.end()) {
         f_out << word_1st << " 0" << std::endl;
      }
      // word_1st already exists in the hash table
      else {
         f_out << word_1st << " " << map_read_occurrence[word_1st] << std::endl;
      }

      // skip remaining lines
      getline(f_in_original, line_header);
      getline(f_in_original, line_header);
      getline(f_in_original, line_header);
      getline(f_in_original, line_header);
   }

   f_in_original.close();
   f_out.close();
}
