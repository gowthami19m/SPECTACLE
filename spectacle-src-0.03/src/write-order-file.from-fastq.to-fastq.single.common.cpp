// CONTACT: yunheo1@illinois.edu

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <regex>
#include <unordered_map>
#include <sstream>
#include <iterator>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

int main (int argc, char** argv) {
   // check the number of arguments
   if (argc != 4) {
      std::cout << std::endl << "USAGE: " << argv[0] << " <original fastq file> <fastq file to be reordered> <output order file>" << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // check whether the input is gzipped
   std::string in_fastq(argv[1]);

   bool is_gzip(false);

   std::smatch smatch1;

   std::regex rx_gz_ext("\\.gz$");

   if (std::regex_search(in_fastq, smatch1, rx_gz_ext)) {
      is_gzip = true;
   }   

   // open the original fastq file
   std::ifstream f_in_original;

   if (is_gzip) {
      f_in_original.open(argv[1], std::ios_base::binary);
   }   
   else {
      f_in_original.open(argv[1]);
   }   

   if (f_in_original.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << argv[1] << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // set the io filter
   boost::iostreams::filtering_istream f_in_original_filter;
   if (is_gzip) {
      f_in_original_filter.push(boost::iostreams::gzip_decompressor());
      f_in_original_filter.push(f_in_original);
   }   
   else {
      f_in_original_filter.push(f_in_original);
   }   

   //--------------------------------------------------
   // construct a hash table
   //--------------------------------------------------
   std::unordered_map<std::string, std::size_t> map_read_order;

   // iterate reads
   std::string line_header;

   std::size_t num_reads(0);

   bool already_warned(false);

   getline(f_in_original_filter, line_header);

   while (!f_in_original_filter.eof()) {
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

      // increment the number of reads
      num_reads++;

      // get the 1st word from the header
      std::istringstream iss_header(line_header);
      std::string word_1st;
      iss_header >> word_1st;

      // remove "@"
      word_1st.erase(0, 1);

      // word_1st already exists in the hash table
      if (map_read_order.find(word_1st) != map_read_order.end()) {
         std::cout << std::endl << "ERROR: " << word_1st << " exists multiple times in " << argv[1] << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
      // word_1st does not exist in the hash table
      // add it to the table
      else {
         map_read_order[word_1st] = num_reads;
      }

      // read remaining lines of the read
      getline(f_in_original_filter, line_header);
      getline(f_in_original_filter, line_header);
      getline(f_in_original_filter, line_header);

      // new header
      getline(f_in_original_filter, line_header);
   }

   f_in_original.close();

   //--------------------------------------------------
   // reorder a new fastq file
   //--------------------------------------------------
   // open the aligned fastq file
   std::ifstream f_in_aligned;
   f_in_aligned.open(argv[2]);

   if (f_in_aligned.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << argv[2] << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // open the output file
   std::ofstream f_out;
   f_out.open(argv[3]);

   if (f_out.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << argv[3] << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   getline(f_in_aligned, line_header);

   while (!f_in_aligned.eof()) {
      // get the 1st word from the header
      std::istringstream iss_header(line_header);
      std::string word_1st;
      iss_header >> word_1st;

      // remove "@"
      word_1st.erase(0, 1);

      // word_1st does not exist in the hash table
      if (map_read_order.find(word_1st) == map_read_order.end()) {
         std::cout << std::endl << "ERROR: " << word_1st << " does not exist in " << argv[1] << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
      // word_1st already exists in the hash table
      else {
         f_out << map_read_order[word_1st] << std::endl;
      }

      // read remaining lines of the read
      getline(f_in_aligned, line_header);
      getline(f_in_aligned, line_header);
      getline(f_in_aligned, line_header);

      // new header
      getline(f_in_aligned, line_header);
   }

   f_in_aligned.close();
   f_out.close();
}
