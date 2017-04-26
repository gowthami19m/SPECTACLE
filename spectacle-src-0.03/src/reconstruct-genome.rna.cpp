// CONTACT: yunheo1@illinois.edu

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

int main (int argc, char** argv) {
   // check the number of arguments
   if (argc != 3) {
      std::cout << std::endl << "USAGE: " << argv[0] << " <information file> <output file>" << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // open the input file
   std::ifstream f_in;
   f_in.open(argv[1]);

   if (f_in.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << argv[1] << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // read the last index
   // 1-based
   std::size_t last_index;
   f_in >> last_index;

   if (f_in.good() == false) {
      std::cout << std::endl << "ERROR: Cannot read the last index" << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // initialize the genome sequence
   std::string genome_sequence(last_index, 'N');

   // read each exon
   std::size_t start_index;
   std::size_t length;
   std::string exon_seq;
   f_in >> start_index >> length >> exon_seq;

   while (f_in.good() == true) {
      genome_sequence.replace(start_index, length, exon_seq);

      // read the next exon
      f_in >> start_index >> length >> exon_seq;
   }

   f_in.close();

   // open the output file
   std::ofstream f_out;
   f_out.open(argv[2]);

   if (f_out.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << argv[2] << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   f_out << genome_sequence << std::endl;

   f_out.close();
}
