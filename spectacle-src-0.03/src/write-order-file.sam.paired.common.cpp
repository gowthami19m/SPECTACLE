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
      std::cout << std::endl << "USAGE: " << argv[0] << " <error location file> <sam file> <output order file>" << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // open the error location file
   std::ifstream f_location;
   f_location.open(argv[1]);

   if (f_location.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << argv[1] << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   //--------------------------------------------------
   // construct a hash table
   //--------------------------------------------------
   std::unordered_map<std::string, std::size_t> map_read_order;

   // iterate reads
   std::string line_location;

   std::size_t num_reads(0);

   getline(f_location, line_location);

   while (!f_location.eof()) {
      // increment the number of reads
      num_reads++;

      // get the 1st word from the header
      std::istringstream iss_header(line_location);
      std::string read_name;
      iss_header >> read_name;

      // remove "/1"
      std::size_t read_name_length(read_name.length());
      if (read_name_length >= 3) {
         if ((read_name[read_name_length - 2] == '/') && (read_name[read_name_length - 1] == '1')) {
            read_name.erase(read_name_length - 2, 2);
         }
      }

      // read_name already exists in the hash table
      if (map_read_order.find(read_name) != map_read_order.end()) {
         std::cout << std::endl << "ERROR: " << read_name << " exists multiple times in " << argv[1] << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
      // read_name does not exist in the hash table
      // add it to the table
      else {
         map_read_order[read_name] = num_reads;
      }

      // skip the reverse read
      getline(f_location, line_location);

      // new forward read
      getline(f_location, line_location);
   }

   f_location.close();

   // open the sam file
   std::ifstream f_sam;
   f_sam.open(argv[2]);

   if (f_sam.is_open() == false) {
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

   //--------------------------------------------------
   // iterate reads in the sam file
   //--------------------------------------------------
   getline(f_sam, line_location);

   while (!f_sam.eof()) {
      // non-header lines
      if (line_location[0] != '@') {
         std::istringstream iss_header;

         // get a read name and a flag
         std::string read_name_1st;
         unsigned int flag_1st;

         iss_header.str(line_location);
         iss_header >> read_name_1st;
         iss_header >> flag_1st;

         // read its pair
         getline(f_sam, line_location);

         // get a read name and a flag
         std::string read_name_2nd;
         unsigned int flag_2nd;

         iss_header.str(line_location);
         iss_header >> read_name_2nd;
         iss_header >> flag_2nd;

         // determine a 1st read name
         //std::string read_name;

         // remove "/1" and "/2"
         std::size_t read_name_length;

         read_name_length = read_name_1st.length();
         if (read_name_length >= 3) {
            if ((read_name_1st[read_name_length - 2] == '/') && (read_name_1st[read_name_length - 1] == '1')) {
               read_name_1st.erase(read_name_length - 2, 2);
            }
         }

         read_name_length = read_name_2nd.length();
         if (read_name_length >= 3) {
            if ((read_name_2nd[read_name_length - 2] == '/') && (read_name_2nd[read_name_length - 1] == '1')) {
               read_name_2nd.erase(read_name_length - 2, 2);
            }
         }

         // both names are same
         if (read_name_1st == read_name_2nd) {
            // read_name does not exist in the hash table
            if (map_read_order.find(read_name_1st) == map_read_order.end()) {
               std::cout << std::endl << "ERROR: " << read_name_1st << " does not exist in " << argv[1] << std::endl << std::endl;
               exit(EXIT_FAILURE);
            }
            // read_name already exists in the hash table
            else {
               f_out << map_read_order[read_name_1st] << std::endl;
            }
         }
         // both names are different
         // find the 1st read using flags
         else {
            // read_name_1st does not exist in the hash table
            if (map_read_order.find(read_name_1st) == map_read_order.end()) {
               // read_name_2nd does not exist in the hash table
               if (map_read_order.find(read_name_2nd) == map_read_order.end()) {
                  std::cout << std::endl << "ERROR: Both " << read_name_1st << " and " << read_name_2nd << " do not exist in " << argv[1] << std::endl << std::endl;
                  exit(EXIT_FAILURE);
               }
               // read_name_2nd exists in the hash table
               else {
                  f_out << map_read_order[read_name_2nd] << std::endl;
               }

            }
            // read_name_1st already exists in the hash table
            else {
               f_out << map_read_order[read_name_1st] << std::endl;
            }
         }

         // new header
         getline(f_sam, line_location);
      }
      // filter out header lines
      else {
         getline(f_sam, line_location);
      }
   }

   f_sam.close();
   f_out.close();
}
