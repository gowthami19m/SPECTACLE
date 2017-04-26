//----------------------------------------------------------------------
// libraries
//----------------------------------------------------------------------
//
// c libraries
//
#include <cstdlib>
#include <cstring>

//
// c++ libraries
//
#include <iomanip>
#include <iostream>
#include <regex>
#include <string>
#include <unordered_map>
#include <vector>

// DEBUG
#include <unistd.h>
#include <sys/resource.h>



//----------------------------------------------------------------------
// definitins
//----------------------------------------------------------------------
#define SMALL_NUMBER -1000000
#define MIN_OVERLAP  30



//----------------------------------------------------------------------
// variables
//----------------------------------------------------------------------
int matrix_width;
int matrix_height;
int match_gain;
int mismatch_penalty;
int gap_opening_penalty;
int gap_extension_penalty;
int string1_length;
int string2_length;
int outer_length_5_end;
int outer_length_3_end;
int start_index;
int end_index;
int corrected_read_length;
int read_length;
int num_yyns_substitution_local_best;
int num_ynys_substitution_local_best;
int num_nyys_substitution_local_best;
int num_nyns_substitution_local_best;
int num_nnns_substitution_local_best;
int num_yyns_insertion_local_best;
int num_nyys_insertion_local_best;
int num_nyns_insertion_local_best;
int num_nnns_insertion_local_best;
int num_yyns_deletion_local_best;
int num_nyys_deletion_local_best;
int num_nyns_deletion_local_best;
int num_nnns_deletion_local_best;
int num_from_substitution_to_deletion_local_best;
int num_nyys_substitution_trim_local_best;
int num_nyys_insertion_trim_local_best;
int num_nyys_deletion_trim_local_best;
int num_deletions_5_prime_best;
int num_deletions_3_prime_best;
int alignment_score_best;
int alignment_score_new_best;
int num_substitutions;
int num_insertions;
int num_insertions_unit;
int num_deletions;
int matrix_size;
int longest_alignment_length;
int highest_score;
int max_read_length;
int num_not_evaluated_substitution;
int num_not_evaluated_insertion;
int num_not_evaluated_deletion;
int num_total_bases_percent_similarity;
int num_matched_bases_percent_similarity;
int ref_seq_index;

unsigned int max_candidates;

std::string string1;
std::string string2;
std::string outer_5_end;
std::string outer_3_end;
std::string read_name;
std::string alignment_best;
std::string error_index_best;
std::string substitutions;
std::string insertions;
std::string deletions;
std::string random_alignment1;
std::string random_alignment2;
std::string strand;

int* match_matrix;
int* gap_1_matrix;
int* gap_2_matrix;

std::vector<int> position_vector_local_best_tmp;
std::vector<int> corrected_position_vector_local_best_tmp;

std::vector<std::string> alignment1_vector;
std::vector<std::string> alignment2_vector;

std::unordered_map<int, char>        substitution_org_map;
std::unordered_map<int, char>        substitution_err_map;
std::unordered_map<int, std::string> insertion_map;
std::unordered_map<int, char>        deletion_map;

bool no_end_gap_penalty;
bool is_trimmed;
bool is_detail;
bool too_many_candidates;



//----------------------------------------------------------------------
// functions
//----------------------------------------------------------------------
// DEBUG
//
// get_current_rss
//
std::size_t get_current_rss() { 
   long rss = 0L; 

   FILE* fp = NULL;

   if ((fp = fopen("/proc/self/statm", "r")) == NULL) {
      return (std::size_t)0L;
   }

   if ( fscanf(fp, "%*s%ld", &rss ) != 1) {   
      fclose(fp);
      return (std::size_t)0L;
   }   
   fclose(fp);
   return (std::size_t)rss * (std::size_t)sysconf( _SC_PAGESIZE);
}

//
// max3
//
inline int max3(int a, int b, int c) {
   int result(a);

   if (b > result) {
      result = b;
   }

   if (c > result) {
      result = c;
   }

   return result;
}

//
// max4
//
inline int max4(int a, int b, int c, int d) {
   int result(a);

   if (b > result) {
      result = b;
   }

   if (c > result) {
      result = c;
   }

   if (d > result) {
      result = d;
   }

   return result;
}

//
// complement_base
//
inline char complement_base(char in_char) {
   if (in_char == 'A') {
      return 'T';
   }
   else if (in_char == 'C') {
      return 'G';
   }
   else if (in_char == 'G') {
      return 'C';
   }
   else if (in_char == 'T') {
      return 'A';
   }
   else {
      std::cout << "\nERROR: Wrong character " << in_char << "\n\n";
      exit(EXIT_FAILURE);
   }
}

//
// initialize_variables
//
void initialize_variables() {
   alignment_score_best     = SMALL_NUMBER;
   alignment_score_new_best = SMALL_NUMBER;
   num_substitutions        = 0;
   num_insertions           = 0;
   num_insertions_unit      = 0;
   num_deletions            = 0;
   too_many_candidates      = false;
   matrix_width             = string1.length() + 1;
   matrix_height            = string2.length() + 1;
   string1_length           = string1.length();
   string2_length           = string2.length();
   matrix_size              = matrix_width * matrix_height;
   longest_alignment_length = string1_length + string2_length;

   alignment1_vector.clear();
   alignment2_vector.clear();

   // the clear function does not release memory if the memory size is smaller than a threshold
   alignment1_vector.shrink_to_fit();
   alignment2_vector.shrink_to_fit();
}

//
// decode_errors
//
void decode_errors() {
   //
   // substitution
   //
   substitution_org_map.clear();
   substitution_err_map.clear();

   std::regex rx_substitution("([0-9]+):([ACGT])->([ACGT]);");

   if (substitutions != "-") {
      for(auto it_match = std::sregex_iterator(substitutions.begin(), substitutions.end(), rx_substitution); it_match != std::sregex_iterator(); it_match++) {
         substitution_org_map[atoi((*it_match)[1].str().c_str())] = (*it_match)[2].str()[0];
         substitution_err_map[atoi((*it_match)[1].str().c_str())] = (*it_match)[3].str()[0];
         num_substitutions++;
      }
   }

   //
   // insertion
   //
   insertion_map.clear();

   std::regex rx_insertion("([0-9]+):([ACGT]+);");

   if (insertions != "-") {
      for(auto it_match = std::sregex_iterator(insertions.begin(), insertions.end(), rx_insertion); it_match != std::sregex_iterator(); it_match++) {
         insertion_map[atoi((*it_match)[1].str().c_str())] = (*it_match)[2].str();

         num_insertions_unit += (*it_match)[2].str().length();
         num_insertions++;
      }
   }

   //
   // deletion
   //
   deletion_map.clear();

   std::regex rx_deletion("([0-9]+):([ACGT]);");

   if (deletions != "-") {
      for(auto it_match = std::sregex_iterator(deletions.begin(), deletions.end(), rx_deletion); it_match != std::sregex_iterator(); it_match++) {
         deletion_map[atoi((*it_match)[1].str().c_str())] = (*it_match)[2].str()[0];
         num_deletions++;
      }
   }
}

//
// fill_matrixes
//
void fill_matrixes() {
   // resize matrixes
   match_matrix = new int[matrix_size];
   gap_1_matrix = new int[matrix_size];
   gap_2_matrix = new int[matrix_size];

   //----------------------------------------------------------------------
   // initialize matrixes
   //----------------------------------------------------------------------
   //
   // match
   //
   // (0, 0)
   match_matrix[0] = 0;

   // first row
   for(int it_x = 1; it_x < matrix_width; it_x++) {
      match_matrix[it_x] = SMALL_NUMBER;
   }

   // first column
   for(int it_y = 1; it_y < matrix_height; it_y++) {
      match_matrix[matrix_width * it_y] = SMALL_NUMBER;
   }

   //
   // gap 1
   //
   // first row (including origin)
   for(int it_x = 0; it_x < matrix_width; it_x++) {
      gap_1_matrix[it_x] = SMALL_NUMBER;
   }

   // first column (excluding origin)
   for(int it_y = 1; it_y < matrix_height; it_y++) {
      if (no_end_gap_penalty) {
         gap_1_matrix[matrix_width * it_y] = 0;
      }
      else {
         gap_1_matrix[matrix_width * it_y] = gap_opening_penalty + it_y * gap_extension_penalty;
      }
   }

   //
   // gap 2
   //
   // first row (excluding origin)
   for(int it_x = 1; it_x < matrix_width; it_x++) {
      if (no_end_gap_penalty) {
         gap_2_matrix[it_x] = 0;
      }
      else {
         gap_2_matrix[it_x] = gap_opening_penalty + it_x * gap_extension_penalty;
      }
   }

   // first column (including origin)
   for(int it_y = 0; it_y < matrix_height; it_y++) {
      gap_2_matrix[matrix_width * it_y] = SMALL_NUMBER;
   }

   //----------------------------------------------------------------------
   // fill remainings
   //----------------------------------------------------------------------
   // start from (1, 1)
   int matrix_index(matrix_width + 1);
   int index_up(1);
   int index_left(matrix_width);
   int index_up_left(0);

   int match_penalty_tmp;

   for (int index_y = 0; index_y < string2_length; index_y++) {
      for (int index_x = 0; index_x < string1_length; index_x++) {
         //
         // match
         //
         if (string1[index_x] == string2[index_y]) {
            match_penalty_tmp = match_gain;
         }
         else {
            match_penalty_tmp = mismatch_penalty;
         }

         match_matrix[matrix_index] = max3(
                                            match_matrix[index_up_left] + match_penalty_tmp,
                                            gap_1_matrix[index_up_left] + match_penalty_tmp,
                                            gap_2_matrix[index_up_left] + match_penalty_tmp
                                           );

         //
         // gap 1
         //
         // the last column in the no end gap mode
         if (no_end_gap_penalty && (index_x == (string1_length - 1))) {
            gap_1_matrix[matrix_index] = max3(
                                               match_matrix[index_up],
                                               gap_1_matrix[index_up],
                                               gap_2_matrix[index_up]
                                              );
         }
         // other columns
         else {
            gap_1_matrix[matrix_index] = max3(
                                               match_matrix[index_up] + gap_opening_penalty + gap_extension_penalty,
                                               gap_1_matrix[index_up] + gap_extension_penalty,
                                               gap_2_matrix[index_up] + gap_opening_penalty + gap_extension_penalty
                                              );
         }

         //
         // gap 2
         //
         // the last row in the no end gap mode
         if (no_end_gap_penalty && (index_y == (string2_length - 1))) {
            gap_2_matrix[matrix_index] = max3(
                                               match_matrix[index_left],
                                               gap_1_matrix[index_left],
                                               gap_2_matrix[index_left]
                                              );
         }
         // normal situation
         else {
            gap_2_matrix[matrix_index] = max3(
                                               match_matrix[index_left] + gap_opening_penalty + gap_extension_penalty,
                                               gap_1_matrix[index_left] + gap_opening_penalty + gap_extension_penalty,
                                               gap_2_matrix[index_left] + gap_extension_penalty
                                              );
         }

         // skep the first column
         matrix_index++;
         index_up++;
         index_left++;
         index_up_left++;
      }

      matrix_index++;
      index_up++;
      index_left++;
      index_up_left++;
   }
}

//
// print_matrixes
//
void print_matrixes() {
   //
   // match
   //
   std::cout << "\n";
   std::cout << "Match:\n";

   for (int it_matrix = 0; it_matrix < matrix_size; it_matrix++) {
      std::cout << " " << std::setw(12) << match_matrix[it_matrix];

      if ((it_matrix % matrix_width) == (matrix_width - 1)) {
         std::cout << "\n";
      }
   }

   std::cout << "\n";

   //
   // gap 1
   //
   std::cout << "Gap 1:\n";

   for (int it_matrix = 0; it_matrix < matrix_size; it_matrix++) {
      std::cout << " " << std::setw(12) << gap_1_matrix[it_matrix];

      if ((it_matrix % matrix_width) == (matrix_width - 1)) {
         std::cout << "\n";
      }
   }

   std::cout << "\n";

   //
   // gap 2
   //
   std::cout << "Gap 2:\n";

   for (int it_matrix = 0; it_matrix < matrix_size; it_matrix++) {
      std::cout << " " << std::setw(12) << gap_2_matrix[it_matrix];

      if ((it_matrix % matrix_width) == (matrix_width - 1)) {
         std::cout << "\n";
      }
   }

   std::cout << "\n";
}

//
// traceback
//
inline void traceback(std::string alignment1, std::string alignment2, int index_x, int index_y, char current_matrix, int alignment_index, int matrix_index, int current_score) {
   // first row or first column
   if ((index_x == 0) || (index_y == 0)) {
      // gap in string1
      while(index_y > 0) {
         alignment1[alignment_index] = '-';
         alignment2[alignment_index] = string2[index_y - 1];

         alignment_index--;
         index_y--;
      }

      // gap in string2
      while(index_x > 0) {
         alignment1[alignment_index] = string1[index_x - 1];
         alignment2[alignment_index] = '-';

         alignment_index--;
         index_x--;
      }

      // remove unnecessary space in aliment1/2
      int start_index_tmp(alignment_index + 1);
      int real_alignment_len(longest_alignment_length - start_index_tmp);

      // too few overlaps between two sequences
      // length of the overlapped region = start_index_tmp
      // example: 0 overlap
      // alignment1: CAGGAGGCGGAGGTTGCAGTGAGCCGAGATCATGCCACTGCACTCCAGCCTGGGCAACAAGAGCAAAACTCAGTCTCAAAAAAAAAAAAAAAAAAAAAAA--------------------------------------------------
      // alignment2: ----------------------------------------------------------------------------------------------------CAGGTGACAGAGCCAGACTCCATCTCAAAAACAAACAAACAAACAAAAAA
      if (start_index_tmp < MIN_OVERLAP) {
         too_many_candidates = true;
      }
      else {
         alignment1 = alignment1.substr(start_index_tmp, real_alignment_len);
         alignment2 = alignment2.substr(start_index_tmp, real_alignment_len);

         alignment1_vector.push_back(alignment1);
         alignment2_vector.push_back(alignment2);

         if (alignment1_vector.size() > max_candidates) {
            too_many_candidates = true;
         }
      }
   }
   // not the first row or the first column
   else {
      // update alignment1/2
      switch(current_matrix) {
         case 'M':
            alignment1[alignment_index] = string1[index_x - 1];
            alignment2[alignment_index] = string2[index_y - 1];
            break;

         case '1':
            alignment1[alignment_index] = '-';
            alignment2[alignment_index] = string2[index_y - 1];
            break;

         case '2':
            alignment1[alignment_index] = string1[index_x - 1];
            alignment2[alignment_index] = '-';
            break;

         default:
            std::cout << "\nERROR: Wrong matrix specifier" << current_matrix << "\n\n";
            exit(EXIT_FAILURE);
      }

      // set penalties
      int match_penalty_tmp;
      int gap_1_opening_penalty;
      int gap_2_opening_penalty;
      int gap_1_extension_penalty;
      int gap_2_extension_penalty;

      int prev_match_penalty_tmp;
      int prev_gap_1_penalty;
      int prev_gap_2_penalty;

      if (string1[index_x - 1] == string2[index_y - 1]) {
         match_penalty_tmp = match_gain;
      }
      else {
         match_penalty_tmp = mismatch_penalty;
      }

      gap_1_opening_penalty = gap_opening_penalty + gap_extension_penalty;
      gap_2_opening_penalty = gap_opening_penalty + gap_extension_penalty;

      gap_1_extension_penalty = gap_extension_penalty;
      gap_2_extension_penalty = gap_extension_penalty;

      // no end gap
      if (no_end_gap_penalty) {
         // first or last column
         if ((index_x == 0) || (index_x == string1_length)) {
            gap_1_opening_penalty   = 0;
            gap_1_extension_penalty = 0;
         }

         // first or last row
         if ((index_y == 0) || (index_y == string2_length)) {
            gap_2_opening_penalty   = 0;
            gap_2_extension_penalty = 0;
         }
      }

      // update indices
      switch(current_matrix) {
         case 'M':
            prev_match_penalty_tmp = match_penalty_tmp;
            prev_gap_1_penalty     = match_penalty_tmp;
            prev_gap_2_penalty     = match_penalty_tmp;
            index_x--;
            index_y--;
            matrix_index -= (matrix_width + 1);
            break;

         case '1':
            prev_match_penalty_tmp = gap_1_opening_penalty;
            prev_gap_1_penalty     = gap_1_extension_penalty;
            prev_gap_2_penalty     = gap_1_opening_penalty;
            index_y--;
            matrix_index -= matrix_width;
            break;

          case '2':
            prev_match_penalty_tmp = gap_2_opening_penalty;
            prev_gap_1_penalty     = gap_2_opening_penalty;
            prev_gap_2_penalty     = gap_2_extension_penalty;
            index_x--;
            matrix_index--;
            break;

          default:
            std::cout << "\nERROR: Wrong matrix specifier" << current_matrix << "\n\n";
            exit(EXIT_FAILURE);
      }

      alignment_index--;

      //----------------------------------------------------------------------
      // trace recursively
      //----------------------------------------------------------------------
      int current_score_tmp;

      //
      // match
      //
      if (((match_matrix[matrix_index] + prev_match_penalty_tmp) == current_score) && (too_many_candidates == false)) {
         current_matrix    = 'M';
         current_score_tmp = match_matrix[matrix_index];

         traceback(alignment1, alignment2, index_x, index_y, current_matrix, alignment_index, matrix_index, current_score_tmp);
      }

      //
      // gap in string1
      //
      if ((gap_1_matrix[matrix_index] + prev_gap_1_penalty == current_score) && (too_many_candidates == false)) {
         current_matrix    = '1';
         current_score_tmp = gap_1_matrix[matrix_index];

         traceback(alignment1, alignment2, index_x, index_y, current_matrix, alignment_index, matrix_index, current_score_tmp);
      }

      //
      // gap in string2
      //
      if ((gap_2_matrix[matrix_index] + prev_gap_2_penalty == current_score) && (too_many_candidates == false)) {
         current_matrix    = '2';
         current_score_tmp = gap_2_matrix[matrix_index];

         traceback(alignment1, alignment2, index_x, index_y, current_matrix, alignment_index, matrix_index, current_score_tmp);
      }
   }
}

//
// evaluate_each_alignment
//
inline void evaluate_each_alignment(std::string alignment1, std::string alignment2, int* position_vector_local_best, int* corrected_position_vector_local_best) {
   //----------------------------------------------------------------------
   // variables
   //----------------------------------------------------------------------
   int num_yyns_substitution_tmp(0); 
   int num_ynys_substitution_tmp(0); 
   int num_nyys_substitution_tmp(0); 
   int num_nyns_substitution_tmp(0); 
   int num_nnns_substitution_tmp(0); 

   int num_yyns_insertion_tmp(0); 
   int num_nyys_insertion_tmp(0); 
   int num_nyns_insertion_tmp(0); 
   int num_nnns_insertion_tmp(0); 

   int num_yyns_deletion_tmp(0); 
   int num_nyys_deletion_tmp(0); 
   int num_nyns_deletion_tmp(0); 
   int num_nnns_deletion_tmp(0); 

   int num_from_substitution_to_deletion_tmp(0); 

   int num_nyys_substitution_trim_tmp(0); 

   int num_nyys_insertion_trim_tmp(0); 

   int num_nyys_deletion_trim_tmp(0); 

   int num_deletions_5_prime_tmp;
   int num_deletions_3_prime_tmp;

   int alignment_score(0);
   int alignment_score_new(0);

   int how_many_insertions;

   std::vector<int> position_vector_tmp(position_vector_local_best, position_vector_local_best + max_read_length);
   std::vector<int> corrected_position_vector_tmp(corrected_position_vector_local_best, corrected_position_vector_local_best + max_read_length);

   std::smatch smatch1;

   std::string error_index_tmp("");

   char base_ref_seq;
   char base_org_read;
   char base_mod_read;

   bool flag_matched_prev;

   std::string alignment1_keep(alignment1);
   std::string alignment2_keep(alignment2);

   //----------------------------------------------------------------------
   // insertions at the 5'-end
   //----------------------------------------------------------------------
   int num_insertions_5_prime(0);

   std::regex rx_start_gap("^(-+)[ACGT]+-*");

   // count the number of insertions at the 5'-end
   if (std::regex_match(alignment1, smatch1, rx_start_gap)) {
      num_insertions_5_prime = smatch1[1].length();
   }

   how_many_insertions = num_insertions_5_prime;
   flag_matched_prev   = true;

   for (int it_5_end = 0; it_5_end < num_insertions_5_prime; it_5_end++) {
      // REF: base in the reference
      // not sufficient outer length
      if ((outer_length_5_end - num_insertions_5_prime + it_5_end) < 0) {
         base_ref_seq = '-';
      }
      // sufficient outer length
      else {
         base_ref_seq = outer_5_end[outer_length_5_end - num_insertions_5_prime + it_5_end];
      }

      // MOD: base in the modified read
      base_mod_read = alignment2[it_5_end];

      // REF: -
      if (base_ref_seq == '-') {
         num_yyns_insertion_tmp++;

         if (flag_matched_prev) {
            alignment_score     += (gap_opening_penalty + gap_extension_penalty);
            alignment_score_new += (gap_opening_penalty + gap_extension_penalty);
         }
         else {
            alignment_score     += gap_extension_penalty;
            alignment_score_new += gap_extension_penalty;
         }

         flag_matched_prev = false;
      }
      // REF: A/C/G/T
      else {
         // REF: A/C/G/T
         // MOD: A/C/G/T
         if (base_mod_read != '-') {
            // REF: A/C/G/T
            // MOD: = REF
            if (base_mod_read == base_ref_seq) {
               num_ynys_substitution_tmp++;

               flag_matched_prev = true;
            }
            // REF: A/C/G/T
            // MOD: != REF
            else {
               num_yyns_insertion_tmp++;

               if (flag_matched_prev) {
                  alignment_score     += (gap_opening_penalty + gap_extension_penalty);
                  alignment_score_new += (gap_opening_penalty + gap_extension_penalty);
               }
               else {
                  alignment_score     += gap_extension_penalty;
                  alignment_score_new += gap_extension_penalty;
               }

               flag_matched_prev = false;
            }
         }
         // REF: A/C/G/T
         // MOD: -
         else {
            std::cout << "\nERROR: Check this " << alignment1 << alignment2 << "\n\n";
         }
      }
   }

   // remove the first num_insertion_5_prime bases and the last num_insertion_3_prime bases ("-"s in alignment1)
   // from alignment1 and alignment2
   alignment1 = alignment1.substr(num_insertions_5_prime, alignment1.length() - num_insertions_5_prime);
   alignment2 = alignment2.substr(num_insertions_5_prime, alignment2.length() - num_insertions_5_prime);

   // main region
   int error_index_adjust(0);

   std::unordered_map<int, std::string> insertion_map_tmp(insertion_map);

   // iterate gaps in alignment1 (i.e. insertions)
   // ignore the gap at the 3'-end: add [ACGT] to the end of the regular expression
   // it_gap->position(): start index
   // it_gap->length()  : gap length + 1
   std::regex rx_gap("-+[ACGT]");

   for(auto it_gap = std::sregex_iterator(alignment1.begin(), alignment1.end(), rx_gap); it_gap != std::sregex_iterator(); ++it_gap) {
      int index = it_gap->position() + error_index_adjust;

      // insertion in the same position of an original read
      if (insertion_map.find(index) != insertion_map.end()) {
         // insertion length of the modified read <= insertion length of the original read
         if ((it_gap->length() - 1) <= (int)insertion_map_tmp[index].length()) {
            for (int it_each = it_gap->position(); it_each < (it_gap->position() + (it_gap->length() - 1)); it_each++) {
               base_org_read = insertion_map_tmp[index][it_each - it_gap->position()];
               base_mod_read = alignment2[it_each];

               // two insertions are same
               // the insertion in the original read is not modified
               // before
               //    ref     : A-A
               //    org read: AAA
               //    mod read: AAA
               if (base_mod_read == base_org_read) {
                  num_nnns_insertion_tmp++;
               }
               // two insertions are not same
               // the insertion in the original read is wrongly modified
               //    ref     : A-A
               //    org read: AAA
               //    mod read: ACA
               else {
                  num_nyns_insertion_tmp++;
               }
            }

            // this insertion became shorter
            int length_tmp(insertion_map_tmp[index].length() - (it_gap->length() - 1));
            num_nyys_insertion_tmp += length_tmp;

            alignment_score -= (gap_extension_penalty * length_tmp);
         }

         // insertion length of the modified read > insertion length of the original read
         else {
            for (std::size_t it_each = 0; it_each < insertion_map_tmp[index].length(); it_each++) {
               base_org_read = insertion_map_tmp[index][it_each];
               base_mod_read = alignment2[it_each + it_gap->position()];

               // two insertions are same
               // the insertion in the original read is not modified
               //    ref     : A-A
               //    org read: AAA
               //    mod read: AAA
               if (base_mod_read == base_org_read) {
                  num_nnns_insertion_tmp++;
               }
               // two insertions are not same
               // the insertion in the original read is wrongly modified
               //    ref     : A-A
               //    org read: AAA
               //    mod read: ACA
               else {
                  num_nyns_insertion_tmp++;
               }
            }

            // this insertion became longer
            int length_tmp((it_gap->length() - 1) - insertion_map_tmp[index].length());
            num_yyns_insertion_tmp += length_tmp;

            alignment_score += (gap_extension_penalty * length_tmp);
         }

         insertion_map_tmp.erase(index);
      }
      // newly generated insertion
      else {
         num_yyns_insertion_tmp += (it_gap->length() - 1);

         alignment_score     += (gap_opening_penalty + gap_extension_penalty * (it_gap->length() - 1));
         alignment_score_new += (gap_opening_penalty + gap_extension_penalty * (it_gap->length() - 1));
      }

      error_index_adjust += (it_gap->length() - 1);
   }

   //--------------------------------------------------
   // generate temporary strings for counting deletions and substitutions
   // alignment1_tmp: removing all "-"s from alignment1
   // alignment2_tmp: removing all the inserted bases from alignment2
   //--------------------------------------------------
   std::string alignment1_tmp;
   std::string alignment2_tmp;

   for (std::size_t it_base = 0; it_base < alignment1.length(); it_base++) {
      if (alignment1[it_base] == '-') {
         how_many_insertions++;
      }
      else {
         alignment1_tmp.push_back(alignment1[it_base]);
         alignment2_tmp.push_back(alignment2[it_base]);
      }
   }

   // original_read_tmp: original_read without insertions
   // prev_end: 0-based
   std::string original_read_tmp("");

   // there are substitutions or deletions in the original read
   if ((deletion_map.size() > 0) || ((substitution_err_map.size() > 0))) {
      for (std::size_t it_base = 1; it_base <= alignment1_tmp.length(); it_base++) {
         if (deletion_map.find(it_base) != deletion_map.end()) {
            original_read_tmp.push_back('-');
         }
         else if (substitution_err_map.find(it_base) != substitution_err_map.end()) {
            original_read_tmp.push_back(substitution_err_map[it_base]);
         }
         else {
            original_read_tmp.push_back(alignment1_tmp[it_base - 1]);
         }
      }
   }
   // there are no substitutions or deletions in the original read
   else {
      original_read_tmp = alignment1_tmp;
   }

   // check length of sequences without insertions
   int alignment1_tmp_length(alignment1_tmp.length());

   // count the number of insertions at the 5'-end
   rx_start_gap = "^-+";

   if (std::regex_search(alignment2_tmp, smatch1, rx_start_gap)) {
      num_deletions_5_prime_tmp = smatch1[0].length();
   }
   else {
      num_deletions_5_prime_tmp = 0;
   }

   // count the number of insertions at the 3'-end
   std::regex rx_end_gap("-+$");

   if (std::regex_search(alignment2_tmp, smatch1, rx_end_gap)) {
      num_deletions_3_prime_tmp = smatch1[0].length();
   }
   else {
      num_deletions_3_prime_tmp = 0;
   }

   // count corrected/trimmed insertions
   for (auto it_insertion = insertion_map_tmp.begin(); it_insertion != insertion_map_tmp.end(); it_insertion++) {
      int length_tmp(it_insertion->second.length());

      if (is_trimmed) {
         // corrected insertion
         if (it_insertion->first <= (alignment1_tmp_length - num_deletions_3_prime_tmp)) {
            num_nyys_insertion_tmp += length_tmp;
         }
         // trimmed insertion
         else {
            num_nyys_insertion_trim_tmp += length_tmp;
         }
      }
      else {
         num_nyys_insertion_tmp += length_tmp;
      }

      alignment_score -= (gap_opening_penalty + gap_extension_penalty * length_tmp);
   }

   //--------------------------------------------------
   // SUBSTITUTION & DELETION
   //--------------------------------------------------
   // $a = "000---000";
   // /[\-]+/
   // $-[0]: 3 (0-based)
   // $+[0]: 6 (not 5)
   //--------------------------------------------------

   //----------------------------------------------------------------------
   // substitutions/deletions in boundaries
   //----------------------------------------------------------------------
   // trimmed read
   if (is_trimmed) {
      // deletion
      for (auto it_deletion = deletion_map.begin(); it_deletion != deletion_map.end(); it_deletion++) {
         // 5'-end
         if (it_deletion->first <= num_deletions_5_prime_tmp) {
            num_nyys_deletion_trim_tmp++;

            if (it_deletion->first == 1) {
               alignment_score -= (gap_opening_penalty + gap_extension_penalty);
            }
            else {
               if (original_read_tmp[it_deletion->first - 1] == '-') {
                  alignment_score -= gap_extension_penalty;
               }
               else {
                  alignment_score -= (gap_opening_penalty + gap_extension_penalty);
               }
            }
         }
         // 3'-end
         else if (it_deletion->first > (alignment1_tmp_length - num_deletions_3_prime_tmp)) {
            num_nyys_deletion_trim_tmp++;

            if (it_deletion->first == alignment1_tmp_length) {
               alignment_score -= (gap_opening_penalty + gap_extension_penalty);
            }
            else {
               if (original_read_tmp[it_deletion->first - 1] == '-') {
                  alignment_score -= gap_extension_penalty;
               }
               else {
                  alignment_score -= (gap_opening_penalty + gap_extension_penalty);
               }
            }
         }
      }

      // substitution
      // you don't have to worry about the order of the keys
      // error_index will be sorted later using the 1st and 2nd columns
      for (auto it_substitution = substitution_org_map.begin(); it_substitution != substitution_org_map.end(); it_substitution++) {
         // 5'-end
         if (it_substitution->first <= num_deletions_5_prime_tmp) {
            num_nyys_substitution_trim_tmp++;

            alignment_score -= mismatch_penalty;

            // update position vectors
            position_vector_tmp[it_substitution->first]++;
            corrected_position_vector_tmp[it_substitution->first]++;

            if (is_detail) {
               int current_index;

               char org_base;
               char err_base;

               // calculate the index of the error in the reference sequence
               if (strand == "+") {
                  // both start_index and it_substitution->first are 1-based
                  current_index = start_index + it_substitution->first - 1;

                  org_base = it_substitution->second;
                  err_base = substitution_err_map[it_substitution->first];
               }
               else if (strand == "-") {
                  // both end_index and it_substitution->first are 1-based
                  // end_index of the "-" strand is not changed even if a corrected read is longer than an original read
                  current_index = end_index - it_substitution->first + 1;

                  org_base = complement_base(it_substitution->second);
                  err_base = complement_base(substitution_err_map[it_substitution->first]);
               }
               else {
                  std::cout << "\nERROR: Illegal strand " << strand << "\n\n";
                  exit(EXIT_FAILURE);
               }

               error_index_tmp += 
                                 (std::to_string(ref_seq_index) + " " +
                                  std::to_string(current_index) + " " +
                                  org_base                      + " " +
                                  err_base                      + " " +
                                  "Y CASE05"                    + " " +
                                  read_name                     + " " +
                                  std::to_string(start_index)   + " " +
                                  std::to_string(end_index)     + "\n"
                                 );

            }
         }
         // 3'-end
         else if (it_substitution->first > (alignment1_tmp_length - num_deletions_3_prime_tmp)) {
            num_nyys_substitution_trim_tmp++;

            // $in_mismatch_penalty = -1
            alignment_score -= mismatch_penalty;

            // update position vectors
            position_vector_tmp[it_substitution->first]++;
            corrected_position_vector_tmp[it_substitution->first]++;

            if (is_detail) {
               int current_index;

               char org_base;
               char err_base;

               // calculate the index of the error in the reference sequence
               if (strand == "+") {
                  // both start_index and it_substitution->first are 1-based
                  current_index = start_index + it_substitution->first - 1;

                  org_base = it_substitution->second;
                  err_base = substitution_err_map[it_substitution->first];
               }
               else if (strand == "-") {
                  // both end_index and it_substitution->first are 1-based
                  // end_index of the "-" strand is not changed even if a corrected read is longer than an original read
                  current_index = end_index - it_substitution->first + 1;

                  org_base = complement_base(it_substitution->second);
                  err_base = complement_base(substitution_err_map[it_substitution->first]);
               }
               else {
                  std::cout << "\nERROR: Illegal strand " << strand << "\n\n";
                  exit(EXIT_FAILURE);
               }

               error_index_tmp += 
                                 (std::to_string(ref_seq_index) + " " +
                                  std::to_string(current_index) + " " +
                                  org_base                      + " " +
                                  err_base                      + " " +
                                  "Y CASE06"                    + " " +
                                  read_name                     + " " +
                                  std::to_string(start_index)   + " " +
                                  std::to_string(end_index)     + "\n"
                                 );
            }
         }
      }
   }
   // non-trimmed read
   else {
      // 5'-end
      // deletion
      for (auto it_deletion = deletion_map.begin(); it_deletion != deletion_map.end(); it_deletion++) {
         if (it_deletion->first <= num_deletions_5_prime_tmp) {
            num_nyys_deletion_tmp++;

            if (it_deletion->first == 1) {
               alignment_score -= gap_opening_penalty + gap_extension_penalty;
            }
            else {
               if (original_read_tmp[it_deletion->first - 1] == '-') {
                  alignment_score -= gap_extension_penalty;
               }
               else {
                  alignment_score -= gap_opening_penalty + gap_extension_penalty;
               }
            }
         }
      }

      // substitution
      // a modified read starts from the second index
      // REF: AAAAA
      // ORG: CAAAA
      // MOD: -AAAA
      for (auto it_substitution = substitution_org_map.begin(); it_substitution != substitution_org_map.end(); it_substitution++) {
         // 5'-end
         // 3'-end will be handled later
         if (it_substitution->first <= num_deletions_5_prime_tmp) {
            num_nyys_substitution_tmp++;

            alignment_score -= mismatch_penalty;

            // update position vectors
            position_vector_tmp[it_substitution->first]++;
            corrected_position_vector_tmp[it_substitution->first]++;

            if (is_detail) {
               int current_index;

               char org_base;
               char err_base;

               // calculate the index of the error in the reference sequence
               if (strand == "+") {
                  // both start_index and it_substitution->first are 1-based
                  current_index = start_index + it_substitution->first - 1;

                  org_base = it_substitution->second;
                  err_base = substitution_err_map[it_substitution->first];
               }
               else if (strand == "-") {
                  // both end_index and it_substitution->first are 1-based
                  // end_index of the "-" strand is not changed even if a corrected read is longer than an original read
                  current_index = end_index - it_substitution->first + 1;

                  org_base = complement_base(it_substitution->second);
                  err_base = complement_base(substitution_err_map[it_substitution->first]);
               }
               else {
                  std::cout << "\nERROR: Illegal strand " << strand << "\n\n";
                  exit(EXIT_FAILURE);
               }

               error_index_tmp += 
                                 (std::to_string(ref_seq_index) + " " +
                                  std::to_string(current_index) + " " +
                                  org_base                      + " " +
                                  err_base                      + " " +
                                  "Y CASE07"                    + " " +
                                  read_name                     + " " +
                                  std::to_string(start_index)   + " " +
                                  std::to_string(end_index)     + "\n" 
                                 );
            }
         }
      }
   }

   //--------------------------------------------------
   // main region
   //--------------------------------------------------
   char base_org_read_prev = 'A';
   char base_mod_read_prev = 'A';

   // trimmed read
   if (is_trimmed) {
      // find the end index
      // how_many_insertions: 1
      // read_length: 10

      //      000000000011111
      //      012345678901234
      // REF: AAAAAAAAAAAAAAA
      // MOD:   AAAAAAA-AAA
      // end_index: 11

      // where_read_ends: 0-based
      int where_read_ends(-1);
      int sum_tmp(0);

      for (std::size_t it_base = num_deletions_5_prime_tmp; it_base < alignment1.length(); it_base++) {
         // alignment2_tmp: A/C/G/T
         if (alignment2_tmp[it_base] != '-') {
            sum_tmp++;

            if (sum_tmp == (corrected_read_length - how_many_insertions)) {
               where_read_ends = it_base;
               break;
            }
         }
      }

      // where_read_ends was not found
      if (where_read_ends == -1) {
         std::cout << "\nERROR: " << alignment2 << "\n\n";
         exit(EXIT_FAILURE);
      }

      for (int it_base = num_deletions_5_prime_tmp; it_base <= where_read_ends; it_base++) {
         // REF(base_ref_seq) : A/C/G/T
         // ORG(base_org_read): A/C/G/T/-
         // MOD(base_mod_read): A/C/G/T/-
         base_org_read = original_read_tmp[it_base];
         base_ref_seq  = alignment1_tmp[it_base];
         base_mod_read = alignment2_tmp[it_base];

         // ORG: -
         if (base_org_read == '-') {
            // untouched deletion
            // ORG: -
            // MOD: -
            if (base_mod_read == '-') {
               num_nnns_deletion_tmp++;
            }
            // ORG: -
            // MOD: A/C/G/T
            else {
               // the deletion was correctly modified
               // REF: A/C/G/T
               // ORG: -
               // MOD: REF
               if (base_mod_read == base_ref_seq) {
                  num_nyys_deletion_tmp++;

                  if (base_org_read_prev == '-') {
                     alignment_score -= gap_extension_penalty;
                  }
                  else {
                     alignment_score -= (gap_opening_penalty + gap_extension_penalty);
                  }

                  // substitution in the same position
                  if (substitution_org_map.find(it_base + 1) != substitution_org_map.end()) {
                     num_nyys_substitution_tmp++;

                     alignment_score -= mismatch_penalty;

                     // update position vectors
                     position_vector_tmp[it_base + 1]++;
                     corrected_position_vector_tmp[it_base + 1]++;

                     if (is_detail) {
                        int current_index;

                        char org_base;
                        char err_base;

                        // calculate the index of the error in the reference sequence
                        if (strand == "+") {
                           // both start_index and it_base are 1-based
                           current_index = start_index + it_base;

                           org_base = substitution_org_map[it_base + 1];
                           err_base = substitution_err_map[it_base + 1];
                        }
                        else if (strand == "-") {
                           // both end_index and it_base are 1-based
                           // end_index of the "-" strand is not changed even if a corrected read is longer than an original read
                           current_index = end_index - it_base;

                           org_base = complement_base(substitution_org_map[it_base + 1]);
                           err_base = complement_base(substitution_err_map[it_base + 1]);
                        }
                        else {
                           std::cout << "\nERROR: Illegal strand " << strand << "\n\n";
                           exit(EXIT_FAILURE);
                        }

                        error_index_tmp += 
                                          (std::to_string(ref_seq_index) + " " +
                                           std::to_string(current_index) + " " +
                                           org_base                      + " " +
                                           err_base                      + " " +
                                           "Y CASE08"                    + " " +
                                           read_name                     + " " +
                                           std::to_string(start_index)   + " " +
                                           std::to_string(end_index)     + "\n" 
                                          );
                     }
                  }
               }
               // the deletion was removed
               // but a new substitution was introduced
               // REF: A/C/G/T
               // ORG: -
               // MOD: != REF
               else {
                  num_nyys_deletion_tmp++;
                  //num_yyns_substitution_tmp++;

                  if (base_org_read_prev == '-') {
                     alignment_score -= gap_extension_penalty;
                  }
                  else {
                     alignment_score -= (gap_opening_penalty + gap_extension_penalty);
                  }

                  //alignment_score -= mismatch_penalty;
               }
            }
         }
         // ORG: A/C/G/T
         else {
            // newly generated deletion
            // ORG: A/C/G/T
            // MOD: -
            if (base_mod_read == '-') {
               // REF: A/C/G/T
               // ORG: REF
               // MOD: -
               if (base_org_read == base_ref_seq) {
                  //num_yyns_substitution_tmp++;
               }
               // REF: A/C/G/T
               // ORG: != REF
               // MOD: -
               else {
                  num_from_substitution_to_deletion_tmp++;
               }

               num_yyns_deletion_tmp++;

               if (base_mod_read_prev == '-') {
                  alignment_score     += gap_extension_penalty;
                  alignment_score_new += gap_extension_penalty;
               }
               else {
                  alignment_score     += (gap_opening_penalty + gap_extension_penalty);
                  alignment_score_new += (gap_opening_penalty + gap_extension_penalty);
               }
            }
            //--------------------------------------------------
            // SUBSTITUTION
            //--------------------------------------------------
            // ORG: A/C/G/T
            // MOD: A/C/G/T
            else {
               // base_org_read: error-free base
               // REF: A/C/G/T
               // ORG: REF
               if (base_org_read == base_ref_seq) {
                  // no modification was made
                  // REF: A/C/G/T
                  // ORG: REF
                  // MOD: REF
                  if (base_mod_read == base_ref_seq) {
                     num_ynys_substitution_tmp++;
                  }
                  // newly generated substitution
                  // REF: A/C/G/T
                  // ORG: REF
                  // MOD: != REF
                  else {
                     num_yyns_substitution_tmp++;

                     alignment_score     += mismatch_penalty;
                     alignment_score_new += mismatch_penalty;
                  }
               }
               // base_org_read: erroneous base
               // REF: A/C/G/T
               // ORG: != REF
               else {
                  // no modification was made
                  // base_org_read: erroneous base
                  // REF: A/C/G/T
                  // ORG: != REF
                  // MOD: ORG
                  if (base_mod_read == base_org_read) {
                     num_nnns_substitution_tmp++;

                     // update position vectors
                     position_vector_tmp[it_base + 1]++;

                     if (is_detail) {
                        int current_index;

                        char org_base;
                        char err_base;

                        // calculate the index of the error in the reference sequence
                        if (strand == "+") {
                           // both start_index and it_base are 1-based
                           current_index = start_index + it_base;

                           org_base = substitution_org_map[it_base + 1];
                           err_base = substitution_err_map[it_base + 1];
                        }
                        else if (strand == "-") {
                           // both end_index and it_base are 1-based
                           // end_index of the "-" strand is not changed even if a corrected read is longer than an original read
                           current_index = end_index - it_base;

                           org_base = complement_base(substitution_org_map[it_base + 1]);
                           err_base = complement_base(substitution_err_map[it_base + 1]);
                        }
                        else {
                           std::cout << "\nERROR: Illegal strand " << strand << "\n\n";
                           exit(EXIT_FAILURE);
                        }

                        error_index_tmp += 
                                          (std::to_string(ref_seq_index) + " " +
                                           std::to_string(current_index) + " " +
                                           org_base                      + " " +
                                           err_base                      + " " +
                                           "N CASE09"                    + " " +
                                           read_name                     + " " +
                                           std::to_string(start_index)   + " " +
                                           std::to_string(end_index)     + "\n" 
                                          );
                     }
                  }

                  // a modification was made
                  // REF: A/C/G/T
                  // ORG: != REF
                  // MOD: != ORG
                  else {
                     // correct modification
                     // REF: A/C/G/T
                     // ORG: != REF
                     // MOD: REF
                     if (base_mod_read == base_ref_seq) {
                        num_nyys_substitution_tmp++;

                        alignment_score -= mismatch_penalty;

                        // update position vectors
                        position_vector_tmp[it_base + 1]++;
                        corrected_position_vector_tmp[it_base + 1]++;

                        if (is_detail) {
                           int current_index;

                           char org_base;
                           char err_base;

                           // calculate the index of the error in the reference sequence
                           if (strand == "+") {
                              // both start_index and it_base are 1-based
                              current_index = start_index + it_base;

                              org_base = substitution_org_map[it_base + 1];
                              err_base = substitution_err_map[it_base + 1];
                           }
                           else if (strand == "-") {
                              // both end_index and it_base are 1-based
                              // end_index of the "-" strand is not changed even if a corrected read is longer than an original read
                              current_index = end_index - it_base;

                              org_base = complement_base(substitution_org_map[it_base + 1]);
                              err_base = complement_base(substitution_err_map[it_base + 1]);
                           }
                           else {
                              std::cout << "\nERROR: Illegal strand " << strand << "\n\n";
                              exit(EXIT_FAILURE);
                           }

                           error_index_tmp += 
                                             (std::to_string(ref_seq_index) + " " +
                                              std::to_string(current_index) + " " +
                                              org_base                      + " " +
                                              err_base                      + " " +
                                              "Y CASE10"                    + " " +
                                              read_name                     + " " +
                                              std::to_string(start_index)   + " " +
                                              std::to_string(end_index)     + "\n" 
                                             );
                        }
                     }
                     // wrong  modification
                     // REF: A/C/G/T
                     // ORG: != REF
                     // MOD: != ORG & != REF
                     else {
                        num_nyns_substitution_tmp++;

                        // this is worse than nnn
                        alignment_score += mismatch_penalty;

                        // update position vectors
                        position_vector_tmp[it_base + 1]++;

                        if (is_detail) {
                           int current_index;

                           char org_base;
                           char err_base;

                           // calculate the index of the error in the reference sequence
                           if (strand == "+") {
                              // both $start_index and it_base are 1-based
                              current_index = start_index + it_base;

                              org_base = substitution_org_map[it_base + 1];
                              err_base = substitution_err_map[it_base + 1];
                           }
                           else if (strand == "-") {
                              // both end_index and it_base are 1-based
                              // end_index of the "-" strand is not changed even if a corrected read is longer than an original read
                              current_index = end_index - it_base;

                              org_base = complement_base(substitution_org_map[it_base + 1]);
                              err_base = complement_base(substitution_err_map[it_base + 1]);
                           }
                           else {
                              std::cout << "\nERROR: Illegal strand " << strand << "\n\n";
                              exit(EXIT_FAILURE);
                           }

                           error_index_tmp += 
                                             (std::to_string(ref_seq_index) + " " +
                                              std::to_string(current_index) + " " +
                                              org_base                      + " " +
                                              err_base                      + " " +
                                              "N CASE11"                    + " " +
                                              read_name                     + " " +
                                              std::to_string(start_index)   + " " +
                                              std::to_string(end_index)     + "\n" 
                                             );
                        }
                     }
                  }
               }
            }
         }

         base_org_read_prev = base_org_read;
         base_mod_read_prev = base_mod_read;
      }
   }

   // not trimmed read
   // in this case, $alignment1_tmp does not have the information for the last base
   // outer_3_end should be used
   // outer_3_end: ATA
   // REF    : ACCCCC-
   // MOD    : -CCCCCA
   // MOD REF: ACCCCCATA
   else {
      std::string alignment1_tmp_tmp(std::regex_replace(alignment1_tmp, rx_end_gap, ""));
      alignment1_tmp_tmp += outer_3_end;

      // outer_3_end is not long enough
      // add "-" to the end
      if (alignment1_tmp_tmp.length() < alignment2_tmp.length()) {
         alignment1_tmp_tmp += std::string(alignment2_tmp.length() - alignment1_tmp_tmp.length(), '-');
      }

      // find the end index
      // how_many_insertions: 1
      // read_length: 10
      //      000000000011111
      //      012345678901234
      // REF: AAAAAAAAAAAAAAA
      // MOD: --AAAAAAA-AAA--
      // $where_read_ends: 11
      int where_read_ends(-1);
      int sum_tmp(0);

      for (int it_base = num_deletions_5_prime_tmp; it_base < alignment1_tmp_length; it_base++) {
         // alignment2_tmp: A/C/G/T
         if (alignment2_tmp[it_base] != '-') {
            sum_tmp++;

            if (sum_tmp == (read_length - how_many_insertions)) {
               where_read_ends = it_base;
               break;
            }
         }
      }

      // $where_read_ends was not found
      if (where_read_ends == -1) {
         std::cout << "\nERROR: The read end point is not found\n\n";
         exit(EXIT_FAILURE);
      }

      // count the number of substitutions out of $where_read_ends
      for (auto it_substitution = substitution_org_map.begin(); it_substitution != substitution_org_map.end(); it_substitution++) {
         if (it_substitution->first > (where_read_ends + 1)) {
            num_nyys_substitution_tmp++;

            alignment_score -= mismatch_penalty;

            // update position vector
            position_vector_tmp[it_substitution->first]++;
            corrected_position_vector_tmp[it_substitution->first]++;

            if (is_detail) {
               int current_index;

               char org_base;
               char err_base;

               // calculate the index of the error in the reference sequence
               if (strand == "+") {
                  // both start_index and it_substitution->first are 1-based
                  current_index = start_index + it_substitution->first - 1;

                  org_base = it_substitution->second;
                  err_base = substitution_err_map[it_substitution->first];
               }
               else if (strand == "-") {
                  // both end_index and it_substitution->first are 1-based
                  // end_index of the "-" strand is not changed even if a corrected read is longer than an original read
                  current_index = end_index - it_substitution->first + 1;

                  org_base = complement_base(it_substitution->second);
                  err_base = complement_base(substitution_err_map[it_substitution->first]);
               }
               else {
                  std::cout << "\nERROR: Illegal strand " << strand << "\n\n";
                  exit(EXIT_FAILURE);
               }

               error_index_tmp += 
                                 (std::to_string(ref_seq_index) + " " +
                                  std::to_string(current_index) + " " +
                                  org_base                      + " " +
                                  err_base                      + " " +
                                  "Y CASE12"                    + " " +
                                  read_name                     + " " +
                                  std::to_string(start_index)   + " " +
                                  std::to_string(end_index)     + "\n" 
                                 );
            }
         }
      }

      // count the number of deletions out of $where_read_ends
      for (auto it_deletion = deletion_map.begin(); it_deletion != deletion_map.end(); it_deletion++) {
         if (it_deletion->first > (where_read_ends + 1)) {
            num_nyys_deletion_tmp++;

            alignment_score -= gap_extension_penalty;
         }
      }

      // main routine
      for (int it_base = num_deletions_5_prime_tmp; it_base <= where_read_ends; it_base++) {
         // REF(base_ref_seq) : A/C/G/T
         // ORG(base_org_read): A/C/G/T/-
         // MOD(base_mod_read): A/C/G/T/-
         base_org_read = original_read_tmp[it_base];
         base_ref_seq  = alignment1_tmp_tmp[it_base];
         base_mod_read = alignment2_tmp[it_base];

         // ORG: -
         if (base_org_read == '-') {
            // untouched deletion
            // ORG: -
            // MOD: -
            if (base_mod_read == '-') {
               num_nnns_deletion_tmp++;
            }
            // ORG: -
            // MOD: A/C/G/T
            else {
               // the deletion was correctly modified
               // REF: A/C/G/T
               // ORG: -
               // MOD: REF
               if (base_mod_read == base_ref_seq) {
                  num_nyys_deletion_tmp++;

                  if (base_org_read_prev == '-') {
                     alignment_score -= gap_extension_penalty;
                  }
                  else {
                     alignment_score -= gap_opening_penalty + gap_extension_penalty;
                  }

                  // substitution in the same position
                  if (substitution_org_map.find(it_base + 1) != substitution_org_map.end()) {
                     num_nyys_substitution_tmp++;

                     alignment_score -= mismatch_penalty;

                     // update position vectors
                     position_vector_tmp[it_base + 1]++;
                     corrected_position_vector_tmp[it_base + 1]++;

                     if (is_detail) {
                        int current_index;

                        char org_base;
                        char err_base;

                        // calculate the index of the error in the reference sequence
                        if (strand == "+") {
                           // both start_index and it_base are 1-based
                           current_index = start_index + it_base;

                           org_base = substitution_org_map[it_base + 1];
                           err_base = substitution_err_map[it_base + 1];
                        }
                        else if (strand == "-") {
                           // both end_index and it_base are 1-based
                           // end_index of the "-" strand is not changed even if a corrected read is longer than an original read
                           current_index = end_index - it_base;

                           org_base = complement_base(substitution_org_map[it_base + 1]);
                           err_base = complement_base(substitution_err_map[it_base + 1]);
                        }
                        else {
                           std::cout << "\nERROR: Illegal strand " << strand << "\n\n";
                           exit(EXIT_FAILURE);
                        }

                        error_index_tmp += 
                                          (std::to_string(ref_seq_index) + " " +
                                           std::to_string(current_index) + " " +
                                           org_base                      + " " +
                                           err_base                      + " " +
                                           "Y CASE13"                    + " " +
                                           read_name                     + " " +
                                           std::to_string(start_index)   + " " +
                                           std::to_string(end_index)     + "\n" 
                                          );
                     }
                  }
               }
               // the deletion was removed
               // but a new substitution was introduced
               // REF: A/C/G/T
               // ORG: -
               // MOD: != REF
               else {
                  num_nyys_deletion_tmp++;
                  //num_yyns_substitution_tmp++;

                  if (base_org_read_prev == '-') {
                     alignment_score -= gap_extension_penalty;
                  }
                  else {
                     alignment_score -= gap_opening_penalty + gap_extension_penalty;
                  }

                  //alignment_score -= mismatch_penalty;
               }
            }
         }
         // ORG: A/C/G/T
         else {
            // newly generated deletion
            // ORG: A/C/G/T
            // MOD: -
            if (base_mod_read == '-') {
               // REF: A/C/G/T
               // ORG: REF
               // MOD: -
               if (base_org_read == base_ref_seq) {
                  //num_yyns_substitution_tmp++;
               }
               // REF: A/C/G/T
               // ORG: != REF
               // MOD: -
               else {
                  num_from_substitution_to_deletion_tmp++;
               }

               num_yyns_deletion_tmp++;

               if (base_mod_read_prev == '-') {
                  alignment_score     += gap_extension_penalty;
                  alignment_score_new += gap_extension_penalty;
               }
               else {
                  alignment_score     += gap_opening_penalty + gap_extension_penalty;
                  alignment_score_new += gap_opening_penalty + gap_extension_penalty;
               }
            }
            //--------------------------------------------------
            // SUBSTITUTION
            //--------------------------------------------------
            // ORG: A/C/G/T
            // MOD: A/C/G/T
            else {
               // base_org_read: error-free base
               // REF: A/C/G/T
               // ORG: REF
               if (base_org_read == base_ref_seq) {
                  // no modification was made
                  // REF: A/C/G/T
                  // ORG: REF
                  // MOD: REF
                  if (base_mod_read == base_ref_seq) {
                     num_ynys_substitution_tmp++;
                  }
                  // newly generated substitution
                  // REF: A/C/G/T
                  // ORG: REF
                  // MOD: != REF
                  else {
                     num_yyns_substitution_tmp++;

                     alignment_score     += mismatch_penalty;
                     alignment_score_new += mismatch_penalty;
                  }
               }
               // base_org_read: erroneous base
               // REF: A/C/G/T
               // ORG: != REF
               else {
                  // no modification was made
                  // base_org_read: erroneous base
                  // REF: A/C/G/T
                  // ORG: != REF
                  // MOD: ORG
                  if (base_mod_read == base_org_read) {
                     num_nnns_substitution_tmp++;

                     // update @array_position*
                     position_vector_tmp[it_base + 1]++;

                     if (is_detail) {
                        int current_index;

                        char org_base;
                        char err_base;

                        // calculate the index of the error in the reference sequence
                        if (strand == "+") {
                           // both start_index and it_base are 1-based
                           current_index = start_index + it_base;

                           org_base = substitution_org_map[it_base + 1];
                           err_base = substitution_err_map[it_base + 1];
                        }
                        else if (strand == "-") {
                           // both end_index and it_base are 1-based
                           // end_index of the '-' strand is not changed even if a corrected read is longer than an original read
                           current_index = end_index - it_base;

                           org_base = complement_base(substitution_org_map[it_base + 1]);
                           err_base = complement_base(substitution_err_map[it_base + 1]);
                        }
                        else {
                           std::cout << "\nERROR: Illegal strand " << strand << "\n\n";
                           exit(EXIT_FAILURE);
                        }

                        error_index_tmp += 
                                          (std::to_string(ref_seq_index) + " " +
                                           std::to_string(current_index) + " " +
                                           org_base                      + " " +
                                           err_base                      + " " +
                                           "N CASE14"                    + " " +
                                           read_name                     + " " +
                                           std::to_string(start_index)   + " " +
                                           std::to_string(end_index)     + "\n" 
                                          );
                     }
                  }

                  // a modification was made
                  // REF: A/C/G/T
                  // ORG: != REF
                  // MOD: != ORG
                  else {
                     // correct modification
                     // REF: A/C/G/T
                     // ORG: != REF
                     // MOD: REF
                     if (base_mod_read == base_ref_seq) {
                        num_nyys_substitution_tmp++;

                        alignment_score -= mismatch_penalty;

                        // update position vectors
                        position_vector_tmp[it_base + 1]++;
                        corrected_position_vector_tmp[it_base + 1]++;

                        if (is_detail) {
                           int current_index;

                           char org_base;
                           char err_base;

                           // calculate the index of the error in the reference sequence
                           if (strand == "+") {
                              // both start_index and it_base are 1-based
                              current_index = start_index + it_base;

                              org_base = substitution_org_map[it_base + 1];
                              err_base = substitution_err_map[it_base + 1];
                           }
                           else if (strand == "-") {
                              // both end_index and it_base are 1-based
                              // end_index of the '-' strand is not changed even if a corrected read is longer than an original read
                              current_index = end_index - it_base;

                              org_base = complement_base(substitution_org_map[it_base + 1]);
                              err_base = complement_base(substitution_err_map[it_base + 1]);
                           }
                           else {
                              std::cout << "\nERROR: Illegal strand " << strand << "\n\n";
                              exit(EXIT_FAILURE);
                           }

                           error_index_tmp += 
                                             (std::to_string(ref_seq_index) + " " +
                                              std::to_string(current_index) + " " +
                                              org_base                      + " " +
                                              err_base                      + " " +
                                              "Y CASE15"                    + " " +
                                              read_name                     + " " +
                                              std::to_string(start_index)   + " " +
                                              std::to_string(end_index)     + "\n" 
                                             );
                        }
                     }
                     // wrong  modification
                     // REF: A/C/G/T
                     // ORG: != REF
                     // MOD: != ORG
                     // MOD: != REF
                     else {
                        num_nyns_substitution_tmp++;

                        alignment_score += mismatch_penalty;

                        // update position vectors
                        position_vector_tmp[it_base + 1]++;

                        if (is_detail) {
                           int current_index;

                           char org_base;
                           char err_base;

                           // calculate the index of the error in the reference sequence
                           if (strand == "+") {
                              // both start_index and it_base are 1-based
                              current_index = start_index + it_base;

                              org_base = substitution_org_map[it_base + 1];
                              err_base = substitution_err_map[it_base + 1];
                           }
                           else if (strand == "-") {
                              // both end_index and it_base are 1-based
                              // end_index of the '-' strand is not changed even if a corrected read is longer than an original read
                              current_index = end_index - it_base;

                              org_base = complement_base(substitution_org_map[it_base + 1]);
                              err_base = complement_base(substitution_err_map[it_base + 1]);
                           }
                           else {
                              std::cout << "\nERROR: Illegal strand " << strand << "\n\n";
                              exit(EXIT_FAILURE);
                           }

                           error_index_tmp += 
                                             (std::to_string(ref_seq_index) + " " +
                                              std::to_string(current_index) + " " +
                                              org_base                      + " " +
                                              err_base                      + " " +
                                              "N CASE16"                    + " " +
                                              read_name                     + " " +
                                              std::to_string(start_index)   + " " +
                                              std::to_string(end_index)     + "\n" 
                                             );
                        }
                     }
                  }
               }
            }
         }

         base_org_read_prev = base_org_read;
         base_mod_read_prev = base_mod_read;
      }
   }

   // if
   //    1. a new total alignment score is higher than or equal to the current best and
   //    2. a new newly generated alignment score is higher than the current best
   // then the current best values are updated
   // local: on a node
   if (alignment_score >= alignment_score_best) {
      if (alignment_score_new > alignment_score_new_best) {
         num_yyns_substitution_local_best = num_yyns_substitution_tmp;
         num_ynys_substitution_local_best = num_ynys_substitution_tmp;
         num_nyys_substitution_local_best = num_nyys_substitution_tmp;
         num_nyns_substitution_local_best = num_nyns_substitution_tmp;
         num_nnns_substitution_local_best = num_nnns_substitution_tmp;

         num_yyns_insertion_local_best = num_yyns_insertion_tmp;
         //num_ynys_insertion_local_best = num_ynys_insertion_tmp;
         num_nyys_insertion_local_best = num_nyys_insertion_tmp;
         num_nyns_insertion_local_best = num_nyns_insertion_tmp;
         num_nnns_insertion_local_best = num_nnns_insertion_tmp;

         num_yyns_deletion_local_best = num_yyns_deletion_tmp;
         //num_ynys_deletion_local_best = num_ynys_deletion_tmp;
         num_nyys_deletion_local_best = num_nyys_deletion_tmp;
         num_nyns_deletion_local_best = num_nyns_deletion_tmp;
         num_nnns_deletion_local_best = num_nnns_deletion_tmp;

         num_from_substitution_to_deletion_local_best = num_from_substitution_to_deletion_tmp;

         num_nyys_substitution_trim_local_best = num_nyys_substitution_trim_tmp;

         num_nyys_insertion_trim_local_best = num_nyys_insertion_trim_tmp;

         num_nyys_deletion_trim_local_best = num_nyys_deletion_trim_tmp;

         num_deletions_5_prime_best = num_deletions_5_prime_tmp;
         num_deletions_3_prime_best = num_deletions_3_prime_tmp;

         alignment_score_best     = alignment_score;
         alignment_score_new_best = alignment_score_new;
         alignment_best           = alignment1_keep + "\n" + alignment2_keep + "\n";

         position_vector_local_best_tmp           = position_vector_tmp;
         corrected_position_vector_local_best_tmp = corrected_position_vector_tmp;

         if (is_detail) {
            error_index_best = error_index_tmp;
         }
      }
      else if (alignment_score_new == alignment_score_new_best) {
         if ((num_deletions_5_prime_tmp + num_deletions_3_prime_tmp) > (num_deletions_5_prime_best + num_deletions_3_prime_best)) {
            num_yyns_substitution_local_best = num_yyns_substitution_tmp;
            num_ynys_substitution_local_best = num_ynys_substitution_tmp;
            num_nyys_substitution_local_best = num_nyys_substitution_tmp;
            num_nyns_substitution_local_best = num_nyns_substitution_tmp;
            num_nnns_substitution_local_best = num_nnns_substitution_tmp;

            num_yyns_insertion_local_best = num_yyns_insertion_tmp;
            //num_ynys_insertion_local_best = num_ynys_insertion_tmp;
            num_nyys_insertion_local_best = num_nyys_insertion_tmp;
            num_nyns_insertion_local_best = num_nyns_insertion_tmp;
            num_nnns_insertion_local_best = num_nnns_insertion_tmp;

            num_yyns_deletion_local_best = num_yyns_deletion_tmp;
            //num_ynys_deletion_local_best = num_ynys_deletion_tmp;
            num_nyys_deletion_local_best = num_nyys_deletion_tmp;
            num_nyns_deletion_local_best = num_nyns_deletion_tmp;
            num_nnns_deletion_local_best = num_nnns_deletion_tmp;

            num_from_substitution_to_deletion_local_best = num_from_substitution_to_deletion_tmp;

            num_nyys_substitution_trim_local_best = num_nyys_substitution_trim_tmp;

            num_nyys_insertion_trim_local_best = num_nyys_insertion_trim_tmp;

            num_nyys_deletion_trim_local_best = num_nyys_deletion_trim_tmp;

            num_deletions_5_prime_best = num_deletions_5_prime_tmp;
            num_deletions_3_prime_best = num_deletions_3_prime_tmp;

            alignment_score_best     = alignment_score;
            alignment_score_new_best = alignment_score_new;
            alignment_best           = alignment1_keep + "\n" + alignment2_keep + "\n";

            position_vector_local_best_tmp           = position_vector_tmp;
            corrected_position_vector_local_best_tmp = corrected_position_vector_tmp;

            if (is_detail) {
               error_index_best = error_index_tmp;
            }
         }
      }
   }
}

//
// find_best_alignment
//
void find_best_alignment(int* position_vector_local_best, int* corrected_position_vector_local_best) {
   // find out the highest alignment score
   highest_score = match_matrix[matrix_size - 1];

   if (gap_1_matrix[matrix_size - 1] > highest_score) {
      highest_score = gap_1_matrix[matrix_size - 1];
   }

   if (gap_2_matrix[matrix_size - 1] > highest_score) {
      highest_score = gap_2_matrix[matrix_size - 1];
   }

   // initialize variables
   int index_x         = matrix_width - 1;
   int index_y         = matrix_height - 1;
   int matrix_index    = matrix_size - 1;
   int alignment_index = longest_alignment_length - 1;

   corrected_read_length = string2.length();
   outer_length_5_end    = outer_5_end.length();
   outer_length_3_end    = outer_3_end.length();

   std::string alignment1_tmp;
   std::string alignment2_tmp;

   alignment1_tmp.resize(longest_alignment_length);
   alignment2_tmp.resize(longest_alignment_length);

   char current_matrix;

   // the match array has the highest score
   if (match_matrix[matrix_size - 1] == highest_score) {
      current_matrix = 'M';
      traceback(alignment1_tmp, alignment2_tmp, index_x, index_y, current_matrix, alignment_index, matrix_index, highest_score);
   }

   // the match array has the highest score
   if ((gap_1_matrix[matrix_size - 1] == highest_score) && (too_many_candidates == false)) {
      current_matrix = '1';
      traceback(alignment1_tmp, alignment2_tmp, index_x, index_y, current_matrix, alignment_index, matrix_index, highest_score);
   }

   // the match array has the highest score
   if ((gap_2_matrix[matrix_size - 1] == highest_score) && (too_many_candidates == false)) {
      current_matrix = '2';
      traceback(alignment1_tmp, alignment2_tmp, index_x, index_y, current_matrix, alignment_index, matrix_index, highest_score);
   }

   // not too many candidates
   if (too_many_candidates == false) {
      // evalute each alignment
      for (unsigned it_vec = 0; it_vec < alignment1_vector.size(); it_vec++) {
         evaluate_each_alignment(alignment1_vector[it_vec], alignment2_vector[it_vec], position_vector_local_best, corrected_position_vector_local_best);
      }

      // update the histograms
      std::copy(position_vector_local_best_tmp.begin(),           position_vector_local_best_tmp.end(),           position_vector_local_best);
      std::copy(corrected_position_vector_local_best_tmp.begin(), corrected_position_vector_local_best_tmp.end(), corrected_position_vector_local_best);
   }
   // too many candidates
   else {
      num_not_evaluated_substitution = num_substitutions;
      num_not_evaluated_insertion    = num_insertions_unit;
      num_not_evaluated_deletion     = num_deletions;
   }

   delete[] match_matrix;
   delete[] gap_1_matrix;
   delete[] gap_2_matrix;
}

//
// give_random_alignment
//
void give_random_alignment() {
   // find out the highest alignment score
   highest_score = match_matrix[matrix_size - 1];

   if (gap_1_matrix[matrix_size - 1] > highest_score) {
      highest_score = gap_1_matrix[matrix_size - 1];
   }

   if (gap_2_matrix[matrix_size - 1] > highest_score) {
      highest_score = gap_2_matrix[matrix_size - 1];
   }

   // initialize variables
   int index_x         = matrix_width - 1;
   int index_y         = matrix_height - 1;
   int matrix_index    = matrix_size - 1;
   int alignment_index = longest_alignment_length - 1;

   corrected_read_length = string2.length();
   outer_length_5_end    = outer_5_end.length();
   outer_length_3_end    = outer_3_end.length();

   std::string alignment1_tmp;
   std::string alignment2_tmp;

   alignment1_tmp.resize(longest_alignment_length);
   alignment2_tmp.resize(longest_alignment_length);

   char current_matrix;

   // the match array has the highest score
   if (match_matrix[matrix_size - 1] == highest_score) {
      current_matrix = 'M';
      traceback(alignment1_tmp, alignment2_tmp, index_x, index_y, current_matrix, alignment_index, matrix_index, highest_score);
   }

   // the match array has the highest score
   if (alignment1_vector.size() == 0) {
      if ((gap_1_matrix[matrix_size - 1] == highest_score) && (too_many_candidates == false)) {
         current_matrix = '1';
         traceback(alignment1_tmp, alignment2_tmp, index_x, index_y, current_matrix, alignment_index, matrix_index, highest_score);
      }
   }

   // the match array has the highest score
   if (alignment1_vector.size() == 0) {
      if ((gap_2_matrix[matrix_size - 1] == highest_score) && (too_many_candidates == false)) {
         current_matrix = '2';
         traceback(alignment1_tmp, alignment2_tmp, index_x, index_y, current_matrix, alignment_index, matrix_index, highest_score);
      }
   }

   // no alignment is made
   // because there are too few overlaps between them
   if (alignment1_vector.size() == 0) {
      random_alignment1 = "";
      random_alignment2 = "";
   }
   else {
      random_alignment1 = alignment1_vector[0];
      random_alignment2 = alignment2_vector[0];
   }

   delete[] match_matrix;
   delete[] gap_1_matrix;
   delete[] gap_2_matrix;
}

//
// calculate_percent_similarity
//
void calculate_percent_similarity() {
   // main region
   // find out the highest alignment score
   highest_score = match_matrix[matrix_size - 1];

   if (gap_1_matrix[matrix_size - 1] > highest_score) {
      highest_score = gap_1_matrix[matrix_size - 1];
   }

   if (gap_2_matrix[matrix_size - 1] > highest_score) {
      highest_score = gap_2_matrix[matrix_size - 1];
   }

   // initialize variables
   int index_x         = matrix_width - 1;
   int index_y         = matrix_height - 1;
   int matrix_index    = matrix_size - 1;
   int alignment_index = longest_alignment_length - 1;

   corrected_read_length = string2.length();
   outer_length_5_end    = outer_5_end.length();
   outer_length_3_end    = outer_3_end.length();

   std::string alignment1_tmp;
   std::string alignment2_tmp;

   alignment1_tmp.resize(longest_alignment_length);
   alignment2_tmp.resize(longest_alignment_length);

   char current_matrix;

   // the match array has the highest score
   if (match_matrix[matrix_size - 1] == highest_score) {
      current_matrix = 'M';
      traceback(alignment1_tmp, alignment2_tmp, index_x, index_y, current_matrix, alignment_index, matrix_index, highest_score);
   }

   // the match array has the highest score
   if (alignment1_vector.size() == 0) {
      if ((gap_1_matrix[matrix_size - 1] == highest_score) && (too_many_candidates == false)) {
         current_matrix = '1';
         traceback(alignment1_tmp, alignment2_tmp, index_x, index_y, current_matrix, alignment_index, matrix_index, highest_score);
      }
   }

   // the match array has the highest score
   if (alignment1_vector.size() == 0) {
      if ((gap_2_matrix[matrix_size - 1] == highest_score) && (too_many_candidates == false)) {
         current_matrix = '2';
         traceback(alignment1_tmp, alignment2_tmp, index_x, index_y, current_matrix, alignment_index, matrix_index, highest_score);
      }
   }

   // no alignment is made
   // because there are too few overlaps between them
   if (alignment1_vector.size() == 0) {
      num_total_bases_percent_similarity   = 0;
      num_matched_bases_percent_similarity = 0;

      alignment_best = "NO ALIGNMENT\n";
   }
   else {
      // find indels at the ends of alignments
      std::regex rx_5_prime_gap("^-+");
      std::regex rx_3_prime_gap("-+$");

      std::smatch smatch_5_prime_insertion;
      std::smatch smatch_5_prime_deletion;
      std::smatch smatch_3_prime_insertion;
      std::smatch smatch_3_prime_deletion;

      int num_insertions_5_prime;
      int num_insertions_3_prime;
      int num_deletions_5_prime;
      int num_deletions_3_prime;

      // 5'-end insertions
      // -----AAA
      // AAAAAAAA
      if (std::regex_search(alignment1_vector[0], smatch_5_prime_insertion, rx_5_prime_gap)) {
         num_insertions_5_prime = smatch_5_prime_insertion[0].length();
      }
      else {
         num_insertions_5_prime = 0;
      }

      // 5'-end deletions
      // AAAAAAAA
      // -----AAA
      if (std::regex_search(alignment2_vector[0], smatch_5_prime_deletion, rx_5_prime_gap)) {
         num_deletions_5_prime = smatch_5_prime_deletion[0].length();
      }
      else {
         num_deletions_5_prime = 0;
      }

      // 3'-end insertions
      // AAA-----
      // AAAAAAAA
      if (std::regex_search(alignment1_vector[0], smatch_3_prime_insertion, rx_3_prime_gap)) {
         num_insertions_3_prime = smatch_3_prime_insertion[0].length();
      }
      else {
         num_insertions_3_prime = 0;
      }

      // 3'-end deletions
      // AAAAAAAA
      // AAA-----
      if (std::regex_search(alignment2_vector[0], smatch_3_prime_deletion, rx_3_prime_gap)) {
         num_deletions_3_prime = smatch_3_prime_deletion[0].length();
      }
      else {
         num_deletions_3_prime = 0;
      }

      //----------------------------------------------------------------------
      // calculate percent similarity
      //----------------------------------------------------------------------
      num_total_bases_percent_similarity   = 0;
      num_matched_bases_percent_similarity = 0;

      //
      // main region
      //
      int start_index;
      int end_index;

      start_index = std::max(num_insertions_5_prime, num_deletions_5_prime);
      end_index   = alignment1_vector[0].length() - std::max(num_insertions_3_prime, num_deletions_3_prime) - 1;

      for (int it_alignment = start_index; it_alignment <= end_index; it_alignment++) {
         num_total_bases_percent_similarity++;

         if (alignment1_vector[0][it_alignment] == alignment2_vector[0][it_alignment]) {
            num_matched_bases_percent_similarity++;
         }
      }

      //
      // insertions at the 3'-end
      //
      // AAA---
      // AAAAAA
      // use the first three bases of outer_3_end
      if (num_insertions_3_prime > 0) {
         for (int it_alignment = (alignment1_vector[0].length() - num_insertions_3_prime); it_alignment < alignment1_vector[0].length(); it_alignment++) {
            num_total_bases_percent_similarity++;

            if (outer_3_end[it_alignment] == alignment2_vector[0][it_alignment]) {
               num_matched_bases_percent_similarity++;
            }
         }
      }

      //
      // insertions at the 5'-end
      //
      // ---AAA
      // AAAAAA
      // use the last three bases of outer_5_end
      if (num_insertions_5_prime > 0) {
         for (int it_alignment = 0; it_alignment < num_insertions_5_prime; it_alignment++) {
            num_total_bases_percent_similarity++;

            if (outer_5_end[outer_5_end.length() - num_insertions_5_prime + it_alignment] == alignment2_vector[0][it_alignment]) {
               num_matched_bases_percent_similarity++;
            }
         }
      }

      alignment_best = alignment1_vector[0] + "\n" + alignment2_vector[0] + "\n";
   }

   delete[] match_matrix;
   delete[] gap_1_matrix;
   delete[] gap_2_matrix;
}
