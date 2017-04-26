/* module name */
%module evaluate
%{
/* headers declarations */
extern int end_index;
extern int read_length;
extern int start_index;
extern int gap_extension_penalty;
extern int gap_opening_penalty;
extern int match_gain;
extern int mismatch_penalty;
extern int max_read_length;
extern int num_yyns_substitution_local_best;
extern int num_ynys_substitution_local_best;
extern int num_nyys_substitution_local_best;
extern int num_nyns_substitution_local_best;
extern int num_nnns_substitution_local_best;
extern int num_yyns_insertion_local_best;
extern int num_nyys_insertion_local_best;
extern int num_nyns_insertion_local_best;
extern int num_nnns_insertion_local_best;
extern int num_yyns_deletion_local_best;
extern int num_nyys_deletion_local_best;
extern int num_nyns_deletion_local_best;
extern int num_nnns_deletion_local_best;
extern int num_from_substitution_to_deletion_local_best;
extern int num_nyys_substitution_trim_local_best;
extern int num_nyys_insertion_trim_local_best;
extern int num_nyys_deletion_trim_local_best;
extern int num_not_evaluated_substitution;
extern int num_not_evaluated_insertion;
extern int num_not_evaluated_deletion;
extern int num_total_bases_percent_similarity;
extern int num_matched_bases_percent_similarity;
extern int ref_seq_index;

extern unsigned int max_candidates;

extern std::string outer_3_end;
extern std::string outer_5_end;
extern std::string read_name;
extern std::string string1;
extern std::string string2;
extern std::string substitutions;
extern std::string insertions;
extern std::string deletions;
extern std::string alignment_best;
extern std::string error_index_best;
extern std::string random_alignment1;
extern std::string random_alignment2;
extern std::string strand;

extern bool is_detail;
extern bool is_trimmed;
extern bool no_end_gap_penalty;
extern bool too_many_candidates;

void initialize_variables();
void decode_errors();
void debug_print_variables();
void fill_matrixes();
void find_best_alignment(int* position_vector_local_best, int* corrected_position_vector_local_best);
void print_matrixes();
void give_random_alignment();
void calculate_percent_similarity();
%}

%include std_string.i
%include carrays.i

%array_functions(int, intp)

int end_index;
int read_length;
int start_index;
int gap_extension_penalty;
int gap_opening_penalty;
int match_gain;
int mismatch_penalty;
int max_read_length;
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
int num_not_evaluated_substitution;
int num_not_evaluated_insertion;
int num_not_evaluated_deletion;
int num_total_bases_percent_similarity;
int num_matched_bases_percent_similarity;
int ref_seq_index;

unsigned int max_candidates;

std::string outer_3_end;
std::string outer_5_end;
std::string read_name;
std::string string1;
std::string string2;
std::string substitutions;
std::string insertions;
std::string deletions;
std::string alignment_best;
std::string error_index_best;
std::string random_alignment1;
std::string random_alignment2;
std::string strand;

bool is_detail;
bool is_trimmed;
bool no_end_gap_penalty;
bool too_many_candidates;

void initialize_variables();
void decode_errors();
void debug_print_variables();
void fill_matrixes();
void find_best_alignment(int* position_vector_local_best, int* corrected_position_vector_local_best);
void give_random_alignment();
void print_matrixes();
void calculate_percent_similarity();
