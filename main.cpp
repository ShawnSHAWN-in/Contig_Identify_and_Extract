#include"Kseq_Cpp.h"
#include<vector>
#include<utility>
#include<regex>
#include<ranges>
#include <filesystem> // use std::filesystem::path
#include <fstream>
#include <numeric>
#include <cstring>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include<seqan3/alphabet/nucleotide/dna4.hpp>
#include<seqan3/alphabet/nucleotide/dna15.hpp>
#include<seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/views/all.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/alphabet/views/translate.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>


#include<iostream>

struct cmd_arguments
{
    std::filesystem::path file_nucleotide{};
    std::string input_type{};
    std::filesystem::path file_amino{};
    std::filesystem::path file_output{};
    std::string output_type{};
    int self_set_threshold=0;
};

void initialise_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{

    parser.info.version = "1.0.0";
    parser.info.app_name="IEContig";
    parser.info.author = "CuiXinyuan";
    parser.info.citation="None any Citations";
    parser.info.date="2022年10月9日";
    parser.info.email="145769685@qq.com";
    //parser.info.examples="IEContig -i <INPUT_CONTIG.fasta> -t <TYPE_OF_INPUT_CONTIG> -I <reference.fasta> -T <AA> -o <RESULT_FILE.fasta> -y <OUT_TYPE> -s <SELF_SET_THRESHOLD>";
    parser.info.short_description = "This Script is a test for arguments!";
    //parser.info.description="The program extract specify region(gene) pre-determined by user ,adepted both to NT and AA sequence.\nNote: input sequence necessiate collect homologous sequence";

    //ADD OPTION
    parser.add_option(args.file_nucleotide, 'i', "input_contig", "This is term of input contig sequence!",
                      seqan3::option_spec::required,seqan3::input_file_validator{{"fa","fas","fna","ffn","ffa","frn","fasta"}});


    seqan3::regex_validator type_Validator{"(AA|NT|aa|nt|Aa|aA|Nt|nT)"};
    parser.add_option(args.input_type,'t',"type_of_input_contig","Input sequence type [NT|AA] ",
                      seqan3::option_spec::required,type_Validator);

    parser.add_option(args.file_amino, 'I', "homologous_region", "collected homologous sequence (genes) from (NCBI or others database)",
                      seqan3::option_spec::required,seqan3::input_file_validator{{"fa","fas","fna","ffn","ffa","frn","fasta"}});

    parser.add_option(args.file_output, 'o', "result_file", "Output file name",
                      seqan3::option_spec::required);

    parser.add_option(args.output_type,'y',"out_type","Output sequence type [NT|AA] ",
                      seqan3::option_spec::required,type_Validator);

    seqan3::arithmetic_range_validator self_set_threshold_validator{-5, 6};
    parser.add_option(args.self_set_threshold, 's', "self_set_threshold", "Ourseleves defined threshhold to filter out the special sequence",
                     seqan3::option_spec::required,self_set_threshold_validator);
}

template <std::ranges::view T>
auto make_vector(T a ){
    return std::vector(a.begin(),a.end());
}


int main(int argc, char ** argv)
{
    using namespace seqan3::literals;

//    seqan3::argument_parser myparser{"IEContig", argc, argv}; // initialise myparser
//
//    //Self defined options or flags
//    cmd_arguments args{};
//
//    // ... add information, options, flags and positional options
//    initialise_argument_parser(myparser, args);
//    try
//    {
//        myparser.parse(); // trigger command line parsing
//    }
//    catch (seqan3::argument_parser_error const & ext) // catch user errors
//    {
//        seqan3::debug_stream << "[Winter has come] " << ext.what() << "\n"; // customise your error message
//        //return -1;
//    }
//    std::cout<<"Thanks";

    using seq_na_t=std::pair<std::string,seqan3::dna15_vector>;
    using seq_aa_t=std::pair<std::string,seqan3::aa27_vector>;
    using seq_na_set_t=std::vector<seq_na_t>;
    using seq_aa_set_t=std::vector<seq_aa_t>;

    seq_na_set_t seq_na_set;
    seq_aa_set_t seq_aa_set;

    Kseq_Cpp file_input(std::string{"My_fasta.fasta"},16384);
    while (file_input.get_seq())
    {
        std::regex tmp_reg1("-");
        std::regex tmp_reg2("[^TVGHCDMKNYSABWR]");
        std::string tmp_string1=std::regex_replace(file_input.seq->seq.s,tmp_reg1,"");
        std::string tmp_string2=std::regex_replace(tmp_string1.c_str(),tmp_reg2,"N");
        //std::cout<<tmp_string2<<"\n";
        //正向错位翻译
        seq_aa_set.push_back(
            std::make_pair(std::string(file_input.seq->name.s),
            make_vector(tmp_string2|seqan3::views::char_to<seqan3::dna15>| seqan3::views::translate_single(seqan3::translation_frames::forward_frame0))
            )
        );
        seq_aa_set.push_back(
            std::make_pair(std::string(file_input.seq->name.s),
            make_vector(tmp_string2|seqan3::views::char_to<seqan3::dna15>| seqan3::views::translate_single(seqan3::translation_frames::forward_frame1))
            )
        );
        seq_aa_set.push_back(
            std::make_pair(std::string(file_input.seq->name.s),
            make_vector(tmp_string2|seqan3::views::char_to<seqan3::dna15>| seqan3::views::translate_single(seqan3::translation_frames::forward_frame2))
            )
        );
        //反向互补错位翻译
        seq_aa_set.push_back(
            std::make_pair(std::string(file_input.seq->name.s),
            make_vector(tmp_string2|seqan3::views::char_to<seqan3::dna15>| seqan3::views::complement | std::views::reverse| seqan3::views::translate_single(seqan3::translation_frames::forward_frame0))
            )
        );       
        seq_aa_set.push_back(
            std::make_pair(std::string(file_input.seq->name.s),
            make_vector(tmp_string2|seqan3::views::char_to<seqan3::dna15>| seqan3::views::complement |std::views::reverse|seqan3::views::translate_single(seqan3::translation_frames::forward_frame1))
            )
        );    
        seq_aa_set.push_back(
            std::make_pair(std::string(file_input.seq->name.s),
            make_vector(tmp_string2|seqan3::views::char_to<seqan3::dna15>| seqan3::views::complement |std::views::reverse|seqan3::views::translate_single(seqan3::translation_frames::forward_frame2))
            )
        );             
     
    }
    //seqan3::debug_stream<< seq_aa_set[0]<<"\n";
    
    Kseq_Cpp file_conserved(std::string{"cd21530.FASTA"},16384);
    seq_aa_set_t cd_aa_set;
    while (file_conserved.get_seq())
    {
        //printf("qual: %s\n", file_conserved.seq->seq.s);
        std::regex tmp_reg1("[-+.]");
        std::regex tmp_reg2("[^ARNDCQEGHILKMFPSTWYVUOBJZX*]");
        std::string tmp_string1=std::regex_replace(file_conserved.seq->seq.s,tmp_reg1,"");
        std::string tmp_string2=std::regex_replace(tmp_string1.c_str(),tmp_reg2,"X");
        //std::cout<<tmp_string2<<"\n";
        cd_aa_set.push_back(
            std::make_pair(
                std::string(file_input.seq->name.s),
                make_vector(tmp_string2|seqan3::views::char_to<seqan3::aa27>)
                )
        );
    }

    //seqan3::debug_stream<< cd_aa_set[0]<<"\n";

    std::vector<int8_t> score(cd_aa_set.size());
    using sequence_pair_t = std::pair<seqan3::aa27_vector, seqan3::aa27_vector>;
    using name_pair_t = std::pair<std::string,std::string>;

    std::pair<std::vector<name_pair_t>,std::vector<sequence_pair_t>> sequences(
        std::vector<name_pair_t>(cd_aa_set.size()*seq_aa_set.size()),
        std::vector<sequence_pair_t>(cd_aa_set.size()*seq_aa_set.size())
    );

    std::cout<<"cd_aa_set:"<<cd_aa_set.size()<<"\n";
    std::cout<<"seq_aa_set:"<<seq_aa_set.size()<<"\n";

    for(size_t i = 0 ; i < cd_aa_set.size(); i++){
        for(size_t n = 0 ; n< seq_aa_set.size() ; n++ ){
            std::cout<<"Cycle:"<<i<<"\t"<<"Colomn:"<<n<<"\n";
            sequences.first[i*seq_aa_set.size() + n]={
                std::ref(cd_aa_set[i].first),
                std::ref(seq_aa_set[n].first)
            };
            sequences.second[i*seq_aa_set.size() + n]={
                std::ref(cd_aa_set[i].second),
                std::ref(seq_aa_set[n].second)
            };        
        }
    }

    auto output_config = 
        seqan3::align_cfg::output_score{} | 
        seqan3::align_cfg::output_begin_position{}|
        seqan3::align_cfg::output_end_position{} | 
        seqan3::align_cfg::output_alignment{};

    auto config =
        seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                         seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
                                         seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                         seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}}
        | seqan3::align_cfg::scoring_scheme{seqan3::aminoacid_scoring_scheme{seqan3::aminoacid_similarity_matrix::blosum62}}
        | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-11}, seqan3::align_cfg::extension_score{-1}}
        | seqan3::align_cfg::parallel{15}
        | output_config;

    std::cout<<"END"<<"\n";
    for (auto const & res : seqan3::align_pairwise(sequences.second, config))
        seqan3::debug_stream << "Score: " << res.score() << '\n';

    std::cout<<"END2"<<"\n";
}