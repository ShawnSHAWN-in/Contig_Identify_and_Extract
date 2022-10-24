#include"Kseq_Cpp.h"
#include<vector>
#include<utility>
#include<regex>
#include<ranges>
#include<thread>//获取最大线程数
#include<algorithm>//获取最大值位置
#include <filesystem> // use std::filesystem::path
#include<cctype>
#include <fstream>
#include <numeric>
#include <cstring>
#include <seqan3/io/sequence_file/all.hpp>
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


    seqan3::regex_validator type_Validator{"(AA|NT)"};
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

    seqan3::argument_parser myparser{"IEContig", argc, argv}; // initialise myparser

    //Self defined options or flags
    cmd_arguments args{};

    // ... add information, options, flags and positional options
    initialise_argument_parser(myparser, args);
    try
    {
        myparser.parse(); // trigger command line parsing
    }
    catch (seqan3::argument_parser_error const & ext) // catch user errors
    {
        seqan3::debug_stream << "[Winter has come] " << ext.what() << "\n"; // customise your error message
        //return -1;
    }
    
    using types_nt = seqan3::type_list< std::string,std::vector<seqan3::dna15>>;
    using fields_nt = seqan3::fields<seqan3::field::id,seqan3::field::seq>;
    using sequence_nt_record_type = seqan3::sequence_record<types_nt, fields_nt>;
    using seq_nt_set_t=std::vector<sequence_nt_record_type>;
    seq_nt_set_t seq_nt_set;
    //输入的Contig核苷酸序列
    Kseq_Cpp file_input(args.file_nucleotide.generic_string(),16384);
    while (file_input.get_seq())
    {
        std::regex tmp_reg1("-");
        std::regex tmp_reg2("[^TVGHCDMKNYSABWR]");
        std::string tmp_string1=std::regex_replace(file_input.seq->seq.s,tmp_reg1,"");
        std::string tmp_string2=std::regex_replace(tmp_string1.c_str(),tmp_reg2,"N");
        //std::cout<<tmp_string2<<"\n";
        //用于输出比对上的氨基酸序列
        sequence_nt_record_type record{
                std::string(file_input.seq->name.s),
                make_vector(tmp_string2|seqan3::views::char_to<seqan3::dna15> )
                };
        seq_nt_set.push_back(record);
    }


    using types_aa = seqan3::type_list< std::string,std::vector<seqan3::aa27>>;
    using fields_aa = seqan3::fields<seqan3::field::id,seqan3::field::seq>;
    using sequence_aa_record_type = seqan3::sequence_record<types_aa, fields_aa>;
    using seq_aa_set_t=std::vector<sequence_aa_record_type>;
    //Conserved AA file input
    Kseq_Cpp file_conserved(args.file_amino.generic_string(),16384);
    seq_aa_set_t cd_aa_set;
    while (file_conserved.get_seq())
    {
        //printf("qual: %s\n", file_conserved.seq->seq.s);
        std::regex tmp_reg1("[-+.]");
        std::regex tmp_reg2("[^ARNDCQEGHILKMFPSTWYVUOBJZX*]");
        std::string tmp_string1=std::regex_replace(file_conserved.seq->seq.s,tmp_reg1,"");
        std::string tmp_string2=std::regex_replace(tmp_string1.c_str(),tmp_reg2,"X");
        //std::cout<<tmp_string2<<"\n";
        sequence_aa_record_type record{
                std::string(file_input.seq->name.s),
                make_vector(tmp_string2|seqan3::views::char_to<seqan3::aa27>)
                };
        cd_aa_set.push_back(record);
    }

    auto output_config = 
        seqan3::align_cfg::output_score{} | 
        seqan3::align_cfg::output_begin_position{}|
        seqan3::align_cfg::output_end_position{} | 
        seqan3::align_cfg::output_alignment{};

    auto config =
        seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},//long local
                                         seqan3::align_cfg::free_end_gaps_sequence2_leading{false},//short global
                                         seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                         seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}}
        | seqan3::align_cfg::scoring_scheme{seqan3::aminoacid_scoring_scheme{seqan3::aminoacid_similarity_matrix::blosum62}}
        | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-11}, seqan3::align_cfg::extension_score{-1}}
        | seqan3::align_cfg::parallel{15}
        | output_config;

    auto fasta_file = std::filesystem::current_path() / "result.fasta";
    //seqan3::sequence_file_output fout{std::cout,seqan3::format_fasta{}};
    std::cout<<args.file_output.generic_string()<<"\n";
    seqan3::sequence_file_output fout{args.file_output.generic_string(),seqan3::fields< seqan3::field::id,seqan3::field::seq>{}};
    
    for(size_t i = 0 ; i < seq_nt_set.size(); i++){
        seq_aa_set_t tmp_translation(6);
        tmp_translation[0]=sequence_aa_record_type{
                seq_nt_set[i].id(),
                make_vector(seq_nt_set[i].sequence()|seqan3::views::translate_single(seqan3::translation_frames::forward_frame0))
                };
        tmp_translation[1]=sequence_aa_record_type{
                seq_nt_set[i].id(),
                make_vector(seq_nt_set[i].sequence()|seqan3::views::translate_single(seqan3::translation_frames::forward_frame1))
                };                
        tmp_translation[2]=sequence_aa_record_type{
                seq_nt_set[i].id(),
                make_vector(seq_nt_set[i].sequence()|seqan3::views::translate_single(seqan3::translation_frames::forward_frame2))
                };
        tmp_translation[3]=sequence_aa_record_type{
                seq_nt_set[i].id(),
                make_vector(seq_nt_set[i].sequence()|seqan3::views::complement | std::views::reverse|seqan3::views::translate_single(seqan3::translation_frames::forward_frame0))
                };
        tmp_translation[4]=sequence_aa_record_type{
                seq_nt_set[i].id(),
                make_vector(seq_nt_set[i].sequence()|seqan3::views::complement | std::views::reverse|seqan3::views::translate_single(seqan3::translation_frames::forward_frame1))
                };
        tmp_translation[5]=sequence_aa_record_type{
                seq_nt_set[i].id(),
                make_vector(seq_nt_set[i].sequence()|seqan3::views::complement | std::views::reverse|seqan3::views::translate_single(seqan3::translation_frames::forward_frame2))
                };

        std::vector<double> tmp_score;
        std::vector<std::pair<int,int>> tmp_positions;
        std::vector<std::pair<seqan3::aa27_vector, seqan3::aa27_vector>> tmp_sequence;
        for (auto &tmp_1 : tmp_translation)
        {
            for (auto &tmp_2 : cd_aa_set)
            {
                tmp_sequence.push_back(
                    std::make_pair(
                        tmp_1.sequence(),//Out_put sequence,first
                        tmp_2.sequence()));//Conserved aa sequence,second  
            }
        }
        for (auto const &res : seqan3::align_pairwise(tmp_sequence, config))
        {
            tmp_score.push_back((double)res.score()/(double)std::get<1>(res.alignment()).size());
            tmp_positions.push_back(std::make_pair(res.sequence1_begin_position(),res.sequence1_end_position()));
        }  
        auto maxPosition=std::max_element(tmp_score.begin(),tmp_score.end());
        //std::cout<<"Threshold:"<<args.self_set_threshold<<"\n";
        //std::cout<<"Output_type:"<<args.output_type<<"\n";
        //std::cout<<"MaxScore:"<<*maxPosition<<"\n";
        if (*maxPosition > args.self_set_threshold)
        {
            std::cout<<"Scores:"<<*maxPosition<<"\n";
            auto index = std::distance(tmp_score.begin(), maxPosition);
            // seqan3::debug_stream<<std::vector(beg,end)<<"\n";
            if (args.output_type=="AA")
            {                
                auto tmp_seq = tmp_sequence[index].first;
                auto beg = tmp_seq.begin() + tmp_positions[index].first;  // Begin position
                auto end = tmp_seq.begin() + tmp_positions[index].second; // End position
                sequence_aa_record_type record{
                    seq_nt_set[i].id(),
                    std::vector(beg, end)};
                fout.push_back(record);
            }else if (args.output_type=="NT")
            {
                //核苷酸序列扣取范围
                auto start_position=tmp_positions[index].first * 3 -2;  // Begin position
                auto stop_position=tmp_positions[index].second*3; // End position
                if (index/3==0)
                {
                    auto out_seq=make_vector(seq_nt_set[i].sequence()|std::views::drop(index%3));
                    sequence_nt_record_type record{
                        seq_nt_set[i].id(),
                        std::vector(out_seq.begin()+start_position,out_seq.begin()+stop_position)
                    };
                    fout.push_back(record);                                            
                }else if (index/3==1)
                {
                    auto out_seq=make_vector(seq_nt_set[i].sequence()|seqan3::views::complement | std::views::reverse|std::views::drop(index%3));
                    sequence_nt_record_type record{
                        seq_nt_set[i].id(),
                        std::vector(out_seq.begin()+start_position,out_seq.begin()+stop_position)
                    };                      
                    fout.push_back(record);               
                }                
            }
        }
    }
}
