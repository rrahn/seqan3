
#include <chrono>
#include <vector>
#include <seqan3/std/algorithm>

#include <cereal/archives/binary.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/views/chunk.hpp>
#include <seqan3/search/dream_index/technical_binning_directory.hpp>

auto load_reads(std::string_view const reads_file)
{
    struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
    {
        using sequence_alphabet = seqan3::dna4;
        using sequence_legal_alphabet = seqan3::dna4;
    };

    std::vector<seqan3::dna4_vector> reads{};
    {
        seqan3::sequence_file_input<dna4_traits> reads_in{reads_file};
        for (auto & read : reads_in)
            reads.push_back(seqan3::get<seqan3::field::seq>(read));
    }

    return reads;
}

auto load_binning_directory(std::string_view genome_ibf_path)
{
    constexpr size_t const bin_count{64};
    constexpr size_t const hash_num{2u};
    constexpr size_t const ibf_size{68'719'476'736}; // 8 GiB
    constexpr size_t bin_size{ibf_size / bin_count};

    seqan3::ibf_config const cfg{.number_of_bins = seqan3::bin_count{bin_count},
                                 .size_of_bin = seqan3::bin_size{bin_size},
                                 .number_of_hash_functions = seqan3::hash_function_count{hash_num},
                                 .threads = 16};

    auto hash_adaptor = seqan3::views::kmer_hash(seqan3::ungapped{19u});

    using tbd_t = seqan3::technical_binning_directory<seqan3::data_layout::uncompressed,
                                                      decltype(hash_adaptor),
                                                      seqan3::dna4>;

    tbd_t tbd{cfg, hash_adaptor};

    {   // serialise the binning dir
        std::ifstream binning_dir_stream{genome_ibf_path.data()};
        cereal::BinaryInputArchive iarch{binning_dir_stream};
        iarch(tbd);
    }

    // seqan3::technical_binning_directory tbd{genome | ranges::views::chunk(chunk_size),
    //                                         ),
    //                                         cfg};

    return tbd;
}

template <typename read_t, typename t>
auto count_reads(std::vector<read_t> && reads, t && binning_dir, size_t const threshold = 5)
{
    std::vector<std::pair<size_t, size_t>> count_result{};
    count_result.reserve(reads.size());

    auto agent = binning_dir.counting_agent();

    size_t query_id{};
    std::ranges::for_each(reads, [&] (auto const & query)
    {
        size_t hit_count{};
        std::ranges::for_each(agent.count_query(query), [&] (auto const count)
        {
            hit_count += count > threshold;
        });
        count_result.push_back(std::pair{query_id++, hit_count});
    });
    return count_result;
}

int main(int const argc, char const * argv[])
{
    if (argc < 3)
        throw std::invalid_argument{"Call: ./prog genome.ibf reads.fa"};

    try {

        auto time_point = std::chrono::high_resolution_clock::now();
        auto reads = ::load_reads(argv[2]);
        auto load_reads_duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - time_point);

        time_point = std::chrono::high_resolution_clock::now();
        auto binning_dir = ::load_binning_directory(argv[1]);
        auto load_tbd_duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - time_point);

        time_point = std::chrono::high_resolution_clock::now();
        auto count_result = ::count_reads(std::move(reads), std::move(binning_dir));
        auto count_reads_duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - time_point);

        std::cout << "Timings:\n"
                  << "Load reads: " << load_reads_duration.count() << "ms\n"
                  << "Load tbdir: " << load_tbd_duration.count() << "ms\n"
                  << "Cnt  reads: " << count_reads_duration.count() << "ms\n";


        std::ofstream result_out{argv[3]};
        std::ranges::for_each(count_result, [&] (auto const & result)
        {
            auto const & [id, count] = result;
            result_out << "query_" << id << "," << count << '\n';
        });

        std::cout << "Done\n";
        // std::cout << "genome: size = " << genome.size() << " reads count = " << reads.size() << "\n";

        // seqan3::ibf_config const cfg{seqan3::bin_count{bin_count},
        //                             seqan3::bin_size{bin_size},
        //                             seqan3::hash_function_count{hash_num}};

        // seqan3::technical_binning_directory tbd{genome | ranges::views::chunk(chunk_size),
        //                                         seqan3::views::kmer_hash(seqan3::ungapped{19u}),
        //                                         cfg};

        // auto agent = tbd.counting_agent();
        // for (auto _ : state)
        //     for (auto && query : reads)
        //         benchmark::DoNotOptimize(agent.count_query(query));

    } catch (...) {
        std::throw_with_nested(std::runtime_error("Something went wrong!") );
    }
}
