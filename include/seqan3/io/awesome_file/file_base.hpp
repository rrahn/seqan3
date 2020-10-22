#pragma once

#include <seqan3/std/concepts>
#include <seqan3/std/filesystem>
#include <fstream>
#include <seqan3/std/iterator>
#include <memory>
#include <seqan3/std/ranges>
#include <seqan3/std/type_traits>
#include <vector>

#include <seqan3/io/stream/iterator.hpp>
#include <seqan3/io/awesome_file/format_base.hpp>

namespace seqan3::awesome
{


// TODOs:
// Refactor the record classes using policies
// Let the wrapper interface be automatically generated by the compiler -> like get the fields the type list.
//

template <typename record_base_t>
class input_file
{
private:

    //!\brief The polymorphic base format type.
    using base_format_type = format_base<record_base_t>;
    using base_format_ptr_type = std::unique_ptr<base_format_type>;

    std::vector<base_format_ptr_type> registered_formats{};
    base_format_type * selected_format{};

    // We actually want to use a stream here.
    std::ifstream input_fstream{}; // We need to have the file stream and maybe something else.
    std::basic_istream<char> * input_stream{nullptr};
public:

    struct iterator;

    input_file() = default;
    input_file(input_file const &) = delete;
    input_file(input_file && other) noexcept : input_file{}
    {
        swap(other);
    }

    input_file & operator=(input_file const &) = delete;
    input_file & operator=(input_file && other) noexcept
    {
        swap(other);
        return *this;
    }

    ~input_file() = default;

    // TODO: Only after rebind we can check the property.
    template <typename ...formats_t>
        requires (std::derived_from<typename std::remove_reference_t<formats_t>::rebind_record<record_base_t>, base_format_type> && ...)
    input_file(std::filesystem::path file_path, formats_t && ...formats)
    { // Register all given formats.
        (registered_formats.emplace_back(
                std::make_unique<typename std::remove_reference_t<formats_t>::rebind_record<record_base_t>>(
                        std::forward<formats_t>(formats))), ...);

        // Read the extension and guess the format.
        // Then set the format.
        std::string file_ext = file_path.extension().string().substr(1);

        for (auto & format : registered_formats)
        {
            if (std::ranges::any_of(format->extensions(), [&] (std::string ext) { return ext == file_ext; }))
            {
                selected_format = format.get();
                break;
            }
        }

        // Next step: open file_handle from the path.
        // Also add the decompression format.
        input_fstream.open(file_path.c_str());

        if (!input_fstream.good())
            throw std::runtime_error{"Could not open the file stream."};

        input_stream = &input_fstream; // Set the input stream to the file stream.
    }

    template <typename format_t>
        requires (std::derived_from<typename std::remove_reference_t<format_t>::rebind_record<record_base_t>, base_format_type>)
    input_file(std::istream & outer_input_stream, format_t && format)
    { // Register the format.
        registered_formats.emplace_back(
                std::make_unique<typename std::remove_reference_t<format_t>::rebind_record<record_base_t>>(
                        std::forward<format_t>(format)));

        selected_format = registered_formats.front().get();

        // Next step: open file_handle from the path.
        // Also add the decompression format.
        if (!outer_input_stream.good())
            throw std::runtime_error{"Could not open the file stream."};

        // Set the streambuf_iterator.
        input_stream = &outer_input_stream;
    }

    iterator begin()
    {
        return iterator{this};
    }

    iterator begin() const = delete;


    std::default_sentinel_t end()
    {
        return std::default_sentinel;
    }

    std::default_sentinel_t end() const = delete;

    void swap(input_file & rhs) noexcept
    {
        using std::swap;

        selected_format = std::addressof(*rhs.selected_format);
        swap(registered_formats, rhs.registered_formats);
        swap(input_fstream, rhs.input_fstream);
    }
};

template <typename record_base_t>
class input_file<record_base_t>::iterator
{
private:

    input_file * host_file{};
    record_base_t * cached_record{nullptr};
    // The idea is to present the base type here.
    // Could we remedy this by making the base a CRTP?
    bool eof{false};

public:

    using value_type = record_base_t;
    using reference = std::add_lvalue_reference_t<value_type>;
    using pointer = record_base_t *;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::input_iterator_tag;
    using iterator_concept = std::input_iterator_tag;

    iterator() = default;
    iterator(iterator const &) = default;
    iterator(iterator &&) = default;
    iterator & operator=(iterator const &) = default;
    iterator & operator=(iterator &&) = default;
    ~iterator() = default;

    iterator(input_file * host_file) : host_file{host_file}
    {
        ++(*this);
    }

    reference operator*() const
    {
        // This can't be a header and a record at the same time? Can it?
        // It could by making it a variant, but do we want to?
        return *cached_record;
    }

    pointer operator->() const
    {
        // This can't be a header and a record at the same time? Can it?
        // It could by making it a variant, but do we want to?
        return cached_record;
    }

    iterator & operator++()
    {
        cached_record = host_file->selected_format->read_record(*(host_file->input_stream));
        return *this;
    }

    void operator++(int)
    {
        ++*this;
    }

    bool operator==(std::default_sentinel_t const &) const noexcept
    {
        return cached_record == nullptr;
    }

    // friend bool operator==(std::default_sentinel_t const & lhs, iterator const & rhs) noexcept
    // {
    //     return rhs == lhs;
    // }

    // bool operator!=(std::default_sentinel_t const & rhs) const noexcept
    // {
    //     return !(*this == rhs);
    // }

    // friend bool operator!=(std::default_sentinel_t const & lhs, iterator const & rhs) noexcept
    // {
    //     return rhs != lhs;
    // }

    // Do we need this???
    // friend constexpr range_rvalue_reference_t<V> iter_move(const iterator& i)
    //   noexcept(noexcept(ranges::iter_move(i.current_)));
    // friend constexpr void iter_swap(const iterator& x, const iterator& y)
    //   noexcept(noexcept(ranges::iter_swap(x.current_, y.current_)))
    //   requires indirectly_­swappable<iterator_t<V>>;
};

// A fixed instance of the base class template.
// using sequence_file = file<sequence_file_format_base>;

} // namespace seqan3::awesome
