#ifndef HASHTABLE_H
#define HASHTABLE_H

#include <string>
#include <vector>
#include <cstddef>
#include <stdexcept>
#include <functional>

class HashTable {
public:
    using size_type = std::size_t;

    explicit HashTable(size_type initial_capacity = 16);
    ~HashTable();

    void insert(const std::string& key, int value);
    
    bool contains(const std::string& key) const noexcept;
    
    int get(const std::string& key) const;
    
    int get_or(const std::string& key, int default_value) const noexcept;
    
    bool erase(const std::string& key) noexcept;
    
    void clear() noexcept;

    size_type size() const noexcept { return size_; }
    
    size_type capacity() const noexcept { return capacity_; }
    
    float load_factor() const noexcept;

private:
    struct Entry {
        std::string key;
        int value;
        bool occupied;

        Entry() : value(0), occupied(false) {}
    };

    std::vector<Entry> table_;
    size_type size_;
    size_type capacity_;

    static size_type hash_string(const std::string& key) noexcept;
    
    size_type find_position(const std::string& key) const noexcept;
    
    size_type find_insert_position(const std::string& key) noexcept;
    
    void resize(size_type new_capacity);
    
    bool should_resize() const noexcept;
};

#endif // HASHTABLE_H
