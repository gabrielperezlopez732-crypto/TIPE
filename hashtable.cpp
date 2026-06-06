#include "hashtable.h"
#include <cstring>
#include <algorithm>

HashTable::HashTable(size_type initial_capacity) 
    : size_(0), capacity_(initial_capacity) {
    if (capacity_ < 16) capacity_ = 16;
    table_.resize(capacity_);
}

HashTable::~HashTable() {
    clear();
}

HashTable::size_type HashTable::hash_string(const std::string& key) noexcept {
    size_type hash = 14695981039346656037ULL;
    for (char c : key) {
        hash ^= static_cast<unsigned char>(c);
        hash *= 1099511628211ULL; // FNV prime
    }
    return hash;
}

HashTable::size_type HashTable::find_position(const std::string& key) const noexcept {
    size_type hash = hash_string(key);
    size_type index = hash % capacity_;
    size_type steps = 0;

    while (steps < capacity_) {
        if (!table_[index].occupied) {
            return capacity_; // Not found
        }
        if (table_[index].key == key) {
            return index; // Found
        }
        index = (index + 1) % capacity_;
        steps++;
    }
    return capacity_; // Not found
}

HashTable::size_type HashTable::find_insert_position(const std::string& key) noexcept {
    size_type hash = hash_string(key);
    size_type index = hash % capacity_;
    size_type steps = 0;
    size_type first_tombstone = capacity_;

    while (steps < capacity_) {
        if (!table_[index].occupied) {
            return (first_tombstone != capacity_) ? first_tombstone : index;
        }
        if (table_[index].key == key) {
            return index; 
        }
        index = (index + 1) % capacity_;
        steps++;
    }
    return capacity_;
}

void HashTable::resize(size_type new_capacity) {
    if (new_capacity <= capacity_) return;

    std::vector<Entry> old_table = std::move(table_);
    capacity_ = new_capacity;
    table_.clear();
    table_.resize(capacity_);
    size_ = 0;

    for (const auto& entry : old_table) {
        if (entry.occupied) {
            insert(entry.key, entry.value);
        }
    }
}

bool HashTable::should_resize() const noexcept {
    return size_ > (capacity_ * 3) / 4;
}

void HashTable::insert(const std::string& key, int value) {
    if (should_resize()) {
        resize(capacity_ * 2);
    }

    size_type pos = find_insert_position(key);
    if (pos < capacity_) {
        if (!table_[pos].occupied) {
            size_++;
        }
        table_[pos].key = key;
        table_[pos].value = value;
        table_[pos].occupied = true;
    }
}

bool HashTable::contains(const std::string& key) const noexcept {
    return find_position(key) < capacity_;
}

int HashTable::get(const std::string& key) const {
    size_type pos = find_position(key);
    if (pos < capacity_) {
        return table_[pos].value;
    }
    throw std::out_of_range("Key not found in HashTable");
}

int HashTable::get_or(const std::string& key, int default_value) const noexcept {
    size_type pos = find_position(key);
    if (pos < capacity_) {
        return table_[pos].value;
    }
    return default_value;
}

bool HashTable::erase(const std::string& key) noexcept {
    size_type pos = find_position(key);
    if (pos < capacity_) {
        table_[pos].occupied = false;
        table_[pos].key.clear();
        size_--;
        return true;
    }
    return false;
}

void HashTable::clear() noexcept {
    for (auto& entry : table_) {
        entry.occupied = false;
        entry.key.clear();
        entry.value = 0;
    }
    size_ = 0;
}

float HashTable::load_factor() const noexcept {
    if (capacity_ == 0) return 0.0f;
    return static_cast<float>(size_) / static_cast<float>(capacity_);
}
