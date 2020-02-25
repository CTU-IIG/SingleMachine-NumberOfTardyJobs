#ifndef B_UTIL_H
#define B_UTIL_H

#include <algorithm>

// Prevzato z https://stackoverflow.com/a/9729747
template<class T>
inline void hash_combine(std::size_t &seed, const T &v) {
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

// Prevzato z https://stackoverflow.com/a/9729747
namespace std {
    template<typename S, typename T>
    struct hash<pair<S, T>> {
        inline size_t operator()(const pair<S, T> &v) const {
            size_t seed = 0;
            ::hash_combine(seed, v.first);
            ::hash_combine(seed, v.second);
            return seed;
        }
    };
}

// Prevzato z https://github.com/nicolausYes/iterator-template-sort-library/blob/master/src/InsertionSort.h
namespace sortings {
    class InsertionSort {
    public:
        template<class RandomAccessIterator, class Compare>
        static void sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            for (auto i = first; i != last; i++) {
                for (auto j = i; j != first; j--) {
                    if (comp(*j, *(j - 1)))
                        std::swap(*j, *(j - 1));
                    else
                        break;
                }
            }
        }

    private:
        InsertionSort(void);

        ~InsertionSort(void);
    };
}

#endif //B_UTIL_H
