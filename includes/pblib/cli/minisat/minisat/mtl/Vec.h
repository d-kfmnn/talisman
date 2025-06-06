/*******************************************************************************************[Vec.h]
Copyright (c) 2003-2007, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#ifndef Minisat_Vec_h
#define Minisat_Vec_h

#include <assert.h>
#include <new>
#include <type_traits>
#include <algorithm>
#include <functional>

#include "minisat/mtl/IntTypes.h"
#include "minisat/mtl/XAlloc.h"

namespace Minisat {

//=================================================================================================
// Automatically resizable arrays
//
// NOTE! Don't use this vector on datatypes that cannot be re-located in memory (with realloc)

template<class T>
class vec {
    using size_type = int;
    T* m_data = nullptr;
    T* m_current = nullptr;
    T* m_end = nullptr;

    size_type adjust_next_size(size_type requested) {
        const auto old_cap = capacity();
        if (std::numeric_limits<size_type>::max() - old_cap < old_cap / 2) {
            // Geometric growth would overflow
            return requested;
        }
        // At least 1.5x growth factor
        const auto geometric = old_cap + old_cap / 2;
        return std::max(geometric, requested);
    }

public:
    // other meta things
    using iterator = T * ;
    using const_iterator = T const*;
    using value_type = T;

    // Constructors:
    vec() = default;
    explicit vec(int size)      { growTo(size); }
    vec(int size, const T& pad) { growTo(size, pad); }
   ~vec()                       { clear(true); }

   // Don't allow copying (error prone):
   vec<T>&  operator = (vec<T> const&) = delete;
   vec(vec<T> const&) = delete;

    // Pointer to first element:
    operator T*() { return m_data; }

    // Size operations:
    int size() const { return static_cast<int>(m_current - m_data); }
    bool empty() const { return m_current == m_data; }
    void     shrink   (int nelems) {
        assert(nelems <= size() && "Attempted to shrink vector by more than it has elements");
        shrink_(nelems, std::is_trivially_destructible<T>{});
    }
    // Truncate the vector from position till the `end()`
    void truncate(const_iterator from) {
        // By using `std::less_equal` instead of `<=`, we have guaranteed
        // the definedness of the comparison even if the iterators (pointers)
        // are from a different allocation.
        assert(std::less_equal<const_iterator>{}(m_data, from) && "Outside of vector");
        assert(std::less_equal<const_iterator>{}(from, m_current) && "Outside of valid elements");
        erase_(from, std::is_trivially_destructible<T>{});
        m_current = const_cast<T*>(from);
    }
    int      capacity () const     { return static_cast<int>(m_end - m_data); }
    void     capacity (int min_cap);
    void     growTo   (int size);
    void     growTo   (int size, const T& pad);
    void     clear    (bool dealloc = false);

    // Stack interface:
    void push() {
        if (m_current == m_end) {
            capacity(size() + 1);
        }
        new (m_current) T();
        ++m_current;
    }
    void push (const T& elem) {
        if (m_current == m_end) {
            capacity(size() + 1);
        }
        *m_current = elem;
        ++m_current;
    }
    // Does not realloc if there is not enough space for the element!
    void     push_ (const T& elem)     { assert(m_current != m_end); *(m_current++) = elem; }
    void     pop() { assert(m_data != m_current && "Attempted vec::pop on an empty vec"); pop_(std::is_trivially_destructible<T>{}); }
    // NOTE: it seems possible that overflow can happen in the 'sz+1' expression of 'push()', but
    // in fact it can not since it requires that 'cap' is equal to INT_MAX. This in turn can not
    // happen given the way capacities are calculated (below). Essentially, all capacities are
    // even, but INT_MAX is odd.

    const T& last  () const { assert(m_data != m_current && "Called vec::last on an empty vector."); return *(m_current - 1); }
    T&       last  ()       { assert(m_data != m_current && "Called vec::last on an empty vector."); return *(m_current - 1); }

    // Vector interface:
    const T& operator [] (int index) const { return m_data[index]; }
    T&       operator [] (int index)       { return m_data[index]; }

    // Duplication (preferred instead):
    void copyTo(vec<T>& copy) const {
        copy.clear();
        copy.capacity(size());
        std::copy(begin(), end(), copy.begin());
        copy.m_current = copy.m_data + size();
    }
    void moveTo(vec<T>& dest) {
        swap(dest);
    }

    // Iterator support
    iterator begin() { return m_data; }
    iterator end() { return m_current; }
    const_iterator begin() const { return m_data; }
    const_iterator end() const { return m_current; }
    const_iterator cbegin() const { return m_data; }
    const_iterator cend() const { return m_current; }

    void swap(vec& rhs) {
        std::swap(m_data, rhs.m_data);
        std::swap(m_current, rhs.m_current);
        std::swap(m_end, rhs.m_end);
    }

private:
    //////
    // Specializations of implementation details based on T's properties

    // Is trivially destructible
    void erase_(const_iterator, std::true_type) {}
    // Is not trivially destructible
    void erase_(const_iterator from, std::false_type) {
        while (from != m_end) {
            from->~T();
            ++from;
        }
    }

    // Is trivially destructible
    void shrink_(int nelems, std::true_type) {
        m_current -= nelems;
    }

    // Is not trivially destructible
    void shrink_(int nelems, std::false_type) {
        while (nelems--) {
            --m_current;
            m_current->~T();
        }
    }

    // Is trivially destructible
    void pop_(std::true_type) { --m_current; }
    // Is not trivially destructible
    void pop_(std::false_type) { --m_current; m_current->~T(); }

    // Is trivially destructible
    void clear_(std::true_type) {}
    // Is trivially destructible
    void clear_(std::false_type) {
        if (m_data != NULL) {
            auto cptr = m_data;
            while (cptr != m_current) {
                cptr->~T();
                ++cptr;
            }
        }
    }

};


template<class T>
void vec<T>::capacity(int min_cap) {
    if (capacity() >= min_cap) return;

    const auto new_cap = adjust_next_size(min_cap);

    auto new_data = static_cast<T*>(realloc(m_data, new_cap * sizeof(T)));
    if (new_data == nullptr) {
        throw OutOfMemoryException("vec::capacity could not allocate enough memory");
    }
    const auto current_size = size();
    m_data = new_data;
    m_current = m_data + current_size;
    m_end = m_data + new_cap;
 }


template<class T>
void vec<T>::growTo(int new_size, const T& pad) {
    // If we are already larger, don't do anything
    if (size() >= new_size) return;

    capacity(new_size);
    const auto target = m_data + new_size;
    while (m_current != target) {
        *(m_current++) = pad;
    }
}


template<class T>
void vec<T>::growTo(int new_size) {
    // If we are already larger, don't do anything
    if (size() >= new_size) return;

    capacity(new_size);
    const auto target = m_data + new_size;
    while (m_current != target) {
        new (m_current++) T{};
    }
}


template<class T>
void vec<T>::clear(bool dealloc) {
    clear_(std::is_trivially_destructible<T>{});
    m_current = m_data;
    if (dealloc) {
        free(m_data);
        m_data = nullptr;
        m_current = nullptr;
        m_end = nullptr;
    }
}

template <class T>
void swap(vec<T>& lhs, vec<T>& rhs) {
    lhs.swap(rhs);
}

//=================================================================================================
}

#endif
