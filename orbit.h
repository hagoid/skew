#ifndef SKEW_ORBIT_H
#define SKEW_ORBIT_H

#include <compare>
#include <cstddef>
#include <vector>

template<typename Scalar>
using TOrbitContainer = std::vector<Scalar>;

template<typename Scalar>
class TOrbit {
public:
    using Container = TOrbitContainer<Scalar>;
    class Iterator {
    public:
        using Container = TOrbitContainer<Scalar>;

        Iterator(const Container &container, std::size_t index, bool first);

        const Scalar &operator*() const;
        Iterator &operator++();

        friend bool operator!=(const Iterator &lhs, const Iterator &rhs) {
            return lhs.index != rhs.index || lhs.first != rhs.first || &lhs.container != &rhs.container;
        }

    private:
        const Container &container;
        std::size_t index;
        bool first;
    };

    TOrbit();
    explicit TOrbit(const Container *container);
    TOrbit(const Container *container, std::size_t position);

    [[nodiscard]] std::size_t size() const;
    [[nodiscard]] bool empty() const;
    const Scalar &operator[](std::size_t index) const;
    Iterator begin() const noexcept;
    Iterator end() const noexcept;

    const Container &container() const;

private:
    const Container *container_ = nullptr;
    std::size_t position_ = 0;
};

template<typename Scalar>
TOrbit<Scalar>::Iterator::Iterator(const Container &container, std::size_t index, bool first)
    : container(container)
    , index(index)
    , first(first) {
}

template<typename Scalar>
const Scalar &TOrbit<Scalar>::Iterator::Iterator::operator*() const {
    return container[index];
}

template<typename Scalar>
typename TOrbit<Scalar>::Iterator &TOrbit<Scalar>::Iterator::operator++() {
    ++index;
    if (index >= container.size()) {
        index -= container.size();
    }
    first = false;
    return *this;
}

template<typename Scalar>
TOrbit<Scalar>::TOrbit() = default;

template<typename Scalar>
TOrbit<Scalar>::TOrbit(const Container *container) : container_(container) {
}

template<typename Scalar>
TOrbit<Scalar>::TOrbit(const Container *container, std::size_t position) : container_(container), position_(position) {
}

template<typename Scalar>
std::size_t TOrbit<Scalar>::size() const {
    return container().size();
}

template<typename Scalar>
bool TOrbit<Scalar>::empty() const {
    return container_ == nullptr;
}

template<typename Scalar>
const Scalar &TOrbit<Scalar>::operator[](std::size_t index) const {
    auto index_ = index + position_;
    if (index_ >= size()) {
        index_ -= size();
    }
    return container()[index_];
}

template<typename Scalar>
typename TOrbit<Scalar>::Iterator TOrbit<Scalar>::begin() const noexcept {
    return {container(), position_, true};
}

template<typename Scalar>
typename TOrbit<Scalar>::Iterator TOrbit<Scalar>::end() const noexcept {
    return {container(), position_, false};
}

template<typename Scalar>
const typename TOrbit<Scalar>::Container &TOrbit<Scalar>::container() const {
    return *container_;
}

#endif  // SKEW_ORBIT_H
