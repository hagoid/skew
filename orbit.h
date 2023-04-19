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
    using Iterator = typename TOrbitContainer<Scalar>::const_iterator;
    using Container = TOrbitContainer<Scalar>;

    TOrbit();
    explicit TOrbit(const Container *container);

    [[nodiscard]] std::size_t size() const;
    [[nodiscard]] bool empty() const;
    const Scalar &operator[](std::size_t index) const;
    Iterator begin() const noexcept;
    Iterator end() const noexcept;

    const Container &container() const;

    friend std::strong_ordering operator<=>(const TOrbit<Scalar> &lhs, const TOrbit<Scalar> &rhs) {
        return lhs.container() <=> rhs.container();
    }

    friend bool operator==(const TOrbit<Scalar> &lhs, const TOrbit<Scalar> &rhs) {
        return lhs.container() == rhs.container();
    }
private:
    const Container *container_;
};

template<typename Scalar>
TOrbit<Scalar>::TOrbit() : TOrbit(nullptr) {
}

template<typename Scalar>
TOrbit<Scalar>::TOrbit(const TOrbit::Container *container) : container_(container) {
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
    return container()[index];
}

template<typename Scalar>
typename TOrbit<Scalar>::Iterator TOrbit<Scalar>::begin() const noexcept {
    return container().begin();
}

template<typename Scalar>
typename TOrbit<Scalar>::Iterator TOrbit<Scalar>::end() const noexcept {
    return container().end();
}

template<typename Scalar>
const typename TOrbit<Scalar>::Container &TOrbit<Scalar>::container() const {
    return *container_;
}

#endif  // SKEW_ORBIT_H
