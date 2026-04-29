#pragma once
// geometry.hpp – Modern C++17 Vector3 math utilities
// Replaces cal.h inline functions + macro-based SQUARE/MIN/MAX

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace mpass {

// ─── Vector3 ────────────────────────────────────────────────────────────────

struct Vector3 {
    double x{}, y{}, z{};

    constexpr Vector3() noexcept = default;
    constexpr Vector3(double x, double y, double z) noexcept : x(x), y(y), z(z) {}

    [[nodiscard]] constexpr Vector3 operator+(const Vector3& rhs) const noexcept {
        return {x + rhs.x, y + rhs.y, z + rhs.z};
    }
    [[nodiscard]] constexpr Vector3 operator-(const Vector3& rhs) const noexcept {
        return {x - rhs.x, y - rhs.y, z - rhs.z};
    }
    [[nodiscard]] constexpr Vector3 operator*(double scalar) const noexcept {
        return {x * scalar, y * scalar, z * scalar};
    }
    [[nodiscard]] constexpr Vector3 operator/(double scalar) const noexcept {
        return {x / scalar, y / scalar, z / scalar};
    }

    [[nodiscard]] constexpr double dot(const Vector3& rhs) const noexcept {
        return x * rhs.x + y * rhs.y + z * rhs.z;
    }
    [[nodiscard]] constexpr Vector3 cross(const Vector3& rhs) const noexcept {
        return {y * rhs.z - z * rhs.y,
                z * rhs.x - x * rhs.z,
                x * rhs.y - y * rhs.x};
    }

    [[nodiscard]] double lengthSquared() const noexcept { return dot(*this); }
    [[nodiscard]] double length()        const noexcept { return std::sqrt(lengthSquared()); }

    [[nodiscard]] Vector3 normalized() const {
        double len = length();
        if (len < 1e-15)
            throw std::domain_error("Cannot normalize zero-length vector");
        return *this / len;
    }

    [[nodiscard]] constexpr Vector3 average(int count) const noexcept {
        return *this / static_cast<double>(count);
    }
};

// ─── Free distance helpers ───────────────────────────────────────────────────

[[nodiscard]] inline double distanceSquared(const Vector3& a, const Vector3& b) noexcept {
    return (a - b).lengthSquared();
}

[[nodiscard]] inline double distance(const Vector3& a, const Vector3& b) noexcept {
    return std::sqrt(distanceSquared(a, b));
}

// Distance from a raw coordinate triple to a Vector3
[[nodiscard]] inline double distanceFromCoords(double px, double py, double pz,
                                               const Vector3& b) noexcept {
    return distance({px, py, pz}, b);
}

// ─── Angle between three points (vertex at p1) ──────────────────────────────

[[nodiscard]] inline double angleDegrees(const Vector3& p1, const Vector3& p2,
                                         const Vector3& p3) noexcept {
    Vector3 v1 = (p2 - p1).normalized();
    Vector3 v2 = (p3 - p1).normalized();
    double  dot = std::clamp(v1.dot(v2), -1.0, 1.0);
    constexpr double PI = 3.141592653589793;
    return std::acos(dot) * 180.0 / PI;
}

// ─── Square helper (replaces SQUARE macro) ──────────────────────────────────

template <typename T>
[[nodiscard]] constexpr T sq(T v) noexcept { return v * v; }

} // namespace mpass
