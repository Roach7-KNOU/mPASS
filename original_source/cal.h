#ifndef CAL_H
#define CAL_H

#include <cmath>
#include "pass.h"

inline std::string convertInt(int number) {
    std::stringstream ss;
    ss << number;
    return ss.str();
}

inline Vector3 vectorSum(const Vector3& a, const Vector3& b) {
    Vector3 c;
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    c.z = a.z + b.z;
    return c;
}

inline Vector3 vectorDiff(const Vector3& a, const Vector3& b) {
    Vector3 c;
    c.x = a.x - b.x;
    c.y = a.y - b.y;
    c.z = a.z - b.z;
    return c;
}

inline double vectorDot(const Vector3& a, const Vector3& b) {
    return (a.x * b.x + a.y * b.y + a.z * b.z);
}

inline Vector3 vectorCross(const Vector3& a, const Vector3& b) {
    Vector3 c;
    c.x = a.y * b.z - a.z * b.y;
    c.y = a.z * b.x - a.x * b.z;
    c.z = a.x * b.y - a.y * b.x;
    return c;
}

inline Vector3 vectorScale(const Vector3& a, double b) {
    Vector3 c;
    c.x = a.x * b;
    c.y = a.y * b;
    c.z = a.z * b;
    return c;
}

inline double vectorDistSq(const Vector3& a) {
    return (a.x * a.x + a.y * a.y + a.z * a.z);
}

inline double vectorDist(const Vector3& a) {
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

inline double vsDistSq(const Vector3& a, const Vector3& b) {
    return ((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

inline double vsDist(const Vector3& a, const Vector3& b) {
    return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

inline Vector3 normalizeVector(const Vector3& a) {
    double b = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    Vector3 c;
    c.x = a.x / b;
    c.y = a.y / b;
    c.z = a.z / b;
    return c;
}

inline Vector3 vectorAverage(const Vector3& a, const int b) {
    Vector3 c;
    c.x = a.x / b;
    c.y = a.y / b;
    c.z = a.z / b;
    return c;
}

inline double distNProbeProbe(Vector3 point, int probeId) {
    return vsDist(point, pAtom(probeId));
}

inline double distPoint(double x, double y, double z, const Vector3& b) {
    return sqrt((x - b.x) * (x - b.x) + (y - b.y) * (y - b.y) + (z - b.z) * (z - b.z));
}

inline double distNProbeProtein(Vector3 point, int proteinId) {
    return vsDist(point, prAtom(proteinId));
}

inline unsigned int getPairKey(unsigned short uIdxA, unsigned short uIdxB) {
    unsigned short uMin, uMax;
    if (uIdxA < uIdxB) {
        uMin = uIdxA; uMax = uIdxB;
    } else {
        uMin = uIdxB; uMax = uIdxA;
    }
    return ((PairKey)((uMin & 0xffff) | ((uMax & 0xffff) << 16)));
}

int initializeCellList();
double getCellId(Vector3 point, double size); 

#endif
