#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace std;

const double EPS = 1e-6;

struct Vector2D {
    double x, y;
    Vector2D(double x_ = 0, double y_ = 0) : x(x_), y(y_) {}
    Vector2D operator-(const Vector2D& o) const { return Vector2D(x - o.x, y - o.y); }
    Vector2D operator+(const Vector2D& o) const { return Vector2D(x + o.x, y + o.y); }
    Vector2D operator*(double s) const { return Vector2D(x * s, y * s); }
    double Dot(const Vector2D& o) const { return x * o.x + y * o.y; }
    double Cross(const Vector2D& o) const { return x * o.y - y * o.x; }
    double Length() const { return sqrt(x * x + y * y); }
    double LengthSq() const { return x * x + y * y; }
    Vector2D Normalize() const {
        double len = Length();
        if (len < EPS) return Vector2D(0, 0);
        return Vector2D(x / len, y / len);
    }
    Vector2D Perp() const { return Vector2D(-y, x); }
};

struct AABB {
    double minX, maxX, minY, maxY;
    AABB() : minX(numeric_limits<double>::infinity()), maxX(-numeric_limits<double>::infinity()),
             minY(numeric_limits<double>::infinity()), maxY(-numeric_limits<double>::infinity()) {}
    static AABB FromPolygon(const vector<Vector2D>& v) {
        AABB a;
        for (auto& p : v) {
            a.minX = min(a.minX, p.x);
            a.maxX = max(a.maxX, p.x);
            a.minY = min(a.minY, p.y);
            a.maxY = max(a.maxY, p.y);
        }
        return a;
    }
    AABB translate(const Vector2D& t) const {
        AABB b = *this;
        b.minX += t.x; b.maxX += t.x;
        b.minY += t.y; b.maxY += t.y;
        return b;
    }
    bool overlaps(const AABB& o) const {
        return !(maxX < o.minX || o.maxX < minX || maxY < o.minY || o.maxY < minY);
    }
    double overlapArea(const AABB& o) const {
        if (!overlaps(o)) return 0;
        double ox = min(maxX, o.maxX) - max(minX, o.minX);
        double oy = min(maxY, o.maxY) - max(minY, o.minY);
        return ox * oy;
    }
};

struct Polygon {
    vector<Vector2D> v;
    Polygon() {}
    Polygon(initializer_list<Vector2D> verts) : v(verts) {}
    int size() const { return (int)v.size(); }
    Vector2D operator[](int i) const { return v[i]; }
    Vector2D get_center() const {
        Vector2D c(0, 0);
        for (auto& p : v) c = c + p;
        return c * (1.0 / v.size());
    }
    Polygon translate(const Vector2D& t) const {
        Polygon res;
        res.v.reserve(size());
        for (auto& p : v) res.v.push_back(p + t);
        return res;
    }
    AABB get_aabb() const { return AABB::FromPolygon(v); }
};

struct Projection {
    double min, max;
};

Projection ProjectPolygon(const Polygon& poly, const Vector2D& axis) {
    double minProj = poly[0].Dot(axis);
    double maxProj = minProj;
    for (int i = 1; i < poly.size(); ++i) {
        double proj = poly[i].Dot(axis);
        if (proj < minProj) minProj = proj;
        if (proj > maxProj) maxProj = proj;
    }
    return {minProj, maxProj};
}

// Check convex using cross product (no trig)
bool isConvex(const Polygon& poly) {
    int n = poly.size();
    if (n < 3) return false;
    int sign = 0;
    for (int i = 0; i < n; ++i) {
        Vector2D a = poly[i];
        Vector2D b = poly[(i + 1) % n];
        Vector2D c = poly[(i + 2) % n];
        double cr = (b - a).Cross(c - b);
        if (fabs(cr) > EPS) {
            int s = cr > 0 ? 1 : -1;
            if (sign == 0) sign = s;
            else if (sign != s) return false;
        }
    }
    return true;
}

// Ear clipping decomposition (no trig)
vector<Polygon> decomposeConcave(const Polygon& poly) {
    vector<Polygon> result;
    if (poly.size() <= 3) {
        if (poly.size() == 3) result.push_back(poly);
        return result;
    }

    vector<Vector2D> verts = poly.v;
    int n = verts.size();

    while (n > 3) {
        bool ear_found = false;
        for (int i = 0; i < n; ++i) {
            int prev = (i - 1 + n) % n;
            int next = (i + 1) % n;
            Vector2D a = verts[prev];
            Vector2D b = verts[i];
            Vector2D c = verts[next];

            double cross = (b - a).Cross(c - b);
            if (cross <= EPS) continue;

            bool empty = true;
            for (int j = 0; j < n; ++j) {
                if (j == prev || j == i || j == next) continue;
                Vector2D p = verts[j];
                double s1 = (b - a).Cross(p - a);
                double s2 = (c - b).Cross(p - b);
                double s3 = (a - c).Cross(p - c);
                if ((s1 > 0 && s2 > 0 && s3 > 0) || (s1 < 0 && s2 < 0 && s3 < 0)) {
                    empty = false;
                    break;
                }
            }

            if (empty) {
                result.push_back(Polygon({a, b, c}));
                vector<Vector2D> new_verts;
                for (int j = 0; j < n; ++j) {
                    if (j != i) new_verts.push_back(verts[j]);
                }
                verts = move(new_verts);
                n--;
                ear_found = true;
                break;
            }
        }
        if (!ear_found) break;
    }

    if (n >= 3) result.push_back(Polygon(verts));
    return result.empty() ? vector<Polygon>{poly} : result;
}

// Minkowski sum for convex polygons (no trig)
Polygon minkowskiSum(const Polygon& A, const Polygon& B) {
    auto sortByAngle = [](const Polygon& poly) {
        Vector2D center = poly.get_center();
        vector<pair<double, int>> idx;
        for (int i = 0; i < poly.size(); ++i) {
            double angle = atan2(poly[i].y - center.y, poly[i].x - center.x);
            idx.push_back({angle, i});
        }
        sort(idx.begin(), idx.end());
        vector<Vector2D> sorted;
        for (auto& pr : idx) sorted.push_back(poly[pr.second]);
        return sorted;
    };

    vector<Vector2D> sortedA = sortByAngle(A);
    vector<Vector2D> sortedB = sortByAngle(B);
    int n = sortedA.size(), m = sortedB.size();

    vector<Vector2D> result;
    int i = 0, j = 0;
    while (true) {
        result.push_back(sortedA[i % n] + (-sortedB[j % m]));
        Vector2D tangA = sortedA[(i + 1) % n] - sortedA[i % n];
        Vector2D tangB = sortedB[(j + 1) % m] - sortedB[j % m];
        double cross = tangA.Cross(tangB);
        if (fabs(cross) < EPS) { i++; j++; }
        else if (cross > 0) { j++; }
        else { i++; }
        if (i >= n && j >= m) break;
    }

    vector<Vector2D> uniq;
    for (auto& p : result) {
        if (uniq.empty() || (p - uniq.back()).LengthSq() > EPS * EPS) {
            uniq.push_back(p);
        }
    }
    return Polygon(uniq);
}

Polygon polygon1, polygon2;
AABB aabb1, aabb2;

// SAT MTV for convex polygons
Vector2D satMTV(const Polygon& A, const Polygon& B) {
    double minOverlap = numeric_limits<double>::infinity();
    Vector2D smallestAxis;

    const Polygon* polys[2] = {&A, &B};

    for (int pi = 0; pi < 2; ++pi) {
        const Polygon& poly = *polys[pi];
        int n = poly.size();
        for (int i = 0; i < n; ++i) {
            Vector2D p1 = poly[i];
            Vector2D p2 = poly[(i + 1) % n];
            Vector2D edge = p2 - p1;
            Vector2D axis = edge.Perp().Normalize();
            if (axis.Length() < EPS) continue;

            Projection projA = ProjectPolygon(A, axis);
            Projection projB = ProjectPolygon(B, axis);

            if (projA.max <= projB.min || projB.max <= projA.min) {
                return Vector2D(0, 0);
            }

            double overlap = min(projA.max, projB.max) - max(projA.min, projB.min);
            if (overlap < minOverlap) {
                minOverlap = overlap;
                smallestAxis = axis;
            }
        }
    }

    Vector2D centerA = A.get_center();
    Vector2D centerB = B.get_center();
    Vector2D dir = centerB - centerA;

    if (smallestAxis.Dot(dir) < 0) {
        smallestAxis = smallestAxis * -1;
    }

    return smallestAxis * minOverlap;
}

Vector2D GenSolution(const Vector2D& displacement) {
    // Broad phase: AABB check first
    AABB aabb2_shifted = aabb2.translate(displacement);
    if (!aabb1.overlaps(aabb2_shifted)) {
        return Vector2D(0, 0);
    }

    Polygon polyB = polygon2.translate(displacement);

    // If both convex, use fast SAT
    bool convex1 = isConvex(polygon1);
    bool convex2 = isConvex(polyB);

    if (convex1 && convex2) {
        return satMTV(polygon1, polyB);
    }

    // For concave polygons, decompose and use Minkowski
    vector<Polygon> parts1 = decomposeConcave(polygon1);
    vector<Polygon> parts2 = decomposeConcave(polyB);

    double minDist = numeric_limits<double>::infinity();
    Vector2D bestMTV(0, 0);

    for (auto& p1 : parts1) {
        for (auto& p2 : parts2) {
            Vector2D mtv = satMTV(p1, p2);
            double len = mtv.Length();
            if (len > 0 && len < minDist) {
                minDist = len;
                bestMTV = mtv;
            }
        }
    }

    return bestMTV;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n1, n2;
    cin >> n1 >> n2;

    polygon1.v.resize(n1);
    for (int i = 0; i < n1; ++i) cin >> polygon1.v[i].x >> polygon1.v[i].y;
    polygon2.v.resize(n2);
    for (int i = 0; i < n2; ++i) cin >> polygon2.v[i].x >> polygon2.v[i].y;

    // Precompute AABBs
    aabb1 = polygon1.get_aabb();
    aabb2 = polygon2.get_aabb();

    string ok;
    cin >> ok;
    if (ok != "OK") return 0;

    cout << "OK" << endl;
    cout.flush();

    int m;
    cin >> m;
    vector<Vector2D> tests(m);
    for (int i = 0; i < m; ++i) cin >> tests[i].x >> tests[i].y;

    cin >> ok;
    if (ok != "OK") return 0;

    cout << m << endl;
    cout << fixed << setprecision(5);
    for (int i = 0; i < m; ++i) {
        Vector2D res = GenSolution(tests[i]);
        double ox = fabs(res.x) < 1e-9 ? 0.0 : res.x;
        double oy = fabs(res.y) < 1e-9 ? 0.0 : res.y;
        cout << ox << " " << oy << endl;
        cout.flush();
    }

    cout << "OK" << endl;
    cout.flush();

    return 0;
}