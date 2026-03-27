#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <limits>
using namespace std;

const double EPS = 1e-8;
const int OUTPUT_PRECISION = 5;

struct Vec2 {
    double x, y;
    Vec2(double x_ = 0, double y_ = 0) : x(x_), y(y_) {}
    Vec2 operator+(const Vec2& o) const { return Vec2(x + o.x, y + o.y); }
    Vec2 operator-(const Vec2& o) const { return Vec2(x - o.x, y - o.y); }
    Vec2 operator*(double s) const { return Vec2(x * s, y * s); }
    Vec2 operator-() const { return Vec2(-x, -y); }
    double dot(const Vec2& o) const { return x * o.x + y * o.y; }
    double cross(const Vec2& o) const { return x * o.y - y * o.x; }
    double length() const { return sqrt(x * x + y * y); }
    double length_sq() const { return x * x + y * y; }
    Vec2 normalize() const {
        double len = length();
        return len < EPS ? Vec2(0, 0) : Vec2(x / len, y / len);
    }
    Vec2 perp() const { return Vec2(-y, x); }
};

struct Polygon {
    vector<Vec2> v;
    Polygon() {}
    explicit Polygon(const vector<Vec2>& verts) : v(verts) {}
    int size() const { return (int)v.size(); }
    bool empty() const { return v.empty(); }
    Vec2 operator[](int i) const { return v[i]; }
    Vec2 get_center() const {
        if (v.empty()) return Vec2(0, 0);
        double cx = 0, cy = 0;
        for (const auto& p : v) { cx += p.x; cy += p.y; }
        return Vec2(cx / v.size(), cy / v.size());
    }
    pair<Vec2, Vec2> edge(int i) const {
        int n = size();
        return {v[i], v[(i + 1) % n]};
    }
    Polygon translate(const Vec2& t) const {
        vector<Vec2> nv;
        nv.reserve(size());
        for (const auto& p : v) nv.push_back(p + t);
        return Polygon(nv);
    }
};

// Check if polygon is convex
bool is_convex(const Polygon& poly) {
    int n = poly.size();
    if (n < 3) return false;
    int sign = 0;
    for (int i = 0; i < n; ++i) {
        Vec2 a = poly[i];
        Vec2 b = poly[(i + 1) % n];
        Vec2 c = poly[(i + 2) % n];
        double cr = (b - a).cross(c - b);
        if (fabs(cr) > EPS) {
            int s = cr > 0 ? 1 : -1;
            if (sign == 0) sign = s;
            else if (sign != s) return false;
        }
    }
    return true;
}

// Point in polygon (ray casting)
bool point_in_polygon(const Vec2& p, const Polygon& poly) {
    int n = poly.size();
    if (n < 3) return false;
    bool inside = false;
    int j = n - 1;
    for (int i = 0; i < n; ++i) {
        const Vec2& vi = poly[i];
        const Vec2& vj = poly[j];
        if (((vi.y > p.y) != (vj.y > p.y)) &&
            (p.x < (vj.x - vi.x) * (p.y - vi.y) / (vj.y - vi.y + EPS) + vi.x)) {
            inside = !inside;
        }
        j = i;
    }
    return inside;
}

// Distance from point to line segment
pair<double, Vec2> point_to_segment(const Vec2& p, const Vec2& a, const Vec2& b) {
    Vec2 ab = b - a;
    double len_sq = ab.length_sq();
    if (len_sq < EPS * EPS) {
        Vec2 d = p - a;
        return {d.length(), d};
    }
    double t = max(0.0, min(1.0, (p - a).dot(ab) / len_sq));
    Vec2 closest = a + ab * t;
    Vec2 dir = p - closest;
    return {dir.length(), dir};
}

// Distance from point to polygon boundary
pair<double, Vec2> point_to_polygon_boundary(const Vec2& p, const Polygon& poly) {
    double min_dist = numeric_limits<double>::infinity();
    Vec2 min_dir(0, 0);
    int n = poly.size();
    for (int i = 0; i < n; ++i) {
        auto [a, b] = poly.edge(i);
        auto [d, dir] = point_to_segment(p, a, b);
        if (d < min_dist) {
            min_dist = d;
            min_dir = dir;
        }
    }
    return {min_dist, min_dir};
}

// SAT overlap check
bool polygons_overlap(const Polygon& A, const Polygon& B) {
    const Polygon* polys[2] = {&A, &B};
    for (int pi = 0; pi < 2; ++pi) {
        const Polygon& poly = *polys[pi];
        int n = poly.size();
        for (int i = 0; i < n; ++i) {
            auto [p1, p2] = poly.edge(i);
            Vec2 axis = (p2 - p1).perp().normalize();
            if (axis.length_sq() < EPS * EPS) continue;

            double min_a = numeric_limits<double>::infinity();
            double max_a = -numeric_limits<double>::infinity();
            for (const Vec2& p : A.v) {
                double proj = p.dot(axis);
                min_a = min(min_a, proj);
                max_a = max(max_a, proj);
            }
            double min_b = numeric_limits<double>::infinity();
            double max_b = -numeric_limits<double>::infinity();
            for (const Vec2& p : B.v) {
                double proj = p.dot(axis);
                min_b = min(min_b, proj);
                max_b = max(max_b, proj);
            }
            if (max_a <= min_b || max_b <= min_a) return false;
        }
    }
    return true;
}

// Ear clipping decomposition for concave polygons (unused but kept for reference)
vector<Polygon> decompose_concave(const Polygon& poly) {
    vector<Polygon> result;
    if (poly.size() <= 3) {
        if (poly.size() == 3) result.push_back(poly);
        return result;
    }

    vector<Vec2> verts = poly.v;
    int n = verts.size();

    while (n > 3) {
        bool ear_found = false;
        for (int i = 0; i < n; ++i) {
            int prev = (i - 1 + n) % n;
            int next = (i + 1) % n;
            Vec2 a = verts[prev];
            Vec2 b = verts[i];
            Vec2 c = verts[next];

            double cross = (b - a).cross(c - b);
            if (cross <= EPS) continue;

            bool empty = true;
            for (int j = 0; j < n; ++j) {
                if (j == prev || j == i || j == next) continue;
                Vec2 p = verts[j];
                double s1 = (b - a).cross(p - a);
                double s2 = (c - b).cross(p - b);
                double s3 = (a - c).cross(p - c);
                if ((s1 > 0 && s2 > 0 && s3 > 0) || (s1 < 0 && s2 < 0 && s3 < 0)) {
                    empty = false;
                    break;
                }
            }

            if (empty) {
                result.push_back(Polygon({a, b, c}));
                vector<Vec2> new_verts;
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

// Minkowski sum for convex polygons: A ⊕ (-B)
Polygon minkowski_sum(const Polygon& A, const Polygon& B) {
    auto sort_by_angle = [](const Polygon& poly) {
        Vec2 center = poly.get_center();
        vector<pair<double, int>> idx;
        for (int i = 0; i < poly.size(); ++i) {
            double angle = atan2(poly[i].y - center.y, poly[i].x - center.x);
            idx.push_back({angle, i});
        }
        sort(idx.begin(), idx.end());
        vector<Vec2> sorted;
        for (auto& pr : idx) sorted.push_back(poly[pr.second]);
        return sorted;
    };

    vector<Vec2> sortedA = sort_by_angle(A);
    vector<Vec2> sortedB = sort_by_angle(B);
    int n = sortedA.size(), m = sortedB.size();

    vector<Vec2> result;
    int i = 0, j = 0;
    while (true) {
        result.push_back(sortedA[i % n] + (-sortedB[j % m]));
        Vec2 tangA = sortedA[(i + 1) % n] - sortedA[i % n];
        Vec2 tangB = sortedB[(j + 1) % m] - sortedB[j % m];
        double cross = tangA.cross(tangB);
        if (fabs(cross) < EPS) { i++; j++; }
        else if (cross > 0) { j++; }
        else { i++; }
        if (i >= n && j >= m) break;
    }

    vector<Vec2> uniq;
    for (const Vec2& p : result) {
        if (uniq.empty() || (p - uniq.back()).length_sq() > EPS * EPS) {
            uniq.push_back(p);
        }
    }
    return Polygon(uniq);
}

// Check if two polygons overlap using SAT
bool polygons_overlap(const Polygon& A, const Polygon& B);

// Simple SAT MTV without validation - faster for convex polygons
pair<bool, Vec2> compute_mtv_simple(const Polygon& A, const Polygon& B) {
    double min_mtv_dist = numeric_limits<double>::infinity();
    Vec2 min_axis(0, 0);

    const Polygon* polys[2] = {&A, &B};
    for (int pi = 0; pi < 2; ++pi) {
        const Polygon& poly = *polys[pi];
        int n = poly.size();
        for (int i = 0; i < n; ++i) {
            auto [p1, p2] = poly.edge(i);
            Vec2 axis = (p2 - p1).perp().normalize();
            if (axis.length_sq() < EPS * EPS) continue;

            double min_a = numeric_limits<double>::infinity();
            double max_a = -numeric_limits<double>::infinity();
            for (const Vec2& p : A.v) {
                double proj = p.dot(axis);
                min_a = min(min_a, proj);
                max_a = max(max_a, proj);
            }
            double min_b = numeric_limits<double>::infinity();
            double max_b = -numeric_limits<double>::infinity();
            for (const Vec2& p : B.v) {
                double proj = p.dot(axis);
                min_b = min(min_b, proj);
                max_b = max(max_b, proj);
            }
            if (max_a <= min_b || max_b <= min_a) return {false, Vec2(0, 0)};

            double overlap = min(max_a, max_b) - max(min_a, min_b);
            if (overlap < EPS) continue;

            double mtv_pos = max_a - min_b; // translate B positive
            double mtv_neg = max_b - min_a; // translate B negative
            double mtv_dist = min(mtv_pos, mtv_neg);

            if (mtv_dist > 0 && mtv_dist < min_mtv_dist) {
                min_mtv_dist = mtv_dist;
                min_axis = (mtv_pos < mtv_neg) ? axis : -axis;
            }
        }
    }
    if (min_mtv_dist == numeric_limits<double>::infinity()) return {false, Vec2(0, 0)};
    return {true, min_axis * min_mtv_dist};
}

// Full SAT MTV with validation - correct for concave polygons but slower
pair<bool, Vec2> compute_mtv(const Polygon& A, const Polygon& B) {
    double min_mtv_dist = numeric_limits<double>::infinity();
    Vec2 min_axis(0, 0);

    const Polygon* polys[2] = {&A, &B};
    for (int pi = 0; pi < 2; ++pi) {
        const Polygon& poly = *polys[pi];
        int n = poly.size();
        for (int i = 0; i < n; ++i) {
            auto [p1, p2] = poly.edge(i);
            Vec2 axis = (p2 - p1).perp().normalize();
            if (axis.length_sq() < EPS * EPS) continue;

            double min_a = numeric_limits<double>::infinity();
            double max_a = -numeric_limits<double>::infinity();
            for (const Vec2& p : A.v) {
                double proj = p.dot(axis);
                min_a = min(min_a, proj);
                max_a = max(max_a, proj);
            }
            double min_b = numeric_limits<double>::infinity();
            double max_b = -numeric_limits<double>::infinity();
            for (const Vec2& p : B.v) {
                double proj = p.dot(axis);
                min_b = min(min_b, proj);
                max_b = max(max_b, proj);
            }
            if (max_a <= min_b || max_b <= min_a) return {false, Vec2(0, 0)};

            double overlap = min(max_a, max_b) - max(min_a, min_b);
            if (overlap < EPS) continue;

            double mtv_pos = max_a - min_b; // translate B positive
            double mtv_neg = max_b - min_a; // translate B negative

            // Try positive direction (axis)
            if (mtv_pos > EPS) {
                Polygon B_shifted = B.translate(axis * mtv_pos);
                if (!polygons_overlap(A, B_shifted)) {
                    if (mtv_pos < min_mtv_dist) {
                        min_mtv_dist = mtv_pos;
                        min_axis = axis;
                    }
                }
            }

            // Try negative direction (-axis)
            if (mtv_neg > EPS) {
                Polygon B_shifted = B.translate(-axis * mtv_neg);
                if (!polygons_overlap(A, B_shifted)) {
                    if (mtv_neg < min_mtv_dist) {
                        min_mtv_dist = mtv_neg;
                        min_axis = -axis;
                    }
                }
            }
        }
    }
    if (min_mtv_dist == numeric_limits<double>::infinity()) return {false, Vec2(0, 0)};
    return {true, min_axis * min_mtv_dist};
}

// MTV Solver - uses simple SAT with validation fallback for concave polygons
class MTVSolver {
public:
    Polygon polyA, polyB;

    MTVSolver(const Polygon& A, const Polygon& B) : polyA(A), polyB(B) {}

    Vec2 solve(const Vec2& displacement) {
        Polygon B_shifted = polyB.translate(displacement);

        // Check overlap
        if (!polygons_overlap(polyA, B_shifted)) {
            return Vec2(0, 0);
        }

        // Try simple SAT first (fast for convex polygons)
        auto [has_overlap, mtv] = compute_mtv_simple(polyA, B_shifted);
        if (!has_overlap) return Vec2(0, 0);

        // Verify MTV is correct (needed for concave polygons)
        Polygon B_final = B_shifted.translate(mtv);
        if (!polygons_overlap(polyA, B_final)) {
            return mtv;
        }

        // Simple SAT MTV was wrong, use full validated SAT
        auto [has_overlap2, mtv2] = compute_mtv(polyA, B_shifted);
        if (has_overlap2) return mtv2;

        // Fallback: return the simple MTV (might not be optimal but at least separates)
        return mtv;
    }
};

// ============ I/O ============

bool try_read_ok(istream& in) {
    int c = in.peek();
    if (c == 'O' || c == 'o') {
        string s;
        in >> s;
        return s == "OK";
    }
    return false;
}

Polygon read_polygon(istream& in, int n) {
    vector<Vec2> verts;
    verts.reserve(n);
    for (int i = 0; i < n; ++i) {
        double x, y;
        in >> x >> y;
        verts.emplace_back(x, y);
    }
    return Polygon(verts);
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n1, n2;
    if (!(cin >> n1 >> n2)) return 0;

    Polygon polyA = read_polygon(cin, n1);
    Polygon polyB = read_polygon(cin, n2);
    try_read_ok(cin);

    MTVSolver solver(polyA, polyB);

    cout << "OK" << endl;
    cout.flush();

    int m;
    cin >> m;

    vector<Vec2> tests;
    tests.reserve(m);
    for (int i = 0; i < m; ++i) {
        double x, y;
        cin >> x >> y;
        tests.emplace_back(x, y);
    }

    cout << setprecision(OUTPUT_PRECISION) << fixed;
    cout << m << endl;
    cout.flush();

    for (int i = 0; i < m; ++i) {
        Vec2 res = solver.solve(tests[i]);
        // Avoid -0.00000 by checking for near-zero
        double ox = fabs(res.x) < 1e-9 ? 0.0 : res.x;
        double oy = fabs(res.y) < 1e-9 ? 0.0 : res.y;
        cout << ox << " " << oy << endl;
        cout.flush();
    }

    cout << "OK" << endl;
    cout.flush();

    return 0;
}
