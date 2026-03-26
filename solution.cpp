#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <limits>
#include <cstddef>
using namespace std;

const double EPS = 1e-9;
const double OUTPUT_PRECISION = 5;

struct Vec2 {
    double x, y;
    Vec2(double x_ = 0, double y_ = 0) : x(x_), y(y_) {}

    Vec2 operator+(const Vec2& o) const { return Vec2(x + o.x, y + o.y); }
    Vec2 operator-(const Vec2& o) const { return Vec2(x - o.x, y - o.y); }
    Vec2 operator*(double s) const { return Vec2(x * s, y * s); }
    Vec2 operator/(double s) const { return Vec2(x / s, y / s); }

    Vec2 operator-() const { return Vec2(-x, -y); }

    Vec2& operator+=(const Vec2& o) { x += o.x; y += o.y; return *this; }
    Vec2& operator-=(const Vec2& o) { x -= o.x; y -= o.y; return *this; }
    Vec2& operator*=(double s) { x *= s; y *= s; return *this; }

    double dot(const Vec2& o) const { return x * o.x + y * o.y; }
    double cross(const Vec2& o) const { return x * o.y - y * o.x; }

    double length() const { return sqrt(x * x + y * y); }
    double length_sq() const { return x * x + y * y; }

    Vec2 normalize() const {
        double len = length();
        if (len < EPS) return Vec2(0, 0);
        return *this / len;
    }

    Vec2 perp() const { return Vec2(-y, x); }

    Vec2 rotate90() const { return Vec2(-y, x); }

    bool operator==(const Vec2& o) const {
        return fabs(x - o.x) < EPS && fabs(y - o.y) < EPS;
    }
};

struct Polygon {
    vector<Vec2> v;

    Polygon() {}
    explicit Polygon(const vector<Vec2>& verts) : v(verts) {}

    int size() const { return (int)v.size(); }
    bool empty() const { return v.empty(); }

    Vec2 operator[](int i) const { return v[i]; }
    Vec2& operator[](int i) { return v[i]; }

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

    bool is_convex() const {
        int n = size();
        if (n < 3) return false;
        int sign = 0;
        for (int i = 0; i < n; ++i) {
            Vec2 a = v[i];
            Vec2 b = v[(i + 1) % n];
            Vec2 c = v[(i + 2) % n];
            double cr = (b - a).cross(c - b);
            if (fabs(cr) > EPS) {
                int s = cr > 0 ? 1 : -1;
                if (sign == 0) sign = s;
                else if (sign != s) return false;
            }
        }
        return true;
    }

    double area() const {
        int n = size();
        if (n < 3) return 0;
        double a = 0;
        for (int i = 0; i < n; ++i) {
            a += v[i].cross(v[(i + 1) % n]);
        }
        return a * 0.5;
    }

    int winding() const {
        return area() > 0 ? 1 : -1;
    }
};

struct Projection {
    double mn, mx;
    Projection(double min_val = 0, double max_val = 0) : mn(min_val), mx(max_val) {}

    double overlap(const Projection& o) const {
        return std::min(mx, o.mx) - std::max(mn, o.mn);
    }

    bool is_separated(const Projection& o) const {
        return mx <= o.mn || o.mx <= mn;
    }
};

Projection project(const Polygon& poly, const Vec2& axis) {
    if (poly.empty()) return Projection(0, 0);
    double mn = poly[0].dot(axis);
    double mx = mn;
    for (int i = 1; i < poly.size(); ++i) {
        double p = poly[i].dot(axis);
        if (p < mn) mn = p;
        if (p > mx) mx = p;
    }
    return Projection(mn, mx);
}

// Ear clipping for concave polygon decomposition
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

            // Check if convex vertex
            double cross = (b - a).cross(c - b);
            if (cross <= EPS) continue;

            // Check if ear is empty
            bool empty = true;
            for (int j = 0; j < n; ++j) {
                if (j == prev || j == i || j == next) continue;
                Vec2 p = verts[j];
                // Point in triangle test
                double s1 = (b - a).cross(p - a);
                double s2 = (c - b).cross(p - b);
                double s3 = (a - c).cross(p - c);
                if ((s1 > 0 && s2 > 0 && s3 > 0) || (s1 < 0 && s2 < 0 && s3 < 0)) {
                    empty = false;
                    break;
                }
            }

            if (empty) {
                // Found ear, add triangle
                result.push_back(Polygon({a, b, c}));
                // Remove vertex i
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

    if (n >= 3) {
        result.push_back(Polygon(verts));
    }

    return result.empty() ? vector<Polygon>{poly} : result;
}

// SAT-based MTV computation for convex polygons
// Returns {has_overlap, mtv}
pair<bool, Vec2> compute_mtv_sat_convex(const Polygon& A, const Polygon& B) {
    double min_overlap = INFINITY;
    Vec2 min_axis(0, 0);
    bool has_overlap = false;

    const Polygon* polys[2] = {&A, &B};

    for (int pi = 0; pi < 2; ++pi) {
        const Polygon& poly = *polys[pi];
        int n = poly.size();
        for (int i = 0; i < n; ++i) {
            auto [p1, p2] = poly.edge(i);
            Vec2 edge = p2 - p1;
            Vec2 axis = edge.perp().normalize();

            if (axis.length_sq() < EPS * EPS) continue;

            Projection pa = project(A, axis);
            Projection pb = project(B, axis);

            if (pa.is_separated(pb)) {
                return {false, Vec2(0, 0)};
            }

            double ov = pa.overlap(pb);
            if (ov < min_overlap) {
                min_overlap = ov;
                min_axis = axis;
                has_overlap = true;
            }
        }
    }

    if (!has_overlap) {
        return {false, Vec2(0, 0)};
    }

    // Ensure MTV direction is from A to B
    Vec2 dir = B.get_center() - A.get_center();
    if (min_axis.dot(dir) < 0) {
        min_axis = min_axis * -1;
    }

    return {true, min_axis * min_overlap};
}

// SAT with concave polygon support via convex decomposition
pair<bool, Vec2> compute_mtv_sat(const Polygon& A_orig, const Polygon& B_orig) {
    // Decompose concave polygons
    vector<Polygon> decompA = decompose_concave(A_orig);
    vector<Polygon> decompB = decompose_concave(B_orig);

    // If both are convex, use the fast path
    if (decompA.size() == 1 && decompB.size() == 1) {
        return compute_mtv_sat_convex(A_orig, B_orig);
    }

    // For concave polygons, check overlap between each pair of convex parts
    double best_overlap = INFINITY;
    Vec2 best_axis(0, 0);
    bool any_overlap = false;

    for (const Polygon& Ai : decompA) {
        for (const Polygon& Bj : decompB) {
            auto [has_ov, mtv] = compute_mtv_sat_convex(Ai, Bj);
            if (!has_ov) {
                // Parts don't overlap, continue
                continue;
            }
            any_overlap = true;
            double ov = mtv.length();
            if (ov < best_overlap) {
                best_overlap = ov;
                best_axis = mtv.length() > EPS ? mtv.normalize() : Vec2(0, 0);
            }
        }
    }

    if (!any_overlap) {
        return {false, Vec2(0, 0)};
    }

    return {true, best_axis * best_overlap};
}

// Point in polygon test (ray casting)
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
        double d = (p - a).length();
        return {d, (p - a).normalize()};
    }
    double t = max(0.0, min(1.0, (p - a).dot(ab) / len_sq));
    Vec2 closest = a + ab * t;
    double dist = (p - closest).length();
    Vec2 dir;
    if (dist > EPS) {
        dir = (p - closest) / dist;
    } else {
        dir = ab.perp().normalize();
    }
    return {dist, dir};
}

// Distance from point to polygon boundary, returns {distance, direction_to_nearest}
pair<double, Vec2> point_to_polygon_boundary(const Vec2& p, const Polygon& poly) {
    double min_dist = INFINITY;
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

// Build NFP for convex polygons using Minkowski sum approach
Polygon build_nfp_convex(const Polygon& A, const Polygon& B) {
    // Sort vertices by polar angle
    auto sort_by_angle = [](const Polygon& poly) {
        Vec2 center = poly.get_center();
        vector<pair<double, int>> idx;
        for (int i = 0; i < poly.size(); ++i) {
            double angle = atan2(poly[i].y - center.y, poly[i].x - center.x);
            idx.push_back({angle, i});
        }
        sort(idx.begin(), idx.end());
        vector<Vec2> sorted;
        for (auto& p : idx) sorted.push_back(poly[p.second]);
        return sorted;
    };

    vector<Vec2> sortedA = sort_by_angle(A);
    vector<Vec2> sortedB = sort_by_angle(B);
    int n = sortedA.size(), m = sortedB.size();

    // Minkowski sum: A ⊕ (-B)
    vector<Vec2> nfp;
    int i = 0, j = 0;
    while (i < n || j < m) {
        Vec2 sum = sortedA[i % n] + (-sortedB[j % m]);
        nfp.push_back(sum);

        Vec2 tangA = sortedA[(i + 1) % n] - sortedA[i % n];
        Vec2 tangB = sortedB[(j + 1) % m] - sortedB[j % m];
        double cross = tangA.cross(tangB);

        if (fabs(cross) < EPS) {
            i++; j++;
        } else if (cross > 0) {
            j++;
        } else {
            i++;
        }

        if (i >= n && j >= m) break;
    }

    // Remove duplicate consecutive points
    vector<Vec2> unique_nfp;
    for (const Vec2& p : nfp) {
        if (unique_nfp.empty() || (p - unique_nfp.back()).length_sq() > EPS * EPS) {
            unique_nfp.push_back(p);
        }
    }

    return Polygon(unique_nfp);
}

// Main solver class
class MTVSolver {
public:
    Polygon polyA, polyB;
    Polygon nfp;  // Precomputed NFP
    bool has_nfp = false;

    MTVSolver(const Polygon& A, const Polygon& B) : polyA(A), polyB(B) {
        build_nfp();
    }

    void build_nfp() {
        if (polyA.is_convex() && polyB.is_convex()) {
            nfp = build_nfp_convex(polyA, polyB);
            has_nfp = nfp.size() >= 3;
        }
    }

    Vec2 solve(const Vec2& displacement) {
        // Translate B by displacement
        Polygon B_shifted = polyB.translate(displacement);

        // SAT overlap check
        auto [has_overlap, mtv] = compute_mtv_sat(polyA, B_shifted);

        if (!has_overlap) {
            return Vec2(0, 0);
        }

        return mtv;
    }

    Vec2 solve_with_nfp(const Vec2& displacement) {
        // Translate B by displacement
        Polygon B_shifted = polyB.translate(displacement);

        // Check if reference point (displacement vector tip) is inside NFP
        // The reference point for B is at 'displacement' relative to B's original position
        // When we compute NFP = A ⊕ (-B), the reference point is the origin
        // If displacement is inside NFP, overlap occurs

        if (!has_nfp) {
            return solve(displacement);
        }

        if (!point_in_polygon(displacement, nfp)) {
            return Vec2(0, 0);
        }

        // Find closest point on NFP boundary
        auto [dist, dir] = point_to_polygon_boundary(displacement, nfp);

        // MTV is from displacement point to nearest boundary
        return dir * dist;
    }
};

// ============ I/O ============

bool try_read_ok(istream& in) {
    // Peek at next character to see if it's a number or "OK"
    // If it's a digit or '-', there's no OK marker (practice data)
    // If it's 'O', it's an OK marker
    int c = in.peek();
    if (c == 'O' || c == 'o') {
        string s;
        in >> s;
        return s == "OK";
    }
    // No OK marker, just continue
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

    // =========== 1. Read polygons ===========
    int n1, n2;
    if (!(cin >> n1 >> n2)) {
        return 0;
    }

    Polygon polyA = read_polygon(cin, n1);
    Polygon polyB = read_polygon(cin, n2);

    // Try to read OK marker (works for judge, harmless for practice data)
    try_read_ok(cin);

    // =========== 2. Preprocess ===========
    MTVSolver solver(polyA, polyB);

    cout << "OK" << endl;
    cout.flush();

    // =========== 3. Read test cases ===========
    int m;
    cin >> m;

    vector<Vec2> tests;
    tests.reserve(m);
    for (int i = 0; i < m; ++i) {
        double x, y;
        cin >> x >> y;
        tests.emplace_back(x, y);
    }

    // Try to read OK marker
    try_read_ok(cin);

    // =========== 4. Solve and output ===========
    cout << setprecision(OUTPUT_PRECISION) << fixed;
    cout << m << endl;
    cout.flush();

    for (int i = 0; i < m; ++i) {
        Vec2 res = solver.solve(tests[i]);
        cout << res.x << " " << res.y << endl;
        cout.flush();
    }

    cout << "OK" << endl;
    cout.flush();

    return 0;
}
