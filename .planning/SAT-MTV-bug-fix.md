# SAT MTV 算法 Bug 修复记录

**日期**: 2026-03-27
**问题**: Wrong Answer (格式/算法双重问题)
**状态**: 已修复并推送

## 发现的 Bug

### Bug 1: 输出格式 `-0.00000`
**位置**: `Solution.cpp` main() 输出循环
**原因**: 浮点数计算产生 -0.0 而非 0.0
**修复**: clamp 接近零的值 (`fabs(x) < 1e-9 -> 0.0`)

### Bug 2: SAT MTV 计算错误 (主要问题)
**位置**: `compute_mtv()` 函数
**原因**: 错误使用 `overlap` 而非 `gap + overlap`

错误公式:
```cpp
double overlap = min(max_a, max_b) - max(min_a, min_b);
MTV = axis * overlap;  // 错误！
```

正确公式:
```cpp
double mtv_pos = max_a - min_b;  // 沿正轴分离的距离
double mtv_neg = max_b - min_a;  // 沿负轴分离的距离
double mtv_dist = min(mtv_pos, mtv_neg);
MTV = axis * mtv_dist;
```

**验证**: 用测试用例 (0.44714, -0.10372) 验证：
- 修复前 MTV: (0, -0.08603) → 分离后仍重叠 (错误)
- 修复后 MTV: (0, -0.11349) → 一次性完全分离 (正确)

## 提交记录
- `e4056dd`: Fix -0.00000 output by clamping near-zero values
- `9cd6119`: Fix SAT MTV algorithm: use gap+overlap instead of overlap only

## 关键代码 (C++)

```cpp
pair<bool, Vec2> compute_mtv(const Polygon& A, const Polygon& B) {
    double min_mtv_dist = numeric_limits<double>::infinity();
    Vec2 min_axis(0, 0);
    bool has_overlap = false;

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
                has_overlap = true;
            }
        }
    }
    if (!has_overlap) return {false, Vec2(0, 0)};
    return {true, min_axis * min_mtv_dist};
}
```

## 教训
- SAT MTV 不能简单用 overlap，要用 gap + overlap
- Python 验证发现了 C++ 算法错误
- 测试用例验证是调试算法的关键

## Bug 3: SAT 对凹多边形MTV方向错误 (已修复)
**发现日期**: 2026-03-27
**问题**: 简单SAT对凹多边形返回的MTV方向可能错误（指向内部而非外部）
**原因**: 凹多边形的部分边法向量指向内部，min(mtv_pos, mtv_neg)可能选择错误方向
**验证**: practice_3 (凹多边形B, 172顶点) 测试10000个位移，简单SAT有2432个失败
**修复**: 添加MTV验证 - 对每个MTV候选，验证其是否真正分离多边形
**修复后**: 10000个测试全部通过

### Python验证代码
```python
def compute_mtv_sat(poly_a, poly_b, validate=True):
    # ... 对每个轴 ...
    mtv_pos = proj_a.max - proj_b.min  # 沿正轴分离
    mtv_neg = proj_b.max - proj_a.min  # 沿负轴分离

    if not validate:
        # 简单SAT：对凹多边形可能错误
        mtv_dist = min(mtv_pos, mtv_neg)
        ...
    else:
        # 验证MTV是否真正分离
        if mtv_pos > EPS:
            test_b = poly_b.translate(axis * mtv_pos)
            if not polygons_overlap(poly_a, test_b):
                # 验证通过，使用此MTV
                ...

        if mtv_neg > EPS:
            test_b = poly_b.translate(-axis * mtv_neg)
            if not polygons_overlap(poly_a, test_b):
                ...
```

### 两阶段策略 (更高效)
1. 先用简单SAT (validate=False) 计算MTV
2. 验证MTV是否真正分离
3. 如果验证失败，用完整验证SAT重新计算

### 性能
- Python: ~209秒处理10000个凹多边形测试 (48 tests/sec)
- C++: 预计<1秒 (原生编译)
- 优化策略有效: 大部分测试用例(凸多边形)使用快速简单SAT