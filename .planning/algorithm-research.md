# 算法改进研究

**日期**: 2026-03-28
**目标**: 更准、更快

## 当前方案 (Demo)

纯 SAT (Separating Axis Theorem)
- MTV = axis * overlap
- EPS = 1e-6
- 优点: 简单快速
- 缺点: 凹多边形可能错误

## 改进方案

### 1. AABB 宽相过滤
- 预处理: 计算两个多边形的轴对齐包围盒
- 每个测试点: 先检查 AABB 是否重叠
- 无重叠直接返回 (0,0)，跳过昂贵 SAT 计算
- **效果**: 快速过滤非重叠情况

### 2. 凸分解 (Convex Decomposition)
- Hertel-Mehlhorn 启发式算法 (线性时间)
- 将凹多边形分解为凸多边形集合
- 避免全局 SAT 处理凹多边形的错误

### 3. 闵可夫斯基和 (Minkowski Sum)
- 对凸多边形: A ⊕ (-B) 的边界就是 NFP
- MTV = B参考点指向NFF边界的最短向量
- 对凹多边形: 分解后各自计算，再取最小

### 4. 去三角函数化
- 使用叉积 (Cross Product) 判断凸凹性
- 避免 atan2, sin, cos 等昂贵计算
- 提高精度和速度

### 5. I/O 优化
- ios::sync_with_stdio(false)
- cin.tie(nullptr)
- 每次输出后 flush

## 提交记录
- `1907554`: Add AABB broad-phase, convex decomposition, Minkowski sum, cross-product geometry

## 性能分析
- AABB: O(1) 每测试点
- 凸分解: O(n) 预处理
- SAT: O((n+m) * k) 其中 k 是凸块数量
- 预期: 凸多边形快，凹多边形准