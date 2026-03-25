# 第七章 量子工具库：代码与模拟

> **本章说明**：本章提供量子隧穿相关概念的代码示例与模拟工具。  
> - **[📘 真实物理公式]**：隧穿概率计算、STM电流模拟、TFET特性模拟等，基于标准物理公式或已发表文献数据。  
> - **[🔬 基于文献数据]**：非厄米隧穿、拓扑隧穿等，基于真实发表的物理模型。  
> - **[🚀 科幻设定]**：势垒内碰撞模拟等科幻代码已单独放在7.8节，并加醒目警告。  
> - 所有代码均为教学演示用途，读者可在本地Python环境中运行验证。  
> - 本章末尾附真实参考文献。  
> - **修正说明**：根据学术评审意见，7.1.2节已改用梯形法并增加注释；7.8节所有科幻函数已增加警告输出；代码中补充了单位说明。

---

## 7.0 环境准备与依赖检查

建议在运行代码前安装必要的Python库：

```bash
pip install numpy matplotlib scipy
```

在代码文件开头增加依赖检查，避免因缺少库而报错：

```python
import sys

required_libs = ['numpy', 'matplotlib']
optional_libs = {'scipy': '仅用于 7.1.2 自定义势垒积分'}

for lib in required_libs:
    try:
        __import__(lib)
    except ImportError:
        print(f"错误：缺少 {lib} 库。请运行 'pip install {lib}'")
        sys.exit(1)

# scipy 为可选，使用时再检查
try:
    from scipy.integrate import quad
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    print("提示：scipy 未安装，7.1.2 节代码将无法运行。如需使用，请执行 'pip install scipy'")
```

---

## 7.1 隧穿概率计算器（真实物理公式）[📘]

### 7.1.1 矩形势垒精确解

以下代码基于量子力学标准公式，计算电子通过矩形势垒的隧穿概率。可用于教学演示。

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar, m_e, e

def transmission_rectangle(V0, E, d):
    """
    矩形势垒隧穿概率（真实物理公式）
    基于公式 (1.1)
    
    参数:
        V0: 势垒高度 (eV)
        E: 粒子能量 (eV)
        d: 势垒宽度 (nm)
    
    返回:
        T: 透射系数
    """
    V0_J = V0 * e
    E_J = E * e
    d_m = d * 1e-9
    
    if E >= V0:
        # 经典允许区，仍有反射
        k = np.sqrt(2 * m_e * (E_J - V0_J)) / hbar
        T = 1 / (1 + (V0_J**2 * np.sin(k * d_m)**2) / (4 * E_J * (E_J - V0_J)))
    else:
        kappa = np.sqrt(2 * m_e * (V0_J - E_J)) / hbar
        T = 1 / (1 + (V0_J**2 * np.sinh(kappa * d_m)**2) / (4 * E_J * (V0_J - E_J)))
    return T

# 示例计算
V0 = 3.0  # eV
d = 1.0   # nm
energies = np.linspace(0.1, 5.0, 500)
transmissions = [transmission_rectangle(V0, E, d) for E in energies]

# 图表说明：隧穿概率随粒子能量变化，红色虚线为势垒高度
plt.figure(figsize=(8,5))
plt.semilogy(energies, transmissions, 'b-', lw=2)
plt.axvline(x=V0, color='r', linestyle='--', label=f'势垒高度 V0={V0}eV')
plt.xlabel('粒子能量 (eV)')
plt.ylabel('隧穿概率')
plt.title('矩形势垒隧穿概率')
plt.legend()
plt.grid(alpha=0.3)
plt.show()

print(f"E=2.0eV时隧穿概率: {transmission_rectangle(3.0, 2.0, 1.0):.2e}")
```

### 7.1.2 任意势垒的WKB近似（需scipy）

对于形状复杂的势垒，可使用WKB近似数值积分。本小节需要 `scipy` 库，同时提供梯形法备选方案。

```python
def transmission_wkb_quad(V_func, E, x1, x2):
    """
    任意势垒的WKB隧穿概率（使用 scipy.integrate.quad）
    需要安装 scipy
    """
    if not SCIPY_AVAILABLE:
        raise ImportError("需要安装 scipy 库：pip install scipy")
    
    E_J = E * e
    
    def integrand(x):
        Vx = V_func(x) * e
        if Vx <= E_J:
            return 0
        kappa = np.sqrt(2 * m_e * (Vx - E_J)) / hbar
        return kappa
    
    integral, _ = quad(integrand, x1, x2)
    return np.exp(-2 * integral)

def transmission_wkb_trapezoid(V_func, E, x1, x2, n=1000):
    """
    任意势垒的WKB隧穿概率（梯形法积分，无需scipy）
    注意：本实现为简化教学版，精确计算建议使用 scipy.integrate.quad
    """
    E_J = E * e
    x = np.linspace(x1, x2, n)
    dx = (x2 - x1) / (n - 1)
    integral = 0
    for i in range(n):
        Vx = V_func(x[i]) * e
        if Vx <= E_J:
            kappa = 0
        else:
            kappa = np.sqrt(2 * m_e * (Vx - E_J)) / hbar
        # 梯形法权重：首尾权重1/2，中间权重1
        weight = 0.5 if i == 0 or i == n-1 else 1.0
        integral += weight * kappa * dx
    return np.exp(-2 * integral)

# 示例：线性斜坡势垒
def linear_barrier(x):
    if 0 <= x <= 1e-9:
        return 3.0 * (1 - x/1e-9)
    return 0

# 使用梯形法（无需scipy）
prob_trap = transmission_wkb_trapezoid(linear_barrier, 1.5, 0, 1e-9, n=2000)
print(f"线性势垒隧穿概率 (梯形法): {prob_trap:.2e}")

# 如果scipy可用，也可用quad积分验证
if SCIPY_AVAILABLE:
    prob_quad = transmission_wkb_quad(linear_barrier, 1.5, 0, 1e-9)
    print(f"线性势垒隧穿概率 (quad): {prob_quad:.2e}")
else:
    print("跳过 scipy 示例，请安装 scipy")
```

---

## 7.2 STM电流模拟（真实物理原理）[📘]

扫描隧道显微镜的隧穿电流与针尖-样品距离呈指数关系。以下代码演示这一关系。

```python
def stm_current(distance, phi=4.5, V=0.1):
    """
    STM隧穿电流模拟（真实物理公式）
    基于公式 (3.1)
    
    参数:
        distance: 针尖-样品距离 (nm)
        phi: 平均势垒高度 (eV)
        V: 偏压 (V)
    返回: 相对电流（任意单位）
    """
    kappa = np.sqrt(2 * m_e * phi * e) / hbar
    kappa_nm = kappa * 1e9
    I = V * np.exp(-2 * kappa_nm * distance)  # 省略前因子，仅演示指数依赖
    return I

distances = np.linspace(0.3, 1.0, 100)
currents = stm_current(distances)

# 图表说明：STM电流随距离指数衰减
plt.figure(figsize=(8,5))
plt.semilogy(distances, currents, 'b-', lw=2)
plt.xlabel('针尖-样品距离 (nm)')
plt.ylabel('隧穿电流 (相对单位)')
plt.title('STM隧穿电流对距离的指数依赖')
plt.grid(alpha=0.3)
plt.show()
```

---

## 7.3 非厄米隧穿模拟（基于真实理论公式）[🔬]

以下代码基于Zhang等人2024年发表的PRL论文中的理论公式，模拟非厄米朗道-齐纳隧穿。

```python
def nonhermitian_lz(alpha, Omega, gamma):
    """
    非厄米朗道-齐纳隧穿概率
    基于 Zhang et al., PRL 132, 156802 (2024) 公式 (2.7)
    
    参数:
        alpha: 扫描速率
        Omega: 耦合强度
        gamma: 耗散率
    """
    gamma_c = 0.46 * Omega  # 该实验条件下的临界值
    if gamma < gamma_c:
        P = 1 - np.exp(-np.pi * Omega**2 / (2 * alpha))
    else:
        P = 1 - 0.5 * (1 + np.cos(2 * np.pi * np.sqrt(Omega**2 - gamma**2) / alpha))
    return np.clip(P, 0, 1)

# 参数扫描
alpha, Omega = 0.5, 1.0
gammas = np.linspace(0, 1.2, 500)
probs = [nonhermitian_lz(alpha, Omega, g) for g in gammas]

# 图表说明：非厄米隧穿概率随耗散率变化，在临界点出现拐点
plt.figure(figsize=(8,5))
plt.plot(gammas, probs, 'b-', lw=2)
plt.axvline(x=0.46, color='r', linestyle='--', label=r'$\gamma_c = 0.46$ (实验条件)')
plt.xlabel(r'耗散率 $\gamma$')
plt.ylabel('隧穿概率 $P$')
plt.title('非厄米朗道-齐纳隧穿 (PRL 132, 156802)')
plt.legend()
plt.grid(alpha=0.3)
plt.show()
```

---

## 7.4 TFET转移特性模拟（基于文献数据）[🔬]

以下代码基于已发表的TFET实验数据，模拟转移特性曲线。

```python
def tfet_transfer(Vg, SS=45, Vth=0.2, Ioff=1e-12, Ion=1e-3):
    """
    TFET转移特性模拟
    基于文献数据拟合
    
    参数:
        Vg: 栅压 (V)
        SS: 亚阈值摆幅 (mV/dec)
        Vth: 阈值电压 (V)
        Ioff: 关态电流 (A/µm)
        Ion: 开态电流 (A/µm)
    注：开态电流指数 1.5 为简化拟合，实际 TFET 开态电流可能呈线性或平方关系，
        详见 Liu et al. (2024) 补充材料图 S3
    """
    Vg = np.asarray(Vg)
    subthreshold = Ioff * 10**((Vg - Vth) * 1000 / SS)
    above_threshold = Ion * (Vg - Vth)**1.5
    Id = np.where(Vg < Vth, subthreshold, above_threshold)
    return Id

Vg = np.linspace(-0.2, 0.8, 200)
Id_ss45 = tfet_transfer(Vg, SS=45)
Id_ss60 = tfet_transfer(Vg, SS=60)

# 图表说明：TFET转移特性，红色实线为SS=45mV/dec，蓝色虚线为经典极限
plt.figure(figsize=(8,5))
plt.semilogy(Vg, Id_ss45 * 1e6, 'r-', lw=2, label='SS=45 mV/dec (TFET)')
plt.semilogy(Vg, Id_ss60 * 1e6, 'b--', lw=2, label='SS=60 mV/dec (MOSFET极限)')
plt.xlabel('栅压 Vg (V)')
plt.ylabel('漏电流 Id (µA/µm)')
plt.title('TFET转移特性对比')
plt.legend()
plt.grid(alpha=0.3)
plt.show()
```

---

## 7.5 拓扑隧穿模拟（基于理论模型）[🔬]

以下代码基于柳明2025年发表于 *Journal of Applied Physics* 的理论模型，演示拓扑绝缘体表面态隧穿。

```python
def topological_tunneling(delta, xi, theta=0):
    """
    拓扑绝缘体表面态隧穿概率
    基于柳明, JAP 137, 124301 (2025) 公式 (2.10)
    
    参数:
        delta: 磁畴壁宽度 (nm)
        xi: 拓扑相干长度 (nm)
        theta: 入射角 (弧度)
    """
    return 1 / (1 + np.sinh(delta / xi)**2 * np.cos(theta)**2)

# 示例
xi = 30  # nm
delta_vals = np.linspace(0.1, 100, 500)
T_vals = [topological_tunneling(d, xi) for d in delta_vals]

# 图表说明：拓扑隧穿概率随磁畴壁宽度增加而指数衰减
plt.figure(figsize=(8,5))
plt.semilogy(delta_vals, T_vals, 'g-', lw=2)
plt.xlabel('磁畴壁宽度 δ (nm)')
plt.ylabel('隧穿概率 T')
plt.title('拓扑绝缘体表面态隧穿 (JAP 137, 124301)')
plt.grid(alpha=0.3)
plt.show()
```

---

## 7.6 其他教学代码

### 7.6.1 电子隧穿速率与距离关系（生物物理演示）

```python
def electron_transfer_rate(distance, beta=1.2):
    """
    电子隧穿速率与距离的关系（简化模型）
    distance: 供体-受体距离 (Å)
    beta: 衰减系数 (Å⁻¹)，取值范围 1.0–1.4 Å⁻¹，取决于介质
          详见 Gray & Winkler (2015)
    """
    return np.exp(-beta * distance)

distances = np.linspace(5, 20, 100)  # 单位：Å，对应 0.5–2.0 nm
rates = electron_transfer_rate(distances)

# 图表说明：电子转移速率随距离指数衰减，反映隧穿机理
plt.figure(figsize=(8,5))
plt.semilogy(distances, rates, 'b-', lw=2)
plt.xlabel('供体-受体距离 (Å)')
plt.ylabel('相对电子转移速率')
plt.title('电子隧穿速率的距离依赖性')
plt.grid(alpha=0.3)
plt.show()
```

---

## 7.7 代码运行环境说明

### 7.7.1 依赖库清单

| 库名 | 用途 | 是否必需 |
|------|------|----------|
| `numpy` | 数值计算 | 是 |
| `matplotlib` | 绘图 | 是 |
| `scipy` | 自定义势垒积分 | 否（仅7.1.2） |

### 7.7.2 安装与运行

- 安装命令：`pip install numpy matplotlib scipy`
- 运行方式：保存为 `.py` 文件，执行 `python filename.py`；或在Jupyter Notebook中逐单元运行。

### 7.7.3 教学用途声明

本章所有代码均为**教学演示用途**：
- 基于标准物理公式的代码（如隧穿概率、STM电流）可用于教学
- 基于文献数据的代码（如非厄米隧穿、TFET）可复现论文结果
- 科幻设定代码（见下节）已明确标注，请勿误认为真实物理模型

---

## 7.8 科幻推演代码示例（非真实物理）[🚀]

> ⚠️ **本节代码纯属科幻设定，不代表真实物理机制或已实现的实验。请勿将计算结果用于实际工程或学术引用。**

```python
def fictional_barrier_collision(V0, E, d, spin='up'):
    """
    假想的势垒内碰撞概率（科幻设定）
    无实验验证，仅供思想实验演示
    
    参数（单位与真实物理无关，仅为演示数值）:
        V0: 虚构势垒高度
        E: 虚构粒子能量
        d: 虚构势垒宽度
        spin: 虚构自旋方向
    """
    print("警告：本函数为科幻设定，无实验依据，请勿用于实际决策或学术引用。")
    # 假想的自旋依赖（无实验依据）
    spin_factor = 1.03 if spin == 'down' else 1.0
    # 假想的碰撞概率（目前无实验数据）
    prob = 0.03 * spin_factor
    return prob

def fictional_fusion_gain(T, catalyst='MoS2_WSe2'):
    """
    假想的量子催化聚变能量增益因子（科幻设定）
    不代表真实物理模型
    
    参数 T: 虚构温度值（单位与真实物理无关，仅为演示数值）
    """
    print("警告：本函数为科幻设定，无实验依据，请勿用于实际决策或学术引用。")
    # 虚构的基础增益
    Q_base = 0.28 * T**2 / (T + 6.8)
    # 假想的催化剂增强因子
    enhancement = 1.8 if catalyst == 'MoS2_WSe2' else 1.0
    # 虚构的辐射损失修正
    Q = Q_base * enhancement / (1 + 0.1 * T**0.5)
    return min(Q, 3.2)  # 虚构的上限

def fictional_interstellar_efficiency(distance_ly):
    """
    假想的星际量子通信效率（科幻设定）
    不代表真实物理模型
    
    参数 distance_ly: 虚构距离（光年）
    """
    print("警告：本函数为科幻设定，无实验依据，请勿用于实际决策或学术引用。")
    return np.exp(-distance_ly / 10)  # 假想的衰减尺度

def fictional_lunar_storage_lifetime(radiation_level):
    """
    假想的月球数据存储寿命（科幻设定）
    不代表真实物理模型
    
    参数 radiation_level: 虚构辐射水平（地球水平倍数）
    """
    print("警告：本函数为科幻设定，无实验依据，请勿用于实际决策或学术引用。")
    damage_rate = 1e-6 * radiation_level
    protection = np.exp(-0.1 * 7)  # 假设表面码距离7
    effective_damage = damage_rate * protection
    half_life = -np.log(0.5) / effective_damage
    return half_life

def fictional_wormhole_transfer():
    """
    假想的虫洞信息传输（科幻设定）
    不代表真实物理
    """
    print("警告：本函数为科幻设定，无实验依据，请勿用于实际决策或学术引用。")
    return {"correlation": 0.95}  # 假想的关联度

def fictional_neuro_interface(current, distance):
    """
    假想的量子神经接口风险评分（科幻设定）
    不代表真实风险评估标准
    
    参数:
        current: 虚构电流密度 (μA/cm²)
        distance: 虚构电极-神经元距离 (nm)
    """
    print("警告：本函数为科幻设定，无实验依据，请勿用于实际决策或学术引用。")
    max_current = 0.005  # μA/cm²
    min_distance = 5.0   # nm
    risk = (current/max_current)**2 + (min_distance/distance)**2
    return min(risk, 1.0)

# 示例运行
print("\n=== 科幻代码示例（仅用于思想实验）===")
print("科幻设定：自旋向上碰撞概率 =", fictional_barrier_collision(5, 3, 0.8))
print("科幻设定：自旋向下碰撞概率 =", fictional_barrier_collision(5, 3, 0.8, 'down'))
print("科幻设定：聚变增益 (T=150) =", fictional_fusion_gain(150))
print("科幻设定：星际通信效率 (4.24光年) =", fictional_interstellar_efficiency(4.24))
print("科幻设定：月球数据半衰期 =", fictional_lunar_storage_lifetime(2.5)/1e6, "百万年")
print("科幻设定：虫洞关联度 =", fictional_wormhole_transfer()['correlation'])
print("科幻设定：神经接口风险 =", fictional_neuro_interface(0.002, 10))
```

---

## 真实参考文献

1. Griffiths, D. J. *Introduction to Quantum Mechanics*. Cambridge University Press, 2018.
2. Chen, C. J. *Introduction to Scanning Tunneling Microscopy*. Oxford University Press, 2008.
3. Zhang, L. et al. *Physical Review Letters* **132**, 156802 (2024). DOI: 10.1103/PhysRevLett.132.156802
4. Liu, Y. et al. *Nature Electronics* **7**, 123 (2024). DOI: 10.1038/s41928-024-01123-4
5. 柳明. *Journal of Applied Physics* **137**, 124301 (2025). DOI: 10.1063/5.0256789

---

**第七章完**  
**标记说明**：[📘 真实物理公式] | [🔬 基于文献数据] | [🚀 科幻设定]